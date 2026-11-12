"""Parameterization of small molecules (ligands) for GROMACS/AMBER.

Takes SMILES strings, runs a RESP2 charge fit with psiresp/Psi4, derives GAFF2
parameters with acpype/AmberTools and writes GROMACS `.itp`, `.gro`, `.pdb` and
`.rtp` files, optionally injecting the new residues into an existing force field
directory.
"""

import os
import shutil
import subprocess
from collections import defaultdict
from multiprocessing import get_context
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import parmed as pmd

from . import psiresp
from IPython.display import display
from parmed import gromacs
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

from .acpype import topol

from ._rdkit_utils import simplify_molecule
from ._parallel import initialize_worker

psiresp.grid.GridOptions.vdw_radii = {'Br': 3.95}

IPythonConsole.ipython_useSVG = True
IPythonConsole.molSize = (1000, 1000)

#: Default Psi4 executable; resolved through $PATH when not an absolute path.
DEFAULT_PSI4 = os.environ.get('PSI4_EXE', 'psi4')


def _to_path(path):
    """Normalize a user supplied path into an absolute Path.

    Empty string / None (the "not given" convention used by Smi2Gmx) -> None,
    so callers can keep using `if not path` checks.
    """
    if path is None or path == '':
        return None
    return Path(path).expanduser().resolve()


class Smi2Gmx:
    def __init__(self, smiles=[], res_names = [], max_conf_n = [], chiral = [], net_charges=[], ff_path='', atomtypes_filename='', ffnonbonded_filename = '',
                work_dir = '', max_n_cores = 32, *, n_conformer_pool=4000, n_conformer_threads=1):
        self.smiles = smiles
        self.res_names = res_names
        self.n_mol = len(smiles)


        if not max_conf_n:
            self.max_conf_n = [5 for _ in smiles]
        else:
            self.max_conf_n = max_conf_n

        if isinstance(n_conformer_pool, (list, tuple)):
            pools = list(n_conformer_pool)
            if len(pools) != self.n_mol:
                raise ValueError('n_conformer_pool must have one value per SMILES')
        else:
            # Validate scalars even when no molecules were supplied.
            psiresp.ConformerGenerationOptions(n_conformer_pool=n_conformer_pool)
            pools = [n_conformer_pool] * self.n_mol
        psiresp.ConformerGenerationOptions(n_threads=n_conformer_threads)
        conformer_options = [
            psiresp.ConformerGenerationOptions(
                n_conformer_pool=pool,
                n_threads=n_conformer_threads,
                n_max_conformers=self.max_conf_n[i],
                keep_original_conformer=False,
            )
            for i, pool in enumerate(pools)
        ]
        self.n_conformer_pool = pools
        self.n_conformer_threads = n_conformer_threads
        
        if not chiral:
            self.chiral = [False for _ in smiles]
        else:
            self.chiral = chiral
        if work_dir == '':
            self.work_dir = Path.cwd() / 'tmp_Smi2Gmx'
        else:
            self.work_dir = _to_path(work_dir)
        self.ff_path = _to_path(ff_path)
        self.atomtypes_filename = atomtypes_filename
        self.ffnonbonded_filename = ffnonbonded_filename
       
        self.max_n_cores = max_n_cores
        # load smiles to psiresp
        if smiles:
            self.mols = [
                psiresp.Molecule.from_smiles(
                    mol, optimize_geometry=True, conformer_generation_options=conformer_options[i]
                )
                for i, mol in enumerate(self.smiles)
            ]
            self.rdkit_mols = [Chem.MolFromSmiles(smile) for smile in self.smiles]
            if not net_charges:
                self.net_charges = [mol.charge for mol in self.mols]
            else:
                self.net_charges = net_charges
            
            self._check_input()

        self.job = None
        self.charges = None
        self.equivalent_atoms = None
        self.has_rtps = False

    def _check_input(self):
        errors = []
        if self.smiles:
            if self.res_names == []:
                errors.append('List of names of introduced residues NOT GIVEN. Check variable res_names in Smi2Gmx')
            #if self.ff_path == '':
            #    errors.append('Path to GROMACS forcefield NOT GIVEN. Check variable ff_path in  Smi2Gmx')
            #if self.atomtypes_filename == '':
             #   errors.append('Filename, in which atomtypes are written in your forcefield, NOT GIVEN (e.g. atomtypes.atp). Check variable atomtypes_filename in Smi2Gmx')
        #else: # smiles == []
            #if self.xyz_targets == []:
              #  errors.append('Structure files (pdb, gro, ...) of target molecules for aligning to reference coordinates NOT GIVEN. Check variable xyz_targets in Smi2Gmx')
            #if self.xyz_refs == []:
                #errors.append('Reference structure files (pdb, gro, ...) NOT GIVEN. Check variable xyz_refs in Smi2Gmx')
        output = '\n' + '\n'.join(errors)
        if errors:
            raise Exception(output)
        
    def _check_chirality(self):
        pass

    # --- path helpers -------------------------------------------------------

    @property
    def _acpype_dir(self):
        return self.work_dir / 'acpype'

    @property
    def _temporary_pdbs_dir(self):
        return self.work_dir / 'temporary_pdbs'

    def _temporary_pdb(self, res_name):
        return self._temporary_pdbs_dir / f'{res_name}.pdb'

    def _acpype_file(self, res_name, suffix):
        """Path of a file produced by acpype for `res_name` (suffix: 'itp', 'gro', ...)."""
        return self._acpype_dir / f'{res_name}.acpype' / f'{res_name}_GMX.{suffix}'

    @staticmethod
    def _output_dir(dir, default):
        """Resolve an output directory argument and make sure it exists."""
        out_dir = _to_path(dir) or default
        out_dir.mkdir(parents=True, exist_ok=True)
        return out_dir


    def draw_smiles(self, charge = True, index = True, molsPerRow = 3, subImgSize=(500,500)):
        if self.job:
            molecules = []
            for m in self.job.molecules:
        #Chem.CreateAtomDoublePropertyList(m,'PartialCharge')
                for i,atom in enumerate(m.rdmol.GetAtoms()):
                    if index:
                        atom.SetAtomMapNum(atom.GetIdx())
                    if charge:
                        atom.SetDoubleProp('PartialCharge', m.charges[i])
                        atom.SetProp('atomNote',str("%4.3f" %  m.charges[i]))
                molecules.append(m.rdmol)
            img = Draw.MolsToGridImage(molecules, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=self.res_names)
            display(img)
        else:
            img = Draw.MolsToGridImage(self.rdkit_mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=self.res_names)
            display(img)
    
    def draw_refs(self):
        pass

    def _get_equivalent_atoms(self):

        def get_equivalent_atoms(ranks):
            """obtain equivalent atoms based on their ranks"""

            groups = defaultdict(list)
            for idx, rank in enumerate(ranks):
                groups[rank].append(idx)
            
            return [tuple(sorted(indices)) for _, indices in sorted(groups.items()) 
                    if len(indices) >= 2]



        mol_graphs = [ simplify_molecule(mol, add_hs=True) for mol in self.rdkit_mols ]

        ranks = [ list(Chem.CanonicalRankAtoms(mol_graph, breakTies=False)) for mol_graph in mol_graphs ]

        self.equivalent_atoms = [ get_equivalent_atoms(rank) for rank in ranks]
        
        

    @staticmethod
    def _run_psi4_input(psi4, input_path):
        log_path = input_path.with_suffix('.psi4.log')
        with log_path.open('w') as log:
            result = subprocess.run(
                [psi4, '--qcschema', str(input_path), '--nthread', '1'],
                cwd=input_path.parent, stdout=log, stderr=subprocess.STDOUT,
                check=False,
            )
        if result.returncode:
            raise RuntimeError(
                f'Psi4 exited with code {result.returncode} for {input_path}; log: {log_path}'
            )

    def _run_qm_stage(self, action, psi4, workers):
        try:
            action()
        except psiresp.qm.QMTasksPending as pending:
            paths = pending.paths
            print(f'Running {len(paths)} Psi4 calculations with up to {workers} workers')
            if workers == 1:
                for path in paths:
                    self._run_psi4_input(psi4, path)
            else:
                # Threads only wait for independent external executables.
                with ThreadPoolExecutor(max_workers=min(workers, len(paths))) as pool:
                    futures = [pool.submit(self._run_psi4_input, psi4, path) for path in paths]
                    try:
                        for future in as_completed(futures):
                            future.result()
                    except BaseException:
                        for future in futures:
                            future.cancel()
                        raise
            # Validate QCSchema results as well as the process exit codes.
            try:
                action()
            except psiresp.qm.QMTasksPending as incomplete:
                raise RuntimeError(
                    f'Psi4 finished without producing required results: {incomplete.paths}; '
                    'see the corresponding .psi4.log files'
                ) from incomplete

    def run_psi4(self, psi4_path=DEFAULT_PSI4, parallel=True):
        """Run optimization and single points using Python's external-process scheduler."""
        executable = shutil.which(str(Path(psi4_path).expanduser()))
        if executable is None:
            raise FileNotFoundError(f'Psi4 executable not found: {psi4_path}')
        executable = str(Path(executable).resolve())
        self.job = psiresp.Job(
            molecules=self.mols, working_directory=str(self.work_dir),
            n_processes=self.max_n_cores,
        )
        workers = self.max_n_cores if parallel else 1
        self.job.generate_conformers()
        self._run_qm_stage(self.job.optimize_geometries, executable, workers)
        self.job.generate_orientations()
        self._run_qm_stage(self.job.compute_orientation_energies, executable, workers)

    def run_resp(self, constraints=[]):
        if constraints == []:
            constraints = ['auto' for _ in self.res_names]
        self._get_equivalent_atoms()
        constr = psiresp.ChargeConstraintOptions()

        for i, mol in enumerate(self.mols):
            constr.add_charge_sum_constraint_for_molecule(mol, charge=self.net_charges[i])
            if constraints[i] == 'no':
                continue
            elif constraints[i] == 'auto':
                for eq_group in self.equivalent_atoms[i]:
                    constr.add_charge_equivalence_constraint_for_molecule(mol, indices=eq_group)
            else: # custom
                for eq_group in constraints[i]:
                    constr.add_charge_equivalence_constraint_for_molecule(mol, indices=eq_group)

        self.job.charge_constraints = constr

        self.charges = self.job.run()

    def _parallel_acpype(self, i_mol):
        i = i_mol[0]
        mol = i_mol[1]

        
        pdb = self._temporary_pdb(self.res_names[i])

        mol.to_mdanalysis().atoms.write(str(pdb), reindex = False)
        structure = pmd.load_file(str(pdb))
        structure.residues[0].name = self.res_names[i]
        structure.save(str(pdb), overwrite=True, renumber=False )
        new_mol = topol.ACTopol(str(pdb),chargeType='user',chargeVal=self.net_charges[i],basename=self.res_names[i],outTopol='gmx',chiral=self.chiral[i])
        new_mol.createACTopol()
        new_mol.createMolTopol()


    def run_acpype(self, parallel = True):
        previous_dir = Path.cwd()
        if parallel:
            self._temporary_pdbs_dir.mkdir(parents=True, exist_ok=True)
            self._acpype_dir.mkdir(parents=True, exist_ok=True)

            os.chdir(self._acpype_dir)
            try:
                tasks = [ (i,self.mols[i]) for i in range(len(self.res_names))]
                print(tasks)

                n_proc =  len(self.res_names) if len(self.res_names) < self.max_n_cores else self.max_n_cores
                if tasks:
                    with ProcessPoolExecutor(
                        max_workers=n_proc, mp_context=get_context('spawn'), initializer=initialize_worker,
                    ) as pool:
                        list(pool.map(self._parallel_acpype, tasks))
            finally:
                os.chdir(previous_dir)
        else:
            self._temporary_pdbs_dir.mkdir(parents=True, exist_ok=True)
            for i, mol in enumerate(self.job.molecules):
                pdb = self._temporary_pdb(self.res_names[i])
                mol.to_mdanalysis().atoms.write(str(pdb), reindex = False)
                structure = pmd.load_file(str(pdb))
                structure.residues[0].name = self.res_names[i]
                structure.save(str(pdb), overwrite=True, renumber=False )
            self._acpype_dir.mkdir(parents=True, exist_ok=True)
            os.chdir(self._acpype_dir)
            try:
                for i, res_name in enumerate(self.res_names):
                    mol = topol.ACTopol(str(self._temporary_pdb(res_name)),chargeType='user',chargeVal=self.net_charges[i],basename=res_name,outTopol='gmx',chiral=self.chiral[i])
                    mol.createACTopol()
                    mol.createMolTopol()
            finally:
                os.chdir(previous_dir)
            shutil.rmtree(self._temporary_pdbs_dir, ignore_errors=True)


    def get_itps(self, dir = ''):
        print('Saving itps')
        out_dir = self._output_dir(dir, self.work_dir / 'itps')
        for i, name in enumerate(self.res_names):
            top = gromacs.GromacsTopologyFile(str(self._acpype_file(name, 'itp')))
            for j, atom in enumerate(top.atoms):
                atom.charge = self.charges[i][j]
            top.write(str(out_dir / f'{name}_GMX.itp'),combine = 'all', itp=True)
       # for name in self.res_names:
           # shutil.copy(self._acpype_file(name, 'itp'), out_dir / f'{name}_GMX.itp')

    def get_gros(self, dir = ''):
        print('Saving gros')
        out_dir = self._output_dir(dir, self.work_dir / 'gros')
        for name in self.res_names:
            shutil.copy(self._acpype_file(name, 'gro'), out_dir / f'{name}_GMX.gro')
            #structure = pmd.load_file(str(out_dir / f'{name}_GMX.gro'))
           # structure.residues[0].name = name
           #structure.save(str(out_dir / f'{name}_GMX.gro'), )

    def get_pdbs(self, dir = ''):
        print('Saving pdbs')
        out_dir = self._output_dir(dir, self.work_dir / 'pdbs')
        for name in self.res_names:
            structure = pmd.load_file(str(self._acpype_file(name, 'gro')))
            structure.save(str(out_dir / f'{name}.pdb'), overwrite = True, renumber=False)



    def get_rtps(self,dir=''):
        def make_rtp_from_gmx(g_parm,charges,name, c5=[],c3=[],shift = 0,):
            rtptext = []
            rtptext.append('[ bondedtypes ]\n')
            rtptext.append('; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih\n')
            rtptext.append('     1       1          9          4        1         3      1     0\n')
            rtptext.append(f'[ {name} ]\n')
            rtptext.append('[ atoms ]\n')
            for i,a in enumerate(g_parm.atoms):
                rtptext.append("%5s%5s%10.6f%4d\n"% (a.name, a.atom_type, charges[i], a.idx +shift))
            rtptext.append('[ bonds ]\n')
            for b in g_parm.bonds:
                rtptext.append("%5s%5s%8.4f%13.2f\n" % (b.atom1.name, b.atom2.name, 0.1*b.type.req, b.type.k*2*4.184*100))
            rtptext.append('[ angles ]\n')
            for a in g_parm.angles:
            #rtptext += "%5s%5s%5s %8.4f %12.1f\n" % (a.atom1.name, a.atom2.name,a.atom3.name, a.type.theteq, a.type.k*1000)
                rtptext.append("%5s%5s%5s %8.4f %14.3f\n" % (a.atom1.name, a.atom2.name,a.atom3.name, 
                                                    a.type.theteq, a.type.k*2*4.184))
            rtptext.append('[ dihedrals ]\n')
            for d in g_parm.dihedrals:
                if d.funct == 9:
                    rtptext.append("%5s%5s%5s%5s %8.3f %8.4f %5i\n" % (d.atom1.name,d.atom2.name,d.atom3.name,d.atom4.name,
                                                                d.type[0].phase,d.type[0].phi_k*4.184 ,d.type[0].per))
            rtptext.append('[ impropers ]\n')
            for d in g_parm.dihedrals:
                if d.funct == 4:
                    rtptext.append("%5s%5s%5s%5s %8.3f %8.4f %5i\n" % (d.atom1.name,d.atom2.name,d.atom3.name,d.atom4.name,
                                                                d.type.phase,d.type.phi_k,d.type.per ))
            out = ''.join(rtptext)    
            
            return out
        
        print('Saving rtps')
        out_dir = self._output_dir(dir, self.work_dir / 'rtps')

        for i, name in enumerate(self.res_names):
            top = gromacs.GromacsTopologyFile(str(self._acpype_file(name, 'itp')))
            rtp  = make_rtp_from_gmx(top,charges=self.job.molecules[i].charges,name=name)
            with open(out_dir / f'{name}.rtp' , 'w') as f:
                f.write(rtp)

        self.has_rtps = True
        self.rtps_dir = out_dir

    def get_ff_with_rtps(self,ff_path = '', atomtypes_filename = 'atomtypes.atp', ffnonbonded_filename = 'ffnonbonded.itp'):
        print('Intruducing rtps into forcefield')
        # an explicitly given ff_path wins over the one stored on the instance
        ff_path = _to_path(ff_path) or self.ff_path
        if ff_path is None:
            raise Exception('Path to your forcefield NOT GIVEN')
        # `/` keeps absolute filenames absolute, so both names and full paths work
        atomtypes_path = ff_path / atomtypes_filename
        ffnonbonded_path = ff_path / ffnonbonded_filename

        atomtypes_from_atp = set()
        with open(atomtypes_path, 'r') as f:
            for line in f:
                if ';' not in line.strip().split()[0]:
                    atomtypes_from_atp.add(line.strip().split()[0])
        atomtypes_from_ffnonbond = set ()
        with open(ffnonbonded_path, 'r') as f:
            for line in f:
                if line.strip() == '[ atomtypes ]':
                    break
            for line in f:
                if '[' in line.strip():
                    break
                first_field = line.strip().split()[0]
                if ';' in first_field:
                    continue
                else:
                    atomtypes_from_ffnonbond.add(first_field)

        old_atomtypes = atomtypes_from_ffnonbond & atomtypes_from_atp ## maybe union?
        new_atomtypes = {}
        for i, res_name in enumerate(self.res_names):
            itp = pmd.load_file(str(self._acpype_file(res_name, 'itp')))
            for atom in itp.atoms:
                if atom.type not in old_atomtypes:
                    d = {}
                    d['atomic_number'] = atom.atomic_number
                    d['mass'] = atom.mass
                    d['charge'] = float(0)
                    d['ptype'] = 'A'
                    d['sigma'] = atom.sigma/10
                    d['epsilon'] = atom.epsilon*4.184
                    new_atomtypes[atom.type] = d
        
        with open(atomtypes_path, 'a') as f:
            f.write('; Smi2Gmx\n')
            for type in new_atomtypes.keys():
                f.write(f' {type}          {new_atomtypes[type]["mass"]:>10.5f} \n')
        begin_txt = []
        atomtypes_txt = []
        end_txt = []
        with open(ffnonbonded_path, 'r') as f:
            for line in f:
                if line.strip().startswith('[ atomtypes ]'):
                    atomtypes_txt.append(line)
                    break
                begin_txt.append(line)
            for line in f:
                if '[' in line:
                    end_txt.append(line)
                    break
                atomtypes_txt.append(line)
            for line in f:
                end_txt.append(line)

        with open(ffnonbonded_path, 'w') as f:
            atomtypes_txt.append('; Smi2Gmx\n')
            for atomtype in new_atomtypes.keys():
                d = new_atomtypes[atomtype]
                atomtypes_txt.append(f"{atomtype:<8}{d['atomic_number']:>8}{d['mass']:>10.4f}{d['charge']:>10.4f}{d['ptype']:>6}{d['sigma']:>14.5e}{d['epsilon']:>14.5e}\n")

            txt = begin_txt + atomtypes_txt + end_txt
            out = ''.join(txt)
            f.write(out)

        if not self.has_rtps:
            self.get_rtps()

        for i, name in enumerate(self.res_names):
            shutil.copy(self.rtps_dir / f'{name}.rtp', ff_path / f'{name}.rtp')
