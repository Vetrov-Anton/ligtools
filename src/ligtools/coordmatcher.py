"""Place a ligand into the coordinate frame of a reference structure.

Given a probe molecule (the parameterized ligand, typically straight out of
:class:`ligtools.Smi2Gmx`) and a reference (a ligand sitting in a binding site,
a crystal pose, ...), CoordMatcher rewrites the probe coordinates so that it
occupies the reference position. Five strategies are available, from an exact
atom-by-atom copy to a conformer search scored by shape overlap; see
:meth:`CoordMatcher.align`.

The result can then be merged back into the original complex, keeping water and
ions, with :meth:`CoordMatcher.get_complex`.
"""

from pathlib import Path

import numpy as np
import parmed as pmd
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS, rdShapeAlign
from rdkit.Chem.Draw import IPythonConsole

from ._rdkit_utils import simplify_molecule

IPythonConsole.ipython_useSVG = True
IPythonConsole.molSize = (1000, 1000)

#: Alignment strategies accepted by :meth:`CoordMatcher.align`.
EXACT = 0
ONE_CONF = 1
CONFS = 2
ONE_SHAPE = 3
SHAPES = 4


def _to_path(path):
    """Normalize a user supplied path into an absolute Path ('' / None -> None)."""
    if path is None or path == '':
        return None
    return Path(path).expanduser().resolve()


class CoordMatcher:
    """Align probe structures onto reference structures.

    Parameters
    ----------
    xyz_probes : list of path-like
        Structures to move (PDB files).
    xyz_refs : list of path-like
        Reference structures, one per probe, that define the target position.
    smiles : list of str
        SMILES of the probes. Only needed for the conformer-searching modes
        (`CONFS` and `SHAPES`), which rebuild the molecule from SMILES to
        generate conformers.
    """

    def __init__(self, xyz_probes=[], xyz_refs=[], smiles=[]):
        self.smiles = smiles
        self.xyz_probes = [_to_path(p) for p in xyz_probes]
        self.xyz_refs = [_to_path(p) for p in xyz_refs]
        self.pmd_objs = []

    @staticmethod
    def _load_pair(probe, ref, addHs_to_probe, addHs_to_ref):
        """Read probe and reference PDBs into RDKit molecules."""
        mol_probe = Chem.MolFromPDBFile(str(probe), removeHs=addHs_to_probe)
        mol_ref = Chem.MolFromPDBFile(str(ref), removeHs=addHs_to_ref)
        return mol_probe, mol_ref

    def _exact_align(self, probe, ref, addHs_to_probe=False, addHs_to_ref=False):
        """Copy reference coordinates atom by atom (probe must be a substructure of ref)."""
        mol_probe, mol_ref = self._load_pair(probe, ref, addHs_to_probe, addHs_to_ref)

        match = mol_ref.GetSubstructMatch(mol_probe)

        pmd_probe = pmd.load_file(str(probe))
        pmd_ref = pmd.load_file(str(ref))
        new_coord = []
        for i, atom in enumerate(pmd_probe.atoms):
            old_idx = match[i]
            old_xyz = pmd_ref.coordinates[old_idx]
            new_coord.append([old_xyz[0], old_xyz[1], old_xyz[2]])

        pmd_probe.coordinates = np.array(new_coord)
        return pmd_probe

    def _align_one_conf(self, probe, ref, addHs_to_probe=False, addHs_to_ref=False):
        """Rigid-body align the probe onto the reference over their maximum common substructure."""
        mol_probe, mol_ref = self._load_pair(probe, ref, addHs_to_probe, addHs_to_ref)

        mcs_result = rdFMCS.FindMCS([mol_probe, mol_ref])
        core = Chem.MolFromSmarts(mcs_result.smartsString)

        match1 = mol_ref.GetSubstructMatch(core)
        match2 = mol_probe.GetSubstructMatch(core)

        atom_map = list(zip(match2, match1))

        rmsd = AllChem.AlignMol(mol_probe, mol_ref, atomMap=atom_map)
        print(f'RMSD {rmsd}')

        conf = mol_probe.GetConformer()
        new_coord = conf.GetPositions()

        pmd_probe = pmd.load_file(str(probe))
        pmd_probe.coordinates = new_coord
        return pmd_probe

    def _map_probe_to_smiles(self, mol_smile, mol_probe, draw=True):
        """Map probe atom order onto the SMILES-built molecule, ignoring bond orders."""
        graph_smile = simplify_molecule(mol_smile, add_hs=False)
        graph_probe = simplify_molecule(mol_probe, add_hs=False)
        if draw:
            display(Draw.MolToImage(mol_smile, size=(500, 500)))
            display(Draw.MolToImage(mol_probe, size=(300, 300)))
        probe2smile = graph_smile.GetSubstructMatch(graph_probe)
        print(f'probe -> smiles atom map: {probe2smile}')
        return probe2smile

    @staticmethod
    def _write_coords(probe, best_conf_coords, probe2smile, output):
        """Write the selected conformer coordinates back in probe atom order."""
        pmd_probe = pmd.load_file(str(probe))
        new_coords = []
        for i, atom in enumerate(pmd_probe.atoms):
            smile_idx = probe2smile[i]
            xyz_new = [best_conf_coords[smile_idx][0],
                       best_conf_coords[smile_idx][1],
                       best_conf_coords[smile_idx][2]]
            new_coords.append(xyz_new)
        pmd_probe.coordinates = np.array(new_coords)
        if output:
            pmd_probe.save(str(output), overwrite=True)
        return pmd_probe

    def _align_confs(self, probe, ref, smile, n_confs=10000, addHs_to_probe=False,
                     addHs_to_ref=False, seed=119, opt_geom=False, n_cores=16,
                     output=None, draw=True):
        """Generate conformers from SMILES and keep the one with the lowest MCS RMSD."""
        mol_smile = Chem.AddHs(Chem.MolFromSmiles(smile))
        mol_probe, mol_ref = self._load_pair(probe, ref, addHs_to_probe, addHs_to_ref)

        probe2smile = self._map_probe_to_smiles(mol_smile, mol_probe, draw=draw)

        AllChem.EmbedMultipleConfs(mol_smile, numConfs=n_confs, randomSeed=seed, numThreads=n_cores)

        if opt_geom:
            AllChem.MMFFOptimizeMolecule(mol_smile, maxIters=1000)

        mcs_result = rdFMCS.FindMCS([mol_smile, mol_ref])
        core = Chem.MolFromSmarts(mcs_result.smartsString)

        match_smile = mol_smile.GetSubstructMatch(core)
        match_ref = mol_ref.GetSubstructMatch(core)

        atom_map = list(zip(match_smile, match_ref))

        best_rmsd = float('inf')
        best_conf_id = -1
        best_conf_coords = None

        for conf_id in range(mol_smile.GetNumConformers()):
            rmsd = AllChem.AlignMol(mol_smile, mol_ref,
                                    prbCid=conf_id,
                                    refCid=-1,
                                    atomMap=atom_map)

            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_conf_id = conf_id
                best_conf_coords = mol_smile.GetConformer(best_conf_id).GetPositions()
        print(f'Best rmsd {best_rmsd} for conformer {best_conf_id}')

        return self._write_coords(probe, best_conf_coords, probe2smile, output)

    def _align_one_shape(self, probe, ref, addHs_to_probe=False, addHs_to_ref=False,
                         opt_param=0.5, output=None):
        """Align the probe onto the reference by maximizing shape overlap."""
        mol_probe, mol_ref = self._load_pair(probe, ref, addHs_to_probe, addHs_to_ref)
        shape_score, color_score = rdShapeAlign.AlignMol(mol_ref, mol_probe, opt_param=opt_param)
        print(f'shape {shape_score}, color {color_score}')

        conf = mol_probe.GetConformer()
        new_coord = conf.GetPositions()

        pmd_probe = pmd.load_file(str(probe))
        pmd_probe.coordinates = new_coord
        if output:
            pmd_probe.save(str(output), overwrite=True)
        return pmd_probe

    def _align_shapes(self, probe, ref, smile, n_confs=1000, addHs_to_probe=False,
                      addHs_to_ref=False, opt_param=0.5, seed=119, opt_geom=False,
                      n_cores=16, output=None, draw=True):
        """Generate conformers from SMILES and keep the best shape+color overlap."""
        mol_smile = Chem.AddHs(Chem.MolFromSmiles(smile))
        mol_probe, mol_ref = self._load_pair(probe, ref, addHs_to_probe, addHs_to_ref)

        probe2smile = self._map_probe_to_smiles(mol_smile, mol_probe, draw=draw)

        AllChem.EmbedMultipleConfs(mol_smile, numConfs=n_confs, randomSeed=seed, numThreads=n_cores)

        if opt_geom:
            AllChem.MMFFOptimizeMolecule(mol_smile, maxIters=1000)

        best_metric = float('-inf')
        best_conf_id = -1
        best_conf_coords = None

        for conf_id in range(mol_smile.GetNumConformers()):
            shape_score, color_score = rdShapeAlign.AlignMol(
                mol_ref, mol_smile, opt_param=opt_param, probeConfId=conf_id, refConfId=0)

            metric = shape_score + color_score

            if metric > best_metric:
                best_metric = metric
                best_conf_id = conf_id
                best_conf_coords = mol_smile.GetConformer(best_conf_id).GetPositions()
        print(f'Best metric {best_metric} for conformer {best_conf_id}')

        return self._write_coords(probe, best_conf_coords, probe2smile, output)

    def align(self, modes=[], n_confs=1000, n_cores=16, opt_param=0.5, opt_geom=False,
              seed=80, addHs_to_probe=[], addHs_to_ref=[], outputs=None, draw=True):
        """Align every probe onto its reference and collect the results.

        Parameters
        ----------
        modes : list of int
            One strategy per probe. Defaults to `SHAPES` (4) when SMILES were
            given, `ONE_SHAPE` (3) otherwise.

            ==== ============================================================
            mode strategy
            ==== ============================================================
            0    `EXACT`     - copy reference coordinates atom by atom
            1    `ONE_CONF`  - rigid align over the maximum common substructure
            2    `CONFS`     - conformer search, best MCS RMSD
            3    `ONE_SHAPE` - rigid align by shape overlap
            4    `SHAPES`    - conformer search, best shape + color overlap
            ==== ============================================================
        outputs : list of path-like, optional
            Where to write each aligned probe. Nothing is written when omitted;
            the structures are kept in `self.pmd_objs` either way.

        The aligned structures are appended to `self.pmd_objs`.
        """
        if not self.smiles:
            if not modes:
                modes = [ONE_SHAPE for _ in self.xyz_probes]
        else:
            if not modes:
                modes = [SHAPES for _ in self.xyz_probes]

        if not addHs_to_probe:
            addHs_to_probe = [False for _ in self.xyz_probes]

        if not addHs_to_ref:
            addHs_to_ref = [False for _ in self.xyz_probes]

        if outputs is None:
            outputs = [None for _ in self.xyz_probes]
        else:
            outputs = [_to_path(o) for o in outputs]

        for i, mode in enumerate(modes):
            if mode == EXACT:
                self.pmd_objs.append(self._exact_align(
                    probe=self.xyz_probes[i], ref=self.xyz_refs[i],
                    addHs_to_probe=addHs_to_probe[i], addHs_to_ref=addHs_to_ref[i]))

            elif mode == ONE_CONF:
                self.pmd_objs.append(self._align_one_conf(
                    probe=self.xyz_probes[i], ref=self.xyz_refs[i],
                    addHs_to_probe=addHs_to_probe[i], addHs_to_ref=addHs_to_ref[i]))

            elif mode == CONFS:
                self.pmd_objs.append(self._align_confs(
                    probe=self.xyz_probes[i], ref=self.xyz_refs[i], smile=self.smiles[i],
                    n_confs=n_confs, addHs_to_probe=addHs_to_probe[i],
                    addHs_to_ref=addHs_to_ref[i], seed=seed, opt_geom=opt_geom,
                    n_cores=n_cores, output=outputs[i], draw=draw))

            elif mode == ONE_SHAPE:
                self.pmd_objs.append(self._align_one_shape(
                    probe=self.xyz_probes[i], ref=self.xyz_refs[i],
                    addHs_to_probe=addHs_to_probe[i], addHs_to_ref=addHs_to_ref[i],
                    opt_param=opt_param, output=outputs[i]))

            elif mode == SHAPES:
                self.pmd_objs.append(self._align_shapes(
                    probe=self.xyz_probes[i], ref=self.xyz_refs[i], smile=self.smiles[i],
                    n_confs=n_confs, addHs_to_probe=addHs_to_probe[i],
                    addHs_to_ref=addHs_to_ref[i], seed=seed, opt_geom=opt_geom,
                    n_cores=n_cores, opt_param=opt_param, output=outputs[i], draw=draw))

            else:
                raise ValueError(f'Unknown alignment mode {mode} for probe {self.xyz_probes[i]}')

    def get_pdbs(self, filenames):
        """Write every aligned structure to the corresponding file name."""
        for i, pmd_obj in enumerate(self.pmd_objs):
            pmd_obj.save(str(_to_path(filenames[i])), overwrite=True)

    def get_complex(self, output='new_start.pdb', scaffold='', pdb_names=[], water=[], ions=[]):
        """Rebuild a complex: scaffold + aligned ligands + ions + water.

        The old ligand copies named in `pdb_names` are dropped from the
        scaffold and replaced by the aligned structures in `self.pmd_objs`.
        """
        water = water + ['HOH']
        ions = ions + ['NA', 'CL', 'CA', 'MG']
        xyz = pmd.load_file(str(_to_path(scaffold)))

        # One comma-separated Amber mask per group; joining the groups
        # separately would concatenate the last residue of one group with the
        # first of the next ("HOHNA,CL").
        scaffold_mask = '!(:' + ','.join(water + ions + list(pdb_names)) + ')'
        water_mask = '(:' + ','.join(water) + ')'
        ions_mask = '(:' + ','.join(ions) + ')'

        xyz_scaffold = xyz[scaffold_mask]
        xyz_ions = xyz[ions_mask]
        xyz_water = xyz[water_mask]

        for pmd_obj in self.pmd_objs:
            xyz_scaffold += pmd_obj
        xyz_scaffold += xyz_ions
        xyz_scaffold += xyz_water
        xyz_scaffold.save(str(_to_path(output)), overwrite=True)
        return xyz_scaffold
