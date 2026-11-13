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
from numbers import Integral

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

    def __init__(self, xyz_probes=None, xyz_refs=None, smiles=None):
        self.smiles = list(smiles or [])
        self.xyz_probes = [_to_path(p) for p in (xyz_probes or [])]
        self.xyz_refs = [_to_path(p) for p in (xyz_refs or [])]
        if len(self.xyz_probes) != len(self.xyz_refs):
            raise ValueError('xyz_refs must have one reference per probe')
        if self.smiles and len(self.smiles) != len(self.xyz_probes):
            raise ValueError('smiles must have one entry per probe')
        self.pmd_objs = []
        self._initial_pmd_objs = None

    @staticmethod
    def _load_pair(probe, ref):
        """Read PDB atoms in their original order without adding or removing atoms."""
        molecules = []
        for path in (probe, ref):
            if path is None or not Path(path).is_file():
                raise FileNotFoundError(f'PDB file not found: {path}')
            mol = Chem.MolFromPDBFile(str(path), removeHs=False)
            if mol is None or not mol.GetNumAtoms():
                raise ValueError(f'Could not read a molecule from PDB: {path}')
            molecules.append(mol)
        return tuple(molecules)

    @staticmethod
    def _structure_with_coords(probe, coordinates):
        structure = pmd.load_file(str(probe))
        coordinates = np.asarray(coordinates, dtype=float)
        n_atoms = len(structure.atoms)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3 or len(coordinates) != n_atoms:
            raise ValueError('Alignment coordinates do not match probe atoms')
        if not np.isfinite(coordinates).all():
            raise ValueError('Alignment produced non-finite coordinates')
        structure.coordinates = coordinates
        return structure

    @staticmethod
    def _save_structure(structure, output):
        output = _to_path(output)
        if output is None:
            raise ValueError('An output filename is required')
        output.parent.mkdir(parents=True, exist_ok=True)
        structure.save(str(output), overwrite=True)

    @staticmethod
    def _mcs_atom_map(probe, reference):
        heavy = []
        for mol in (probe, reference):
            mol = Chem.Mol(mol)
            for atom in mol.GetAtoms():
                atom.SetIntProp('_alignment_original_index', atom.GetIdx())
            heavy.append(Chem.RemoveHs(mol))
        # Fit the common heavy-atom scaffold, preserving all original atoms in output.
        result = rdFMCS.FindMCS(heavy, timeout=30)
        if result.canceled:
            raise ValueError('Common-substructure search timed out')
        if not result.numAtoms:
            raise ValueError('Probe and reference have no common substructure')
        core = Chem.MolFromSmarts(result.smartsString)
        prb = heavy[0].GetSubstructMatch(core)
        ref = heavy[1].GetSubstructMatch(core)
        if not prb or not ref:
            raise ValueError('Could not match the common substructure')
        return [
            (heavy[0].GetAtomWithIdx(i).GetIntProp('_alignment_original_index'),
             heavy[1].GetAtomWithIdx(j).GetIntProp('_alignment_original_index'))
            for i, j in zip(prb, ref)
        ]

    @staticmethod
    def _generate_conformers(smile, n_confs, n_cores, seed):
        if not isinstance(smile, str) or not smile.strip():
            raise ValueError('A valid SMILES is required for conformer search')
        mol = Chem.MolFromSmiles(smile)
        if mol is None or not mol.GetNumAtoms():
            raise ValueError(f'Invalid SMILES: {smile}')
        # Internal SMILES working molecule only; input PDB atoms are untouched.
        mol = Chem.AddHs(mol)
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=seed, numThreads=n_cores)
        if not len(ids):
            raise RuntimeError('RDKit could not generate any conformers')
        return mol

    def _exact_align(self, probe, ref):
        """Copy matched reference coordinates, preserving the probe atom set."""
        mol_probe, mol_ref = self._load_pair(probe, ref)
        match = mol_ref.GetSubstructMatch(mol_probe)
        if len(match) != mol_probe.GetNumAtoms():
            raise ValueError('Exact alignment requires the probe to match a reference substructure')
        coordinates = mol_ref.GetConformer().GetPositions()[list(match)]
        return self._structure_with_coords(probe, coordinates)

    def _align_one_conf(self, probe, ref):
        """Rigid-body align the probe onto the reference over their maximum common substructure."""
        mol_probe, mol_ref = self._load_pair(probe, ref)

        atom_map = self._mcs_atom_map(mol_probe, mol_ref)

        rmsd = AllChem.AlignMol(mol_probe, mol_ref, atomMap=atom_map)
        print(f'RMSD {rmsd}')

        conf = mol_probe.GetConformer()
        new_coord = conf.GetPositions()

        return self._structure_with_coords(probe, new_coord)

    def _map_probe_to_smiles(self, mol_smile, mol_probe, draw=False):
        """Map original probe atoms, preferring chemistry over a topology-only fallback."""
        graph_smile = simplify_molecule(mol_smile, add_hs=False)
        graph_probe = simplify_molecule(mol_probe, add_hs=False)
        if draw:
            display(Draw.MolToImage(mol_smile, size=(500, 500)))
            display(Draw.MolToImage(mol_probe, size=(300, 300)))
        if mol_smile.GetNumHeavyAtoms() != mol_probe.GetNumHeavyAtoms():
            raise ValueError('SMILES must describe the same heavy-atom graph as the probe')
        # Preserve chemical identities when PDB bond orders are available.
        # Topological matching is a fallback for PDBs lacking bond orders.
        probe2smile = mol_smile.GetSubstructMatch(mol_probe)
        if not probe2smile:
            probe2smile = graph_smile.GetSubstructMatch(graph_probe)
        if len(probe2smile) != mol_probe.GetNumAtoms():
            raise ValueError('SMILES does not match the probe atom graph')
        print(f'probe -> smiles atom map: {probe2smile}')
        return probe2smile

    @staticmethod
    def _write_coords(probe, best_conf_coords, probe2smile, output):
        """Write the selected conformer coordinates back in probe atom order."""
        n_atoms = len(pmd.load_file(str(probe)).atoms)
        if len(probe2smile) < n_atoms:
            raise ValueError('SMILES atom map does not cover every probe atom')
        coordinates = best_conf_coords[list(probe2smile)]
        structure = CoordMatcher._structure_with_coords(probe, coordinates)
        if output:
            CoordMatcher._save_structure(structure, output)
        return structure

    def _align_confs(self, probe, ref, smile, n_confs=10000,
                     seed=119, n_cores=16,
                     output=None, draw=False):
        """Generate conformers from SMILES and keep the one with the lowest MCS RMSD."""
        mol_smile = self._generate_conformers(smile, n_confs, n_cores, seed)
        mol_probe, mol_ref = self._load_pair(probe, ref)

        probe2smile = self._map_probe_to_smiles(mol_smile, mol_probe, draw=draw)

        atom_map = self._mcs_atom_map(mol_smile, mol_ref)

        best_rmsd = float('inf')
        best_conf_id = -1
        best_conf_coords = None

        for conf_id in [conf.GetId() for conf in mol_smile.GetConformers()]:
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

    def _align_one_shape(self, probe, ref, opt_param=0.5, output=None):
        """Align the probe onto the reference by maximizing shape overlap."""
        mol_probe, mol_ref = self._load_pair(probe, ref)
        shape_score, color_score = rdShapeAlign.AlignMol(mol_ref, mol_probe, opt_param=opt_param)
        print(f'shape {shape_score}, color {color_score}')

        conf = mol_probe.GetConformer()
        new_coord = conf.GetPositions()

        structure = self._structure_with_coords(probe, new_coord)
        if output:
            self._save_structure(structure, output)
        return structure

    def _align_shapes(self, probe, ref, smile, n_confs=1000,
                      opt_param=0.5, seed=119,
                      n_cores=16, output=None, draw=False):
        """Generate conformers from SMILES and keep the best shape+color overlap."""
        mol_smile = self._generate_conformers(smile, n_confs, n_cores, seed)
        mol_probe, mol_ref = self._load_pair(probe, ref)

        probe2smile = self._map_probe_to_smiles(mol_smile, mol_probe, draw=draw)

        best_metric = float('-inf')
        best_conf_id = -1
        best_conf_coords = None

        for conf_id in [conf.GetId() for conf in mol_smile.GetConformers()]:
            shape_score, color_score = rdShapeAlign.AlignMol(
                mol_ref, mol_smile, opt_param=opt_param, probeConfId=conf_id, refConfId=0)

            metric = shape_score + color_score

            if metric > best_metric:
                best_metric = metric
                best_conf_id = conf_id
                best_conf_coords = mol_smile.GetConformer(best_conf_id).GetPositions()
        print(f'Best metric {best_metric} for conformer {best_conf_id}')

        return self._write_coords(probe, best_conf_coords, probe2smile, output)

    def align(self, modes=None, n_confs=1000, n_cores=16, opt_param=0.5,
              seed=80, outputs=None, draw=False):
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

        draw : bool, optional
            Show RDKit depictions when True. Disabled by default.

        The aligned structures are appended to `self.pmd_objs`.
        """
        n_probes = len(self.xyz_probes)
        modes = list(modes) if modes else [SHAPES if self.smiles else ONE_SHAPE] * n_probes
        if len(modes) != n_probes:
            raise ValueError('modes must have one entry per probe')
        if any(isinstance(mode, bool) or not isinstance(mode, Integral) or mode not in range(5) for mode in modes):
            raise ValueError('Unknown alignment mode: expected integers 0 through 4')
        for name, value, minimum in [('n_confs', n_confs, 1), ('n_cores', n_cores, 0)]:
            if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}')
        if not np.isfinite(opt_param) or not 0 <= opt_param <= 1:
            raise ValueError('opt_param must be between 0 and 1')
        if any(mode in (CONFS, SHAPES) and (not self.smiles or not self.smiles[i])
               for i, mode in enumerate(modes)):
            raise ValueError('SMILES is required for conformer-search modes 2 and 4')
        outputs = [None] * n_probes if outputs is None else list(outputs)
        if len(outputs) != n_probes:
            raise ValueError('outputs must have one entry per probe')

        results = []
        for i, mode in enumerate(modes):
            pair = dict(probe=self.xyz_probes[i], ref=self.xyz_refs[i])
            if mode == EXACT:
                result = self._exact_align(**pair)
            elif mode == ONE_CONF:
                result = self._align_one_conf(**pair)
            elif mode == ONE_SHAPE:
                result = self._align_one_shape(**pair, opt_param=opt_param)
            else:
                search = self._align_confs if mode == CONFS else self._align_shapes
                options = dict(smile=self.smiles[i], n_confs=int(n_confs), n_cores=int(n_cores),
                               seed=seed, draw=draw)
                if mode == SHAPES:
                    options['opt_param'] = opt_param
                result = search(**pair, **options)
            if outputs[i] is not None:
                self._save_structure(result, outputs[i])
            results.append(result)
        self.pmd_objs.extend(results)

    def get_pmd(self) -> list[pmd.Structure]:
        """Return aligned ParmEd structures, or the original probes before alignment.

        Original probes are loaded and cached on first access. The returned
        list is separate, but its Structure objects are not copied. Calling
        this method before align() does not add originals to alignment results.
        """
        if self.pmd_objs:
            return list(self.pmd_objs)
        if self._initial_pmd_objs is None:
            self._initial_pmd_objs = [pmd.load_file(str(path)) for path in self.xyz_probes]
        return list(self._initial_pmd_objs)

    def get_pdbs(self, filenames):
        """Write every aligned structure to the corresponding file name."""
        if len(filenames) != len(self.pmd_objs):
            raise ValueError('filenames must have one path per aligned structure')
        for pmd_obj, filename in zip(self.pmd_objs, filenames):
            self._save_structure(pmd_obj, filename)

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
        self._save_structure(xyz_scaffold, output)
        return xyz_scaffold
