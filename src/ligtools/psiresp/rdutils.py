
import itertools
import logging
from typing import TYPE_CHECKING, Set, Tuple, Optional, List
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import multiprocessing

from typing_extensions import Literal
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdmolfiles import SmilesParserParams
except ImportError:
    raise ImportError("rdkit is necessary for this functionality. "
                      "Please install it with "
                      "`conda install -c conda-forge rdkit`")


if TYPE_CHECKING:
    import rdkit
    import qcelemental

BONDTYPES = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
}

logger = logging.getLogger(__name__)

OPT_FF_NAMES = Literal["uff", "mmff94"]
OFF_SMILES_ATTRIBUTE = "canonical_isomeric_explicit_hydrogen_mapped_smiles"


def rdmol_from_smiles(smiles, order_by_map_number: bool = False):
    smiles_parser = SmilesParserParams()
    smiles_parser.removeHs = False
    rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles, smiles_parser))
    if order_by_map_number:
        map_numbers = [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()]
        if 0 not in map_numbers:
            new_order = tuple(map(int, np.argsort(map_numbers)))
            rdmol = Chem.RenumberAtoms(rdmol, new_order)
            for atom in rdmol.GetAtoms():
                atom.SetAtomMapNum(0)
    return rdmol


def rdmol_to_smiles(rdmol, mapped=True):
    rdmol = Chem.Mol(rdmol)
    if not any(at.GetAtomMapNum() for at in rdmol.GetAtoms()) and mapped:
        for i, at in enumerate(rdmol.GetAtoms(), 1):
            at.SetAtomMapNum(i)

    return Chem.MolToSmiles(rdmol, allBondsExplicit=True, allHsExplicit=True)


def assign_connectivity_from_proximity(rdmol, scale_factor=0.55):
    table = Chem.GetPeriodicTable()
    conformer = rdmol.GetConformer(0)
    n_atoms = rdmol.GetNumAtoms()
    for i, atom1 in enumerate(rdmol.GetAtoms()):
        xyz1 = np.array(conformer.GetAtomPosition(i))
        rad1 = table.GetRvdw(atom1.GetSymbol())
        for j in range(i + 1, n_atoms):
            atom2 = rdmol.GetAtomWithIdx(j)
            xyz2 = np.array(conformer.GetAtomPosition(j))
            rad2 = table.GetRvdw(atom2.GetSymbol())

            max_distance = (rad1 + rad2) * scale_factor
            distance = np.linalg.norm(xyz1 - xyz2)
            if distance < max_distance:
                rdmol.AddBond(i, j, Chem.BondType.SINGLE)


def n_unpaired_e(atom):
    table = Chem.GetPeriodicTable()
    current = atom.GetTotalValence() - atom.GetFormalCharge()
    valences = list(table.GetValenceList(atom.GetSymbol()))
    return [valence - current for valence in valences]


def assign_bond_order_and_charges_from_electrons(rdmol):
    """
    Modelled after Cedric Bouysset's method for assigning bond order and charges
    as detailed in https://cedric.bouysset.net/blog/2020/07/22/rdkit-converter-part2 .

    In short:
        - if atom has no unpaired electrons, no action
        - if atom has too many paired electrons, turn them into positive charge
        - if atom has unpaired electrons:
             - search for neighbors with unpaired electrons
             - set bond between them to smallest one that satisfies the lowest
               number of electrons in any state in either atom
             - repeat for all neighbors
             - if atom has remaining unpaired electrons, turn them into negative charge
    """

    # sort by # unpaired electrons to deal with unpaired electrons first
    for atom in sorted(rdmol.GetAtoms(), reverse=True, key=n_unpaired_e):
        n_unpaired = n_unpaired_e(atom)
        if len(n_unpaired) == 1:  # single valence state
            n = n_unpaired[0]
            # set excess electrons to positive charge
            if n < 0:
                atom.SetFormalCharge(-n)
                rdmol.UpdatePropertyCache(strict=False)
            if n <= 0:
                continue

        # deal with unpaired electrons
        n_unpaired = np.array(n_unpaired_e(atom))
        # compare with unpaired electrons of neighboring atoms
        # again, sort for the likeliest pairs first
        for neighbor in sorted(
            atom.GetNeighbors(),
            reverse=True,
            key=n_unpaired_e,
        ):
            neighbor_unpaired = np.array(n_unpaired_e(neighbor))
            # set bond between neighbor and self to smallest possible
            # to satisfy valences
            smallest_overlap = min(
                min(n_unpaired[n_unpaired >= 0], default=0),
                min(neighbor_unpaired[neighbor_unpaired >= 0], default=0),
            )
            if smallest_overlap == 0:
                continue
            bond = rdmol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            bond.SetBondType(BONDTYPES[smallest_overlap + 1])
            # update valence states
            rdmol.UpdatePropertyCache(strict=False)
            n_unpaired = np.array(n_unpaired_e(atom))

        # any residual unpaired electrons are turned into negative charge
        n_unpaired = n_unpaired_e(atom)[0]
        if n_unpaired > 0:
            atom.SetFormalCharge(-n_unpaired)
            atom.SetNumRadicalElectrons(0)
            rdmol.UpdatePropertyCache(strict=False)
    Chem.Kekulize(rdmol)


def rdmol_from_qcelemental(qcmol: "qcelemental.models.Molecule",
                           guess_connectivity: bool = True):
    import qcelemental as qcel
    smiles = qcmol.model_dump().get("extras", {}).get(OFF_SMILES_ATTRIBUTE)
    if smiles:
        rdmol = rdmol_from_smiles(smiles, order_by_map_number=True)
        coordinates = qcmol.geometry * qcel.constants.conversion_factor("bohr", "angstrom")
        add_conformer_from_coordinates(rdmol, coordinates=coordinates)
    else:
        rdmol = Chem.RWMol(_rdmol_from_qcelemental(qcmol))
    if guess_connectivity:
        if not qcmol.connectivity:
            assign_connectivity_from_proximity(rdmol)
        rdmol.UpdatePropertyCache(strict=False)
        assign_bond_order_and_charges_from_electrons(rdmol)

    rdmol = Chem.Mol(rdmol)
    Chem.SanitizeMol(rdmol)
    return rdmol


def rdmol_from_psi4(psi4mol):
    molstr = psi4mol.format_molecule_for_mol()
    rdmol = Chem.MolFromMolBlock(molstr, removeHs=False, sanitize=False)
    return rdmol


def get_rdkit_connectivity(rdmol):
    return np.array([
        [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble()]
        for bond in rdmol.GetBonds()
    ])


def rdmol_to_qcelemental(rdmol, multiplicity=1, random_seed=-1):
    import qcelemental as qcel
    rdmol = Chem.AddHs(rdmol)
    if not rdmol.GetNumConformers():
        Chem.rdDistGeom.EmbedMolecule(rdmol, useRandomCoords=True,
                                      randomSeed=random_seed)

    connectivity = get_rdkit_connectivity(rdmol)

    for i, atom in enumerate(rdmol.GetAtoms(), 1):
        atom.SetAtomMapNum(i)
    extras = {}
    extras[OFF_SMILES_ATTRIBUTE] = Chem.MolToSmiles(rdmol,
                                                    allBondsExplicit=True,
                                                    allHsExplicit=True,
                                                    )
    geometry = np.array(rdmol.GetConformer(0).GetPositions())
    schema = dict(symbols=[a.GetSymbol() for a in rdmol.GetAtoms()],
                  geometry=geometry * qcel.constants.conversion_factor("angstrom", "bohr"),
                  connectivity=connectivity if len(connectivity) else None,
                  molecular_charge=Chem.GetFormalCharge(rdmol),
                  molecular_multiplicity=multiplicity,
                  extras=extras
                  )
    qcmol = qcel.models.v2.Molecule.from_data(schema)
    return qcmol


def _rdmol_from_qcelemental(qcmol: "qcelemental.models.Molecule"):
    import qcelemental as qcel
    rwmol = Chem.RWMol()
    for symbol in qcmol.symbols:
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(True)
        rwmol.AddAtom(atom)
    if qcmol.connectivity is not None:
        for i, j, d in qcmol.connectivity:
            if np.isclose(d, 1.5):
                bondtype = Chem.BondType.AROMATIC
            else:
                bondtype = BONDTYPES.get(int(d), Chem.BondType.UNSPECIFIED)
            rwmol.AddBond(i, j, bondtype)
    coordinates = qcmol.geometry * qcel.constants.conversion_factor("bohr", "angstrom")
    add_conformer_from_coordinates(rwmol, coordinates=coordinates)
    return rwmol


def add_conformer_from_coordinates(rdmol: "rdkit.Chem.Mol",
                                   coordinates: np.ndarray):
    """Add conformer to RDKit from coordinate array

    Parameters
    ----------
    rdmol: rdkit.Chem.Mol
        RDKit Molecule
    coordinates: numpy.ndarray of floats
        Coordinates in angstrom
    """
    from rdkit import Geometry
    n_atoms = rdmol.GetNumAtoms()
    if coordinates.shape != (n_atoms, 3):
        raise ValueError("Shape of coordinates must be (n_atoms, 3)")
    conformer = Chem.Conformer(n_atoms)
    for i, xyz in enumerate(coordinates):
        x, y, z = map(float, xyz)
        conformer.SetAtomPosition(i, Geometry.Point3D(x, y, z))
    conformer.SetId(rdmol.GetNumConformers())
    rdmol.AddConformer(conformer)


def generate_conformers(rdmol: "rdkit.Chem.Mol",
                        n_conformers: int = 0,
                        rms_tolerance: float = 1.5,
                        random_seed=1,
                        minimize: Optional[OPT_FF_NAMES] = None,
                        n_threads: int = 1,
                        ):
    """Generate conformers for an RDKit molecule.

    This does not clear existing conformers.

    Parameters
    ----------
    rdmol: rdkit.Chem.Mol
        RDKit Molecule
    n_conformers: int (optional)
        Maximum number of conformers to generate
    rms_tolerance: float (optional)
        RMSD threshold used to prune new conformers
    n_threads: int (optional)
        RDKit threads for embedding and optional minimization
    """
    Chem.SanitizeMol(rdmol)
    n_existing = rdmol.GetNumConformers()
    AllChem.EmbedMultipleConfs(rdmol, numConfs=n_conformers,
                               pruneRmsThresh=-1,
                               useRandomCoords=True,
                               clearConfs=False,
                               randomSeed=random_seed,
                               numThreads=n_threads,
                               ignoreSmoothingFailures=True)
    if rms_tolerance > 0:
        from ._pruning import prune_conformers
        prune_conformers(rdmol, n_existing, rms_tolerance, n_threads)
    if minimize:
        minimize_conformer_geometries(rdmol, forcefield=minimize, n_threads=n_threads)


def minimize_conformer_geometries(rdmol: "rdkit.Chem.Mol",
                                  forcefield: OPT_FF_NAMES = "mmff94",
                                  minimize_max_iter: int = 2000,
                                  n_threads: int = 1):
    """Minimize conformer geometries of an RDKit molecule

    Parameters
    ----------
    minimize_max_iter: int (optional)
        Maximum number of iterations for minimization
    """
    # TODO: make this more customizable

    forcefield = forcefield.lower()
    if forcefield == "uff":
        AllChem.UFFOptimizeMoleculeConfs(rdmol, numThreads=n_threads,
                                         maxIters=minimize_max_iter)
    elif forcefield == "mmff94":
        AllChem.MMFFOptimizeMoleculeConfs(rdmol, numThreads=n_threads,
                                          maxIters=minimize_max_iter)
    else:
        raise ValueError(
            f"Force field {forcefield} not supported. "
            f"Please use one of {OPT_FF_NAMES}"
        )


def compute_mmff_charges(rdmol: "rdkit.Chem.Mol",
                         forcefield: Literal["MMFF94", "MMFF94s"] = "MMFF94",
                         normalize_partial_charges: bool = True,
                         ) -> np.ndarray:
    mps = AllChem.MMFFGetMoleculeProperties(rdmol, forcefield)
    if mps is None:
        molstr = Chem.MolToSmiles(rdmol)
        raise ValueError(f"MMFF charges could not be computed for {molstr}")
    n_atoms = rdmol.GetNumAtoms()
    charges = np.array([mps.GetMMFFPartialCharge(i) for i in range(n_atoms)])

    if normalize_partial_charges:
        total_charge = Chem.GetFormalCharge(rdmol)
        partial_charges = charges.sum()
        offset = (total_charge - partial_charges) / n_atoms
        charges += offset
    return charges


def get_exclusions(rdmol: "rdkit.Chem.Mol") -> Set[Tuple[int, int]]:
    exclusions = set()
    for i, atom in enumerate(rdmol.GetAtoms()):
        partners = [b.GetOtherAtomIdx(i) for b in atom.GetBonds()]
        exclusions |= set(tuple(sorted([i, x])) for x in partners)
        exclusions |= set(itertools.combinations(sorted(partners), 2))
    return exclusions


def compute_distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    dist_sq = np.einsum('ijk,ilk->ijl', coordinates, coordinates)
    diag = np.einsum("ijj->ij", dist_sq)
    a, b = diag.shape
    dist_sq += dist_sq - diag.reshape((a, 1, b)) - diag.reshape((a, b, 1))
    diag[:] = -0.0
    return np.sqrt(-dist_sq)


def compute_electrostatic_energy(rdmol: "rdkit.Chem.Mol",
                                 forcefield: Literal["MMFF94", "MMFF94s"] = "MMFF94",
                                 n_threads: int = 1,
                                 ) -> np.ndarray:
    charges = np.abs(compute_mmff_charges(rdmol)).reshape(-1, 1)
    charge_products = charges @ charges.T

    excl_i, excl_j = zip(*get_exclusions(rdmol))
    charge_products[(excl_i, excl_j)] = 0.0
    charge_products[(excl_j, excl_i)] = 0.0

    conformers = list(rdmol.GetConformers())
    # Bound temporary distance arrays independently of the initial pool size.
    def block_energy(start):
        coordinates = np.array([c.GetPositions() for c in conformers[start:start + 32]])
        distances = compute_distance_matrix(coordinates)
        inverse = np.reciprocal(distances, out=np.zeros_like(distances),
                                where=~np.isclose(distances, 0))
        return 0.5 * (inverse * charge_products).sum(axis=(1, 2))

    starts = range(0, len(conformers), 32)
    if not conformers:
        return np.array([])
    if n_threads > 1 and len(conformers) > 32:
        with ThreadPoolExecutor(max_workers=n_threads) as pool:
            return np.concatenate(list(pool.map(block_energy, starts)))
    return np.concatenate([block_energy(start) for start in starts])


_rms_worker_molecule = None


def _initialize_rms_worker(binary):
    from .._parallel import initialize_worker
    initialize_worker()
    global _rms_worker_molecule
    _rms_worker_molecule = Chem.Mol(binary)


def _rms_block(molecule, reference, ids, n_threads=1):
    # GetBestRMS changes probe coordinates; this API leaves them untouched.
    return np.array([
        0.0 if probe == reference else AllChem.GetBestAlignmentTransform(
            molecule, molecule, prbCid=probe, refCid=reference,
            numThreads=n_threads,
        )[0] for probe in ids
    ])


def _rms_worker_block(task):
    reference, ids = task
    return _rms_block(_rms_worker_molecule, reference, ids)


class _RMSRows:
    """Calculate only the rows needed by greedy selection, without mutating poses."""

    def __init__(self, molecule, ids, n_threads):
        self.molecule = Chem.RemoveHs(molecule)
        self.ids = ids
        self.n_threads = n_threads
        self.pool = None

    def __enter__(self):
        if self.n_threads > 1 and len(self.ids) >= 64:
            self.pool = ProcessPoolExecutor(
                max_workers=self.n_threads, mp_context=multiprocessing.get_context('spawn'),
                initializer=_initialize_rms_worker, initargs=(self.molecule.ToBinary(),),
            )
        return self

    def row(self, index):
        reference = self.ids[index]
        if self.pool is None:
            return _rms_block(self.molecule, reference, self.ids, self.n_threads)
        size = max(1, (len(self.ids) + self.n_threads * 4 - 1) // (self.n_threads * 4))
        tasks = [(reference, self.ids[start:start + size])
                 for start in range(0, len(self.ids), size)]
        return np.concatenate(list(self.pool.map(_rms_worker_block, tasks)))

    def __exit__(self, *exc):
        if self.pool is not None:
            self.pool.shutdown(wait=True, cancel_futures=True)


def compute_heavy_rms(rdmol: "rdkit.Chem.Mol",
                      conformer_ids=None, n_threads: int = 1):
    if conformer_ids is None:
        conformer_ids = [c.GetId() for c in rdmol.GetConformers()]

    rdmol = Chem.RemoveHs(rdmol)
    n_conformers = len(conformer_ids)
    subset = Chem.Mol(rdmol)
    subset.RemoveAllConformers()
    for identifier in conformer_ids:
        subset.AddConformer(rdmol.GetConformer(identifier), assignId=True)
    rms = np.zeros((n_conformers, n_conformers))
    if n_conformers > 1:
        indices = np.tril_indices(n_conformers, -1)
        rms[indices] = AllChem.GetAllConformerBestRMS(subset, numThreads=n_threads)
        rms += rms.T
    return rms


def select_elf_conformer_ids(rdmol: "rdkit.Chem.Mol",
                             energy_window: float = 30,
                             limit: int = 10,
                             rms_tolerance: float = 0.05,
                             n_threads: int = 1,
                             ):
    import qcelemental as qcel
    n_conformers = rdmol.GetNumConformers()
    if n_conformers == 0 or limit <= 0:
        return []

    all_conformer_ids = [c.GetId() for c in rdmol.GetConformers()]

    try:
        energies = compute_electrostatic_energy(rdmol, forcefield="MMFF94", n_threads=n_threads)
    except ValueError as e:
        warnings.warn(str(e))
        return all_conformer_ids[:limit]

    sorting = np.argsort(energies)
    conversion = qcel.constants.conversion_factor(
        "(4 * pi * electric_constant) * kcal",
        "e * e / angstrom",
    )
    window = energy_window * conversion / qcel.constants.get("avogadro constant")
    upper_energy = window + energies[sorting[0]]
    cutoff = np.searchsorted(energies[sorting], upper_energy)
    conformer_ids = [all_conformer_ids[i] for i in sorting[:cutoff]]
    n_max_output = min(limit, len(conformer_ids))
    if not n_max_output:
        return []
    selected_indices = [0]
    if n_max_output > 1:
        with _RMSRows(rdmol, conformer_ids, n_threads) as rows:
            total_rms = np.zeros(len(conformer_ids))
            too_close = np.zeros(len(conformer_ids), dtype=bool)
            for _ in range(n_max_output - 1):
                row = rows.row(selected_indices[-1])
                total_rms += row
                too_close |= row < rms_tolerance
                too_close[selected_indices] = True
                if np.all(too_close):
                    break
                selected_indices.append(int(np.where(too_close, -np.inf, total_rms).argmax()))

    return [conformer_ids[i] for i in selected_indices]


def select_elf_conformer_coordinates(rdmol: "rdkit.Chem.Mol",
                                     energy_window: float = 15,
                                     limit: int = 10,
                                     rms_tolerance: float = 0.05,
                                     n_threads: int = 1):
    ids = select_elf_conformer_ids(rdmol, energy_window=energy_window,
                                   limit=limit, rms_tolerance=rms_tolerance,
                                   n_threads=n_threads)
    return np.array([rdmol.GetConformer(i).GetPositions() for i in ids])


def generate_diverse_conformer_coordinates(molecule,
                                           n_conformer_pool: int = 4000,
                                           energy_window: float = 30,
                                           n_max_conformers: int = 10,
                                           rms_tolerance: float = 0.05,
                                           minimize: Optional[OPT_FF_NAMES] = None,
                                           n_threads: int = 1,
                                           ):
    if not isinstance(molecule, Chem.Mol):
        molecule = rdmol_from_qcelemental(molecule)
    else:
        molecule = Chem.RWMol(molecule)

    generate_conformers(molecule, n_conformers=n_conformer_pool,
                        rms_tolerance=rms_tolerance,
                        minimize=minimize,
                        n_threads=n_threads,
                        )
    return select_elf_conformer_coordinates(molecule,
                                            energy_window=energy_window,
                                            limit=n_max_conformers,
                                            rms_tolerance=rms_tolerance,
                                            n_threads=n_threads)


def get_sp3_ch_indices(rdmol):
    symbols = np.array([at.GetSymbol() for at in rdmol.GetAtoms()])
    bonds = get_rdkit_connectivity(rdmol)
    single_bonds = np.isclose(bonds[:, 2], np.ones_like(bonds[:, 2]))

    groups = {}
    for i in np.where(symbols == "C")[0]:
        contains_index = np.any(bonds[:, :2] == i, axis=1)
        c_bonds = bonds[contains_index & single_bonds][:, :2]
        c_partners = (c_bonds[c_bonds != i]).astype(int)
        if len(c_partners) == 4:
            groups[i] = c_partners[symbols[c_partners] == "H"]
    return groups


def molecule_to_rdkit(molecule):
    if molecule._rdmol is None:
        rdmol = rdmol_from_qcelemental(molecule.qcmol,
                                       guess_connectivity=True)
    else:
        rdmol = Chem.Mol(molecule._rdmol)

    for conformer in molecule.conformers:
        add_conformer_from_coordinates(rdmol, conformer.coordinates)

    charges = molecule.charges
    if charges is not None:
        for charge, atom in zip(charges, rdmol.GetAtoms()):
            atom.SetDoubleProp("PartialCharge", charge)
        Chem.CreateAtomDoublePropertyList(rdmol, "PartialCharge")
    rdmol.UpdatePropertyCache(strict=False)
    return Chem.Mol(rdmol)


def molecule_from_rdkit(rdmol, molecule_cls, random_seed=1, **kwargs):
    rdmol = Chem.AddHs(rdmol)
    qcmol = rdmol_to_qcelemental(rdmol, random_seed=random_seed)
    kwargs = dict(**kwargs)
    kwargs["qcmol"] = qcmol
    obj = molecule_cls(**kwargs)
    obj._rdmol = rdmol

    for conformer in rdmol.GetConformers():
        obj.add_conformer_with_coordinates(np.array(conformer.GetPositions()))
    return obj


def get_connectivity(molecule):
    return get_rdkit_connectivity(molecule._rdmol)


def GetSymmetricAtomIndices(
    mol: Chem.Mol,
    maxMatches: int = 10000,
) -> List[Tuple[int, ...]]:
    """Get atom indices of symmetric atoms

    Returns
    -------
    symmetric_indices: List[Tuple[int, ...]]
        In this list, one item is a sorted tuple of indices,
        where each index indicates an atom that is symmetric
        to all the other indices in the tuple.
    """
    # take care of resonance
    matches = [
        match
        for resMol in Chem.ResonanceMolSupplier(mol)
        for match in resMol.GetSubstructMatches(mol, uniquify=False,
                                                maxMatches=maxMatches)
    ]
    atom_symmetries = set(tuple(sorted(set(x))) for x in zip(*matches))
    return sorted([x for x in atom_symmetries if len(x) > 1])


def get_symmetric_atom_indices(molecule):
    # from psiresp.utils import require_package
    # require_package("rdkit_utilities")

    # from rdkit_utilities import GetSymmetricAtomIndices
    rdmol = molecule.to_rdkit()
    return GetSymmetricAtomIndices(rdmol)
