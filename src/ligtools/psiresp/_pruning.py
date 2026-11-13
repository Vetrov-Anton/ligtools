"""Ordered embedding-time RMSD pruning with bounded parallel batches.

Match construction and the strict SSD cutoff follow RDKit Embedder.cpp's
getMolSelfMatches/_isConfFarFromRest (BSD-3-Clause, RDKit contributors).
"""
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Numerics import rdAlignment
from scipy.spatial import cKDTree


def _pruning_maps(molecule):
    heavy = Chem.RemoveHs(molecule, sanitize=False)
    query = Chem.RWMol(heavy)
    # Use the same terminal O/N symmetry query as RDKit's embedding pruner.
    atom = 'O,N;D1'
    pattern = Chem.MolFromSmarts(
        f'[{atom};$([{atom}]-[*]=[{atom}]),$([{atom}]=[*]-[{atom}])]~[*]'
    )
    bond = Chem.BondFromSmarts('-,=')
    if getattr(AllChem.EmbedParameters(), 'symmetrizeConjugatedTerminalGroupsForPruning', False):
        for terminal, center in heavy.GetSubstructMatches(pattern):
            query.GetAtomWithIdx(terminal).SetFormalCharge(0)
            query.ReplaceBond(query.GetBondBetweenAtoms(terminal, center).GetIdx(), bond)
    stripped = molecule.GetSubstructMatch(query)
    matches = heavy.GetSubstructMatches(query, uniquify=False, maxMatches=1000)
    if not stripped or not matches:
        raise ValueError('Could not construct embedding RMSD atom mappings')
    return np.asarray(stripped)[np.asarray(matches)]


class _PruningState:
    def __init__(self, coordinates, maps, threshold):
        self.coordinates = coordinates
        self.maps = maps
        # rdAlignment's NumPy binding expects C-contiguous point arrays.
        self.reference = np.ascontiguousarray(coordinates[:, maps[0]])
        centered = self.reference - self.reference.mean(axis=1, keepdims=True)
        # Singular values are invariant to rotations and atom permutations.
        # Their distance / sqrt(n) is a lower bound on optimal RMSD
        # (Hoffman-Wielandt). Only provably distant pairs are skipped.
        self.descriptors = np.linalg.svd(centered, compute_uv=False) / np.sqrt(len(maps[0]))
        self.tree = cKDTree(self.descriptors)
        self.radius = threshold + 1e-7
        self.ssd_cutoff = len(maps[0]) * threshold * threshold

    def close(self, candidate, earlier):
        reference = self.reference[candidate]
        for mapping in self.maps:
            ssd, _ = rdAlignment.GetAlignmentTransform(reference, self.coordinates[earlier, mapping])
            if ssd < self.ssd_cutoff:
                return True
        return False

    def block(self, task):
        start, stop, batch_start, accepted = task
        results = []
        for candidate in range(start, stop):
            neighbors = self.tree.query_ball_point(self.descriptors[candidate], self.radius)
            prior = sorted(index for index in neighbors if index < candidate)
            # A match to an already accepted conformer is conclusive.
            blocked = any(accepted[index] and self.close(candidate, index)
                          for index in prior if index < batch_start)
            local = () if blocked else tuple(
                index for index in prior if index >= batch_start and self.close(candidate, index)
            )
            results.append((candidate, blocked, local))
        return results


_state = None


def _initialize_pruning_worker(coordinates, maps, threshold):
    from .._parallel import initialize_worker
    initialize_worker()
    global _state
    _state = _PruningState(coordinates, maps, threshold)


def _pruning_block(task):
    return _state.block(task)


def prune_conformers(molecule, n_existing, threshold, n_threads=1):
    """Preserve old conformers and apply RDKit's ordered rule to newly added ones.

    Batches share an immutable snapshot of the accepted prefix. Workers check
    that prefix and all earlier candidates within the batch. The parent resolves
    the batch in order, ignoring matches to candidates that were rejected.
    Neither a full distance matrix nor a quadratic global adjacency list is kept.
    """
    conformers = list(molecule.GetConformers())
    count = len(conformers)
    if threshold <= 0 or count <= n_existing:
        return
    coordinates = np.array([conf.GetPositions() for conf in conformers])
    maps = _pruning_maps(molecule)
    accepted = np.zeros(count, dtype=bool)
    accepted[:n_existing] = True
    pool = None
    try:
        if n_threads > 1 and count - n_existing >= 64:
            pool = ProcessPoolExecutor(
                max_workers=n_threads, mp_context=multiprocessing.get_context('spawn'),
                initializer=_initialize_pruning_worker,
                initargs=(coordinates, maps, threshold),
            )
            state = None
        else:
            state = _PruningState(coordinates, maps, threshold)
        for batch in range(n_existing, count, 128):
            end = min(count, batch + 128)
            tasks = [(start, min(start + 16, end), batch, accepted.copy())
                     for start in range(batch, end, 16)]
            blocks = pool.map(_pruning_block, tasks) if pool else map(state.block, tasks)
            for block in blocks:
                for candidate, blocked, local in block:
                    accepted[candidate] = not blocked and not any(accepted[index] for index in local)
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
    # RDKit assigns IDs only when accepting a conformer. Recreate those IDs,
    # including when existing conformers have non-contiguous IDs.
    retained = [Chem.Conformer(conf) for index, conf in enumerate(conformers) if accepted[index]]
    next_id = max((conf.GetId() for conf in retained[:n_existing]), default=-1) + 1
    for conf in retained[n_existing:]:
        conf.SetId(next_id)
        next_id += 1
    molecule.RemoveAllConformers()
    for conf in retained:
        molecule.AddConformer(conf, assignId=False)
