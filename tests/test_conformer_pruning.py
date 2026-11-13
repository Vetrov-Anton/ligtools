import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from ligtools.psiresp import rdutils
from ligtools.psiresp._pruning import prune_conformers, _PruningState, _pruning_maps


def native_embedding(molecule, count, threshold, threads):
    AllChem.EmbedMultipleConfs(
        molecule, numConfs=count, pruneRmsThresh=threshold, useRandomCoords=True,
        clearConfs=False, randomSeed=1, numThreads=threads, ignoreSmoothingFailures=True,
    )


@pytest.mark.parametrize('smiles', [
    'CCO', 'CCCCCCO', 'c1ccccc1', 'CC(=O)[O-]', '[O-][N+](=O)c1ccccc1',
    'CCOP(=O)([O-])OP(=O)([O-])[O-]', 'C[C@H](O)C(=O)[O-]', '[2H]OC',
])
@pytest.mark.parametrize('threads', [1, 4])
def test_embedding_pruning_matches_rdkit(smiles, threads):
    expected = rdutils.rdmol_from_smiles(smiles)
    actual = Chem.Mol(expected)
    native_embedding(expected, 80, 0.3, threads)
    rdutils.generate_conformers(actual, n_conformers=80, rms_tolerance=0.3, n_threads=threads)
    assert [c.GetId() for c in actual.GetConformers()] == [c.GetId() for c in expected.GetConformers()]
    np.testing.assert_array_equal([c.GetPositions() for c in actual.GetConformers()],
                                  [c.GetPositions() for c in expected.GetConformers()])


def test_existing_conformers_and_ids_are_preserved():
    expected = rdutils.rdmol_from_smiles('CCCCC')
    native_embedding(expected, 2, -1, 1)
    expected.GetConformer(0).SetId(7)
    expected.GetConformer(1).SetId(15)
    actual = Chem.Mol(expected)
    old = [c.GetPositions().copy() for c in actual.GetConformers()]
    native_embedding(expected, 80, 0.3, 4)
    rdutils.generate_conformers(actual, n_conformers=80, rms_tolerance=0.3, n_threads=4)
    assert [c.GetId() for c in actual.GetConformers()] == [c.GetId() for c in expected.GetConformers()]
    np.testing.assert_array_equal([c.GetPositions() for c in actual.GetConformers()],
                                  [c.GetPositions() for c in expected.GetConformers()])
    np.testing.assert_array_equal(old, [c.GetPositions() for c in list(actual.GetConformers())[:2]])


def test_rejected_neighbors_do_not_reject_later_candidates():
    molecule = Chem.MolFromSmiles('CCC')
    # A is close to B; B is close to C; A is far from C. Keep A and C.
    for scale in (1.0, 1.04, 1.08):
        conf = Chem.Conformer(3)
        for index, point in enumerate(np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]]) * scale):
            conf.SetAtomPosition(index, point)
        molecule.AddConformer(conf, assignId=True)
    expected = [list(molecule.GetConformers())[i].GetPositions().copy() for i in (0, 2)]
    prune_conformers(molecule, 0, 0.04)
    np.testing.assert_array_equal([c.GetPositions() for c in molecule.GetConformers()], expected)


def test_spectral_filter_never_skips_a_close_pair():
    molecule = rdutils.rdmol_from_smiles('CCOP(=O)([O-])OP(=O)([O-])[O-]')
    native_embedding(molecule, 25, -1, 2)
    coordinates = np.array([c.GetPositions() for c in molecule.GetConformers()])
    state = _PruningState(coordinates, _pruning_maps(molecule), 0.5)
    for candidate in range(len(coordinates)):
        neighbors = state.tree.query_ball_point(state.descriptors[candidate], state.radius)
        for earlier in range(candidate):
            if state.close(candidate, earlier):
                assert earlier in neighbors


def test_pruning_across_batches_matches_native():
    expected = rdutils.rdmol_from_smiles('CCCCC')
    actual = Chem.Mol(expected)
    native_embedding(expected, 300, 0.25, 4)
    rdutils.generate_conformers(actual, n_conformers=300, rms_tolerance=0.25, n_threads=4)
    np.testing.assert_array_equal([c.GetPositions() for c in actual.GetConformers()],
                                  [c.GetPositions() for c in expected.GetConformers()])


def test_threshold_zero_disables_pruning():
    expected = rdutils.rdmol_from_smiles('CCO')
    actual = Chem.Mol(expected)
    native_embedding(expected, 10, 0, 2)
    rdutils.generate_conformers(actual, n_conformers=10, rms_tolerance=0, n_threads=2)
    np.testing.assert_array_equal([c.GetPositions() for c in actual.GetConformers()],
                                  [c.GetPositions() for c in expected.GetConformers()])
