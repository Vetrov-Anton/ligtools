import numpy as np
import pytest

from ligtools import Smi2Gmx, psiresp
from ligtools.psiresp import rdutils


@pytest.mark.parametrize('pool, expected', [(12, [12, 12]), ([10, 16], [10, 16])])
def test_public_options_reach_rdkit(monkeypatch, pool, expected):
    embedding = rdutils.AllChem.EmbedMultipleConfs
    calls = []

    def record(mol, **kwargs):
        calls.append((kwargs['numConfs'], kwargs['numThreads']))
        return embedding(mol, **kwargs)

    monkeypatch.setattr(rdutils.AllChem, 'EmbedMultipleConfs', record)
    ligands = Smi2Gmx(
        smiles=['CCCC', 'CCCO'], res_names=['BUT', 'PRO'],
        max_conf_n=[2, 3], n_conformer_pool=pool, n_conformer_threads=2,
    )
    for molecule, limit in zip(ligands.mols, [2, 3]):
        molecule.generate_conformers()
        assert 0 < len(molecule.conformers) <= limit
        assert np.isfinite(molecule.conformers[0].coordinates).all()
    assert calls == [(count, 2) for count in expected]


@pytest.mark.parametrize('forcefield, function', [
    ('uff', 'UFFOptimizeMoleculeConfs'), ('mmff94', 'MMFFOptimizeMoleculeConfs'),
])
def test_minimization_respects_thread_limit(monkeypatch, forcefield, function):
    original = getattr(rdutils.AllChem, function)
    threads = []

    def record(mol, **kwargs):
        threads.append(kwargs['numThreads'])
        return original(mol, **kwargs)

    monkeypatch.setattr(rdutils.AllChem, function, record)
    molecule = psiresp.Molecule.from_smiles('CCCC', conformer_generation_options={
        'n_conformer_pool': 10, 'n_threads': 2, 'n_max_conformers': 2,
        'minimize': forcefield,
    })
    molecule.generate_conformers()
    assert threads == [2]


@pytest.mark.parametrize('kwargs', [
    {'n_conformer_pool': 0}, {'n_conformer_pool': -1},
    {'n_conformer_pool': 1.5}, {'n_conformer_pool': True},
    {'n_conformer_pool': []}, {'n_conformer_pool': [10, 20]},
    {'n_conformer_pool': [0]},
    {'n_conformer_threads': 0}, {'n_conformer_threads': -1},
    {'n_conformer_threads': 1.5}, {'n_conformer_threads': True},
])
def test_invalid_options_fail_early(kwargs):
    with pytest.raises(ValueError):
        Smi2Gmx(smiles=['O'], res_names=['HOH'], **kwargs)


def test_defaults_and_separate_qm_limit():
    ligands = Smi2Gmx(smiles=['O'], res_names=['HOH'], max_n_cores=8)
    options = ligands.mols[0].conformer_generation_options
    assert options.n_conformer_pool == 4000
    assert options.n_threads == 1
    assert options.n_max_conformers == 5
