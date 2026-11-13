import numpy as np
import pytest

from ligtools import Smi2Gmx, psiresp
from ligtools.psiresp import rdutils


@pytest.mark.parametrize('pool, expected', [(12, [12, 12]), ([10, 16], [10, 16])])
def test_public_options_reach_rdkit(monkeypatch, pool, expected):
    embedding = rdutils.AllChem.EmbedMultipleConfs
    minimization = rdutils.AllChem.MMFFOptimizeMoleculeConfs
    selection = rdutils.select_elf_conformer_ids
    calls = []
    minimization_threads = []
    selection_threads = []

    def record(mol, **kwargs):
        calls.append((kwargs['numConfs'], kwargs['numThreads']))
        return embedding(mol, **kwargs)

    def record_minimization(mol, **kwargs):
        minimization_threads.append(kwargs['numThreads'])
        properties = rdutils.AllChem.MMFFGetMoleculeProperties(mol)
        before = [rdutils.AllChem.MMFFGetMoleculeForceField(
            mol, properties, confId=conf.GetId()).CalcEnergy() for conf in mol.GetConformers()]
        result = minimization(mol, **kwargs)
        assert all(energy <= initial + 1e-6 for (_, energy), initial in zip(result, before))
        return result

    def record_selection(mol, **kwargs):
        selection_threads.append(kwargs['n_threads'])
        return selection(mol, **kwargs)

    monkeypatch.setattr(rdutils.AllChem, 'EmbedMultipleConfs', record)
    monkeypatch.setattr(rdutils.AllChem, 'MMFFOptimizeMoleculeConfs', record_minimization)
    monkeypatch.setattr(rdutils, 'select_elf_conformer_ids', record_selection)
    ligands = Smi2Gmx(
        smiles=['CCCC', 'CCCO'], res_names=['BUT', 'PRO'],
        max_conf_n=[2, 3], n_conformer_pool=pool, n_conformer_threads=2,
    )
    for molecule, limit in zip(ligands.mols, [2, 3]):
        molecule.generate_conformers()
        assert 0 < len(molecule.conformers) <= limit
        assert np.isfinite(molecule.conformers[0].coordinates).all()
    assert calls == [(count, 2) for count in expected]
    assert minimization_threads == [2, 2]
    assert selection_threads == [2, 2]


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
    {'optimize_geometry': 'false'}, {'optimize_geometry': 1},
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
    assert options.minimize == 'mmff94'
    assert ligands.optimize_geometry is False
    assert ligands.mols[0].optimize_geometry is False


def test_qm_geometry_optimization_can_be_enabled_explicitly():
    ligands = Smi2Gmx(smiles=['O'], res_names=['HOH'], optimize_geometry=True)
    assert ligands.optimize_geometry is True
    assert ligands.mols[0].optimize_geometry is True
