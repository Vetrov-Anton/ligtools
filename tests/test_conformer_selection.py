import os
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from rdkit import Chem

from ligtools.psiresp import rdutils


@pytest.fixture
def pool():
    molecule = rdutils.rdmol_from_smiles('CCCCCCO')
    rdutils.generate_conformers(molecule, n_conformers=80, rms_tolerance=-1, n_threads=2)
    return molecule


def reference_selection(molecule, limit, tolerance):
    energies = rdutils.compute_electrostatic_energy(molecule)
    ids = np.argsort(energies).tolist()
    rms = rdutils.compute_heavy_rms(molecule, ids)
    selected = [0]
    for _ in range(min(limit, len(ids)) - 1):
        rows = rms[selected]
        close = (rows < tolerance).any(axis=0)
        close[selected] = True
        if close.all():
            break
        selected.append(int(np.where(close, -np.inf, rows.sum(axis=0)).argmax()))
    return [ids[index] for index in selected]


@pytest.mark.parametrize('threads', [1, 4])
def test_lazy_selection_matches_full_matrix_without_mutating_coordinates(pool, threads):
    before = np.array([c.GetPositions() for c in pool.GetConformers()])
    expected = reference_selection(pool, 5, 0.05)
    environment = dict(os.environ)
    actual = rdutils.select_elf_conformer_ids(
        pool, limit=5, energy_window=1e9, rms_tolerance=0.05, n_threads=threads,
    )
    assert actual == expected
    np.testing.assert_array_equal(before, [c.GetPositions() for c in pool.GetConformers()])
    assert dict(os.environ) == environment


@pytest.mark.parametrize('threads', [1, 4])
def test_energy_blocks_match_original_formula(pool, threads):
    positions = np.array([c.GetPositions() for c in pool.GetConformers()])
    distances = rdutils.compute_distance_matrix(positions)
    inverse = np.reciprocal(distances, out=np.zeros_like(distances), where=~np.isclose(distances, 0))
    charges = np.abs(rdutils.compute_mmff_charges(pool))
    products = charges[:, None] * charges
    i, j = zip(*rdutils.get_exclusions(pool))
    products[i, j] = products[j, i] = 0
    expected = 0.5 * (inverse * products).sum(axis=(1, 2))
    np.testing.assert_allclose(rdutils.compute_electrostatic_energy(pool, n_threads=threads), expected, atol=1e-12)


def test_rms_subset_ids_and_native_parallel_result(pool):
    ids = [70, 2, 55, 9]
    heavy = Chem.RemoveHs(pool)
    expected = np.zeros((4, 4))
    for i in range(4):
        for j in range(i):
            expected[i, j] = expected[j, i] = rdutils.AllChem.GetBestRMS(heavy, heavy, ids[i], ids[j])
    np.testing.assert_allclose(rdutils.compute_heavy_rms(pool, ids, n_threads=4), expected, atol=1e-6)


def test_only_selected_rows_are_computed(pool, monkeypatch):
    original = rdutils._RMSRows.row
    calls = []

    def record(self, index):
        calls.append(index)
        return original(self, index)

    def unexpected(*args, **kwargs):
        raise AssertionError('Full pairwise RMS matrix must not be built during selection')

    monkeypatch.setattr(rdutils._RMSRows, 'row', record)
    monkeypatch.setattr(rdutils, 'compute_heavy_rms', unexpected)
    ids = rdutils.select_elf_conformer_ids(pool, limit=5, energy_window=1e9)
    assert len(ids) == 5
    assert len(calls) == 4
    calls.clear()
    assert len(rdutils.select_elf_conformer_ids(pool, limit=1)) == 1
    assert calls == []


def test_empty_pool_and_zero_limit():
    mol = rdutils.rdmol_from_smiles('CCO')
    assert rdutils.select_elf_conformer_ids(mol) == []
    assert rdutils.select_elf_conformer_ids(mol, limit=0) == []


@pytest.mark.skipif(not all(importlib.util.find_spec(m) for m in ['nbclient', 'ipykernel']),
                    reason='requires notebook test dependencies')
def test_parallel_selection_in_jupyter(tmp_path):
    import nbclient
    import nbformat
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpecManager

    spec = tmp_path / 'kernels' / 'selection-test'
    spec.mkdir(parents=True)
    (spec / 'kernel.json').write_text(json.dumps({
        'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name': 'Selection test', 'language': 'python',
    }))
    code = '''import os
import numpy as np
from ligtools.psiresp import rdutils
assert get_ipython() is not None
environment = dict(os.environ)
mol = rdutils.rdmol_from_smiles('CCCCCCO')
rdutils.generate_conformers(mol, n_conformers=80, rms_tolerance=-1, n_threads=4)
before = np.array([c.GetPositions() for c in mol.GetConformers()])
serial = rdutils.select_elf_conformer_ids(mol, energy_window=1e9, n_threads=1, limit=5)
parallel = rdutils.select_elf_conformer_ids(mol, energy_window=1e9, n_threads=4, limit=5)
assert serial == parallel
np.testing.assert_array_equal(before, [c.GetPositions() for c in mol.GetConformers()])
assert dict(os.environ) == environment
from rdkit import Chem
from rdkit.Chem import AllChem
expected = rdutils.rdmol_from_smiles('CCCCC')
actual = Chem.Mol(expected)
AllChem.EmbedMultipleConfs(expected, numConfs=160, pruneRmsThresh=0.3,
    useRandomCoords=True, clearConfs=False, randomSeed=1, numThreads=4,
    ignoreSmoothingFailures=True)
rdutils.generate_conformers(actual, n_conformers=160, rms_tolerance=0.3, n_threads=4)
np.testing.assert_array_equal([c.GetPositions() for c in actual.GetConformers()],
    [c.GetPositions() for c in expected.GetConformers()])
assert dict(os.environ) == environment
print('PASS: parallel selection in Jupyter')
'''
    notebook = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(code)])
    manager = KernelManager(kernel_name='selection-test', kernel_spec_manager=KernelSpecManager(
        kernel_dirs=[str(tmp_path / 'kernels')],
    ))
    env = dict(os.environ)
    env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1] / 'src')
    try:
        nbclient.NotebookClient(notebook, km=manager, timeout=60).execute(cwd=str(tmp_path), env=env)
    finally:
        nbformat.write(notebook, tmp_path / 'executed.ipynb')
    assert any('PASS: parallel selection in Jupyter' in output.get('text', '')
               for output in notebook.cells[0].outputs)
