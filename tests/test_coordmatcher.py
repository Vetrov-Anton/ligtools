"""Geometry and file regressions for all CoordMatcher modes."""
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from scipy.spatial.distance import pdist

from ligtools import CoordMatcher

SMILES = 'CCOC(=O)NCCCl'


def make_pair(directory, smiles=SMILES, hydrogens=True, permute=False):
    directory.mkdir(parents=True, exist_ok=True)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMultipleConfs(mol, numConfs=1, randomSeed=80, numThreads=1)
    if not hydrogens:
        mol = Chem.RemoveHs(mol)
    for i, atom in enumerate(mol.GetAtoms()):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f'{atom.GetSymbol()}{i+1}'.rjust(4))
        info.SetResidueName('LIG')
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    conf = mol.GetConformer()
    reference_coords = conf.GetPositions() + [11, -6, 4]
    for i, xyz in enumerate(reference_coords):
        conf.SetAtomPosition(i, xyz)
    ref = directory / 'reference.pdb'
    Chem.MolToPDBFile(mol, str(ref))
    probe = Chem.Mol(mol)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    for i, xyz in enumerate(reference_coords @ rotation + [30, -25, 17]):
        probe.GetConformer().SetAtomPosition(i, xyz)
    if permute:
        probe = Chem.RenumberAtoms(probe, list(reversed(range(probe.GetNumAtoms()))))
    path = directory / 'probe.pdb'
    Chem.MolToPDBFile(probe, str(path))
    return path, ref


@pytest.mark.parametrize('mode', range(5))
@pytest.mark.parametrize('hydrogens', [False, True])
def test_all_modes_geometry_and_outputs(tmp_path, mode, hydrogens):
    probe, ref = make_pair(tmp_path / 'input with spaces', hydrogens=hydrogens)
    before = pmd.load_file(str(probe))
    original_files = [path.read_bytes() for path in (probe, ref)]
    output = tmp_path / f'aligned-{mode}.pdb'
    matcher = CoordMatcher([probe], [ref], [SMILES])
    matcher.align(modes=[mode], n_confs=12, n_cores=2, seed=80,
                  outputs=[output], draw=False)
    assert output.is_file()
    result = matcher.pmd_objs[0]
    assert len(result.atoms) == len(before.atoms)
    assert [a.name for a in result.atoms] == [a.name for a in before.atoms]
    assert np.isfinite(result.coordinates).all()
    if mode in (1, 3):
        np.testing.assert_allclose(pdist(result.coordinates), pdist(before.coordinates), atol=1e-5)
    reference = Chem.MolFromPDBFile(str(ref), removeHs=False)
    aligned = Chem.MolFromPDBFile(str(output), removeHs=False)
    assert aligned is not None
    assert rdShapeHelpers.ShapeTanimotoDist(aligned, reference) < 0.15
    if mode == 0:
        np.testing.assert_allclose(result.coordinates, pmd.load_file(str(ref)).coordinates, atol=1e-8)
    assert original_files == [path.read_bytes() for path in (probe, ref)]
    exported = tmp_path / f'exported-{mode}.pdb'
    matcher.get_pdbs([exported])
    assert len(pmd.load_file(str(exported)).atoms) == len(before.atoms)


@pytest.mark.parametrize('mode', [2, 4])
@pytest.mark.parametrize('smiles', [SMILES, 'c1ccccc1O', 'CC(=O)[O-]', 'N[C@@H](C)C(=O)O'])
def test_chemistry_and_permuted_atom_order(tmp_path, mode, smiles):
    probe, ref = make_pair(tmp_path, smiles, hydrogens=False, permute=True)
    matcher = CoordMatcher([probe], [ref], [smiles])
    matcher.align(modes=[mode], n_confs=12, n_cores=2, seed=80, draw=False)
    result = matcher.pmd_objs[0]
    assert np.isfinite(result.coordinates).all()
    assert len(result.atoms) == len(pmd.load_file(str(probe)).atoms)
    # Atom identities/order must survive mapping back from SMILES.
    expected = pmd.load_file(str(ref))
    ref_by_name = {a.name: xyz for a, xyz in zip(expected.atoms, expected.coordinates)}
    target = np.array([ref_by_name[a.name] for a in result.atoms])
    assert np.sqrt(np.mean(np.sum((result.coordinates - target)**2, axis=1))) < 0.25



def test_batch_defaults_and_repeated_calls(tmp_path):
    probe, ref = make_pair(tmp_path)
    matcher = CoordMatcher([probe, probe], [ref, ref])
    matcher.align(draw=False)
    assert len(matcher.pmd_objs) == 2
    matcher.align(modes=[0, 1], draw=False)
    assert len(matcher.pmd_objs) == 4


@pytest.mark.parametrize('kwargs', [
    {'modes': [9]}, {'modes': [0, 1]}, {'outputs': []},
    {'modes': [2]}, {'modes': [4]},
    {'modes': [2], 'n_confs': 0}, {'modes': [4], 'n_confs': -1},
])
def test_invalid_alignment_arguments(tmp_path, kwargs):
    probe, ref = make_pair(tmp_path)
    with pytest.raises(ValueError):
        CoordMatcher([probe], [ref]).align(draw=False, **kwargs)


def test_incompatible_exact_match(tmp_path):
    probe, _ = make_pair(tmp_path / 'probe', 'CCCl', hydrogens=False)
    _, ref = make_pair(tmp_path / 'ref', 'CCO', hydrogens=False)
    with pytest.raises(ValueError, match='match|substructure'):
        CoordMatcher([probe], [ref]).align(modes=[0], draw=False)


def test_complex_preserves_scaffold_water_ions_and_box(tmp_path):
    probe, ref = make_pair(tmp_path / 'ligand')
    matcher = CoordMatcher([probe], [ref])
    matcher.align(modes=[1], draw=False)
    scaffold = pmd.Structure()
    residues = [('ALA', 6), ('LIG', 6), ('HOH', 8), ('NA', 11), ('SOL', 8), ('K', 19)]
    for i, (name, element) in enumerate(residues):
        scaffold.add_atom(pmd.Atom(name=f'A{i}', atomic_number=element, mass=12), name, i+1)
    scaffold.coordinates = np.arange(18).reshape(6, 3)
    scaffold.box = [30, 31, 32, 90, 90, 90]
    source = tmp_path / 'complex.pdb'
    scaffold.save(str(source))
    output = tmp_path / 'replaced.pdb'
    rebuilt = matcher.get_complex(output=output, scaffold=source, pdb_names=['LIG'], water=['SOL'], ions=['K'])
    assert len(rebuilt.atoms) == 5 + len(matcher.pmd_objs[0].atoms)
    names = [residue.name for residue in rebuilt.residues]
    assert names.count('LIG') == 1
    assert all(names.count(name) == 1 for name in ['ALA', 'HOH', 'NA', 'SOL', 'K'])
    np.testing.assert_allclose(rebuilt.box, scaffold.box)
    np.testing.assert_allclose(rebuilt.coordinates[0], scaffold.coordinates[0])
    for name in ['HOH', 'NA', 'SOL', 'K']:
        got = [a.idx for r in rebuilt.residues if r.name == name for a in r.atoms]
        old = [a.idx for r in scaffold.residues if r.name == name for a in r.atoms]
        np.testing.assert_allclose(rebuilt.coordinates[got], scaffold.coordinates[old])
    assert output.exists()


@pytest.mark.parametrize('mode', [1, 3, 4])
def test_different_ligands(tmp_path, mode):
    probe, _ = make_pair(tmp_path / 'probe', 'CCO', hydrogens=False)
    _, ref = make_pair(tmp_path / 'reference', 'CCCO', hydrogens=False)
    output = tmp_path / 'nested' / 'aligned.pdb'
    matcher = CoordMatcher([probe], [ref], ['CCO'])
    matcher.align(modes=[mode], n_confs=12, n_cores=2, outputs=[output], draw=False)
    result = Chem.MolFromPDBFile(str(output))
    reference = Chem.MolFromPDBFile(str(ref))
    assert rdShapeHelpers.ShapeTanimotoDist(result, reference) < 0.5


def test_exact_substructure(tmp_path):
    probe, _ = make_pair(tmp_path / 'probe', 'CO', hydrogens=False)
    _, ref = make_pair(tmp_path / 'reference', 'CCO', hydrogens=False)
    matcher = CoordMatcher([probe], [ref])
    matcher.align(modes=[0], draw=False)
    np.testing.assert_allclose(matcher.pmd_objs[0].coordinates,
                               pmd.load_file(str(ref)).coordinates[[1, 2]])


def test_default_mode_with_smiles(tmp_path):
    probe, ref = make_pair(tmp_path)
    matcher = CoordMatcher([probe], [ref], [SMILES])
    matcher.align(n_confs=12, n_cores=2, draw=False)
    assert len(matcher.pmd_objs) == 1


@pytest.mark.parametrize('removed_option', ['opt_geom', 'addHs_to_probe', 'addHs_to_ref'])
def test_removed_options_are_not_accepted(removed_option):
    with pytest.raises(TypeError, match=removed_option):
        CoordMatcher().align(**{removed_option: True})


@pytest.mark.parametrize('mode', [2, 4])
def test_embedding_failure_is_explicit(tmp_path, monkeypatch, mode):
    import ligtools.coordmatcher as module

    probe, ref = make_pair(tmp_path)
    monkeypatch.setattr(module.AllChem, 'EmbedMultipleConfs', lambda *a, **k: [])
    with pytest.raises(RuntimeError, match='could not generate'):
        CoordMatcher([probe], [ref], [SMILES]).align(modes=[mode], n_confs=2, draw=False)


@pytest.mark.parametrize('mode', [2, 4])
@pytest.mark.parametrize('smiles', ['not a smiles', 'CCCC'])
def test_invalid_smiles_or_mismatched_graph(tmp_path, mode, smiles):
    probe, ref = make_pair(tmp_path)
    with pytest.raises(ValueError, match='SMILES'):
        CoordMatcher([probe], [ref], [smiles]).align(modes=[mode], n_confs=2, draw=False)


def test_missing_and_empty_files(tmp_path):
    probe, ref = make_pair(tmp_path)
    with pytest.raises(FileNotFoundError):
        CoordMatcher([tmp_path / 'missing.pdb'], [ref]).align(modes=[0], draw=False)
    empty = tmp_path / 'empty.pdb'
    empty.write_text('END\n')
    with pytest.raises(ValueError, match='PDB'):
        CoordMatcher([empty], [ref]).align(modes=[0], draw=False)


def test_pdb_without_bond_orders(tmp_path):
    probe, ref = make_pair(tmp_path)
    # ParmEd's PDB writer drops bond-order information, like Smi2Gmx.get_pdbs.
    untyped = tmp_path / 'from-topology.pdb'
    pmd.load_file(str(probe)).save(str(untyped))
    matcher = CoordMatcher([untyped, untyped], [ref, ref], [SMILES, SMILES])
    matcher.align(modes=[2, 4], n_confs=12, n_cores=2, draw=False)
    assert all(len(s.atoms) == len(pmd.load_file(str(probe)).atoms) for s in matcher.pmd_objs)


def test_mismatched_lists_and_export_count(tmp_path):
    probe, ref = make_pair(tmp_path)
    with pytest.raises(ValueError, match='reference'):
        CoordMatcher([probe], [])
    with pytest.raises(ValueError, match='smiles'):
        CoordMatcher([probe], [ref], [SMILES, SMILES])
    matcher = CoordMatcher([probe], [ref])
    matcher.align(modes=[0], draw=False)
    with pytest.raises(ValueError, match='filenames'):
        matcher.get_pdbs([])
