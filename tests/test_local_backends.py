"""Numerical and format regressions for the maintained backends."""
import sys

import numpy as np
import pytest
import qcelemental as qcel

from ligtools import Smi2Gmx, psiresp
from ligtools.acpype import topol
from ligtools.psiresp.qm import QMEnergyOptions
from ligtools.psiresp.qcutils import QCWaveFunction


def test_wavefunction_loading_does_not_infer_chemical_structure(monkeypatch):
    from ligtools.psiresp import rdutils
    qcmol = psiresp.Molecule.from_smiles(
        'CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCOP(=O)([O-])OP(=O)([O-])[O-]'
    ).qcmol.model_copy(update={'extras': {}})
    count = len(qcmol.symbols)
    wavefunction = qcel.models.v2.WavefunctionProperties(
        basis={'name': 'synthetic', 'atom_map': ['s'] * count, 'center_data': {
            's': {'electron_shells': [{'angular_momentum': [0], 'harmonic_type': 'cartesian',
                                     'exponents': [1.0], 'coefficients': [[1.0]]}]},
        }}, restricted=True, scf_orbitals_a=np.eye(count), orbitals_a='scf_orbitals_a',
    )
    result = qcel.models.v2.AtomicResult(
        input_data=qcel.models.v2.AtomicInput(molecule=qcmol, specification={
            'program': 'psi4', 'driver': 'energy', 'model': {'method': 'hf', 'basis': 'synthetic'},
            'protocols': {'wavefunction': 'all'},
        }), molecule=qcmol, wavefunction=wavefunction,
        properties={'calcinfo_nalpha': 1, 'return_energy': -1.0}, return_result=-1.0,
        provenance={'creator': 'test', 'version': '1', 'routine': 'energy'},
    )

    def unexpected(*args, **kwargs):
        raise AssertionError('Reading a wavefunction must not call RDKit graph reconstruction')

    monkeypatch.setattr(rdutils, 'rdmol_from_qcelemental', unexpected)
    loaded = QCWaveFunction.from_atomicresult(result)
    np.testing.assert_array_equal(loaded.qcmol.geometry, qcmol.geometry)
    density = loaded.reconstruct_density()
    np.testing.assert_allclose(density, np.diag([1.0] + [0.0] * (count - 1)))


def test_local_modules_and_no_numpy_patch():
    assert topol.__name__ == 'ligtools.acpype.topol'
    assert psiresp.__name__ == 'ligtools.psiresp'
    assert 'psiresp' not in sys.modules
    assert 'acpype' not in sys.modules
    if np.lib.NumpyVersion(np.__version__) >= '2.4.0':
        assert not hasattr(np, 'in1d')


def test_schema_boundary(tmp_path):
    mol = psiresp.Molecule.from_smiles('O')
    options = QMEnergyOptions()
    path = options.write_input(mol.qcmol, tmp_path)
    payload = qcel.util.deserialize(path.read_bytes(), 'msgpack')
    assert payload['schema_version'] == 1
    assert payload['model']['method'] == options.method
    assert len(payload['molecule']['symbols']) == 3
    # A resumed run may encounter an input file that Psi4 has not executed yet.
    with pytest.raises(SystemExit, match='QM computations'):
        options.manage_external_output([mol.qcmol], tmp_path)
    # Simulate the QCSchema v1 result produced by the external Psi4 executable.
    payload.update(schema_name="qcschema_output", success=True, return_result=-75.0,
                   properties={'return_energy': -75.0},
                   provenance={'creator': 'test', 'version': '1', 'routine': 'test'})
    path.write_bytes(qcel.util.serialize(payload, 'msgpack-ext'))
    result = options.read_output(mol.qcmol, tmp_path)
    assert isinstance(result, qcel.models.v2.AtomicResult)
    assert result.properties.return_energy == -75.0


@pytest.mark.parametrize('two_stages', [False, True])
def test_recovers_known_charges_from_coulomb_potential(two_stages):
    mol = psiresp.Molecule.from_smiles('O', random_seed=1)
    mol.generate_conformers()
    conf = mol.conformers[0]
    conf.add_orientation_with_coordinates(conf.coordinates)
    orientation = conf.orientations[0]
    expected = np.array([-0.8, 0.4, 0.4])
    grid = np.random.default_rng(42).normal(size=(200, 3)) * 3 + [5, 0, 0]
    distances = np.linalg.norm(grid[:, None, :] - orientation.coordinates, axis=-1)
    orientation.grid = grid
    orientation.esp = qcel.constants.conversion_factor('bohr', 'angstrom') / distances @ expected
    job = psiresp.Job(molecules=[mol], resp_options={
        'restrained_fit': False, 'stage_2': two_stages,
    })
    job.compute_charges()
    np.testing.assert_allclose(mol.charges, expected, atol=1e-6)
    assert mol.charges.sum() == pytest.approx(0, abs=1e-8)


def test_model_roundtrip_and_independent_defaults():
    mol = psiresp.Molecule.from_smiles('CC(=O)C(=O)[O-]')
    restored = psiresp.Molecule.model_validate_json(mol.model_dump_json())
    np.testing.assert_allclose(restored.coordinates, mol.coordinates)
    assert restored.get_hash() == mol.get_hash()
    a, b = psiresp.Job(), psiresp.Job()
    a.molecules.append(mol)
    assert b.molecules == []


def test_public_api():
    ligands = Smi2Gmx(smiles=['CC(=O)C(=O)[O-]'], res_names=['PYR'])
    assert ligands.net_charges == [-1]
