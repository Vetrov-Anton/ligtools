import os
from pathlib import Path

import numpy as np
import pytest
import qcelemental as qcel

from ligtools import Smi2Gmx, psiresp
from ligtools.psiresp.base import Model
from ligtools.psiresp.qm import QMGeometryOptimizationOptions, QMEnergyOptions, QMTasksPending


def distorted_water():
    return qcel.models.v2.Molecule(
        symbols=['O', 'H', 'H'], geometry=[[0, 0, 0], [0, 0, 2.8], [2.5, 0, -0.8]],
        connectivity=[(0, 1, 1), (0, 2, 1)], fix_com=True, fix_orientation=True,
    )


def test_default_geometry_tolerance_and_explicit_overrides(tmp_path):
    options = QMGeometryOptimizationOptions()
    path = options.write_input(distorted_water(), tmp_path)
    data = qcel.util.deserialize(path.read_bytes(), 'msgpack')
    keywords = data['specification']['specification']['keywords']
    names = ('max_force_g_convergence', 'rms_force_g_convergence',
             'max_disp_g_convergence', 'rms_disp_g_convergence')
    assert all(keywords[name] == 1e-3 for name in names)
    strict = QMGeometryOptimizationOptions(convergence_tolerance=None)
    assert not any(name in strict.generate_keywords() for name in names)
    assert strict.get_hash() != options.get_hash()
    custom = QMGeometryOptimizationOptions(keywords={'MAX_FORCE_G_CONVERGENCE': 2e-4})
    assert custom.generate_keywords()['MAX_FORCE_G_CONVERGENCE'] == 2e-4
    assert 'max_force_g_convergence' not in custom.generate_keywords()
    assert custom.generate_keywords()['rms_force_g_convergence'] == 1e-3


@pytest.mark.parametrize('value', [0, -1e-3, float('nan'), float('inf')])
def test_invalid_geometry_tolerance(value):
    with pytest.raises(ValueError):
        QMGeometryOptimizationOptions(convergence_tolerance=value)


def test_old_gradient_cache_is_not_used(tmp_path):
    molecule = distorted_water()
    options = QMGeometryOptimizationOptions()
    current = options.get_job_file_for_molecule(molecule, tmp_path, make_directory=True)
    old = current.with_name(current.name.replace(options.get_hash(), Model.get_hash(options)))
    assert old != current
    old.write_bytes(qcel.util.serialize({'success': True, 'driver': 'gradient'}, 'msgpack-ext'))
    with pytest.raises(QMTasksPending) as pending:
        options.manage_external_output([molecule], tmp_path)
    assert pending.value.paths == (current,)
    data = qcel.util.deserialize(current.read_bytes(), 'msgpack')
    assert data['schema_name'] == 'qcschema_optimization_input'
    assert data['specification']['specification']['model'] == {'method': 'hf', 'basis': '6-31g*'}
    assert current.with_suffix('.opt.py').exists()
    assert '--qcschema' not in options.get_run_file(tmp_path).read_text()


def test_single_gradient_is_not_a_converged_geometry():
    molecule = distorted_water()
    result = qcel.models.v2.AtomicResult(
        input_data=qcel.models.v2.AtomicInput(molecule=molecule, specification={
            'program': 'psi4', 'driver': 'gradient', 'model': {'method': 'hf', 'basis': '6-31g*'},
        }), molecule=molecule, properties={}, return_result=np.ones((3, 3)),
        provenance={'creator': 'test', 'version': '1', 'routine': 'gradient'},
    )
    with pytest.raises(ValueError, match='single gradient'):
        QMGeometryOptimizationOptions()._postprocess_result(result)


def test_resp_uses_prepared_orientations_without_restarting_qm(monkeypatch):
    ligands = Smi2Gmx(smiles=['O'], res_names=['WAT'])
    mol = ligands.mols[0]
    mol.generate_conformers()
    conf = mol.conformers[0]
    conf.add_orientation_with_coordinates(conf.coordinates)
    orientation = conf.orientations[0]
    expected = np.array([-0.8, 0.4, 0.4])
    grid = np.random.default_rng(42).normal(size=(300, 3)) * 3 + [5, 0, 0]
    distances = np.linalg.norm(grid[:, None, :] - orientation.coordinates, axis=-1)
    orientation.grid = grid
    orientation.esp = qcel.constants.conversion_factor('bohr', 'angstrom') / distances @ expected
    ligands.job = psiresp.Job(molecules=[mol], resp_options={'restrained_fit': False})

    def unexpected(*args, **kwargs):
        raise AssertionError('RESP fitting must not regenerate QM inputs or orientations')

    monkeypatch.setattr(psiresp.Job, 'generate_orientations', unexpected)
    ligands.run_resp()
    assert conf.orientations[0] is orientation
    np.testing.assert_allclose(ligands.charges[0], expected, atol=1e-5)


@pytest.mark.skipif(not os.environ.get('LIGTOOLS_PSI4'), reason='requires external Psi4 executable')
@pytest.mark.parametrize('max_iter, converges', [(200, True), (1, False)])
def test_real_optimization_convergence_and_failure(tmp_path, max_iter, converges):
    molecule = psiresp.Molecule(qcmol=distorted_water(), optimize_geometry=True)
    molecule.generate_conformers()
    conformer = molecule.conformers[0]
    initial = conformer.qcmol.geometry.copy()
    options = QMGeometryOptimizationOptions(max_iter=max_iter, convergence_tolerance=None)
    job = psiresp.Job(molecules=[molecule], qm_optimization_options=options, working_directory=tmp_path)
    with pytest.raises(QMTasksPending) as pending:
        job.optimize_geometries()
    assert not conformer.is_optimized
    path, = pending.value.paths
    if not converges:
        with pytest.raises(RuntimeError, match='Psi4 exited'):
            Smi2Gmx._run_psi4_input(os.environ['LIGTOOLS_PSI4'], path)
        with pytest.raises(ValueError, match='QM calculation failed'):
            job.optimize_geometries()
        assert not conformer.is_optimized
        np.testing.assert_allclose(conformer.qcmol.geometry, initial)
        return
    Smi2Gmx._run_psi4_input(os.environ['LIGTOOLS_PSI4'], path)
    result = options.read_output(distorted_water(), tmp_path)
    assert result.properties.optimization_iterations > 1
    assert result.properties.return_energy < result.trajectory_properties[0].return_energy - 0.01
    assert result.properties.final_max_force < 3e-5
    job.optimize_geometries()
    assert conformer.is_optimized
    assert np.max(np.abs(conformer.qcmol.geometry - initial)) > 0.1
    np.testing.assert_allclose(conformer.qcmol.geometry, result.final_molecule.geometry)
    timestamp = path.stat().st_mtime_ns
    job.optimize_geometries()
    assert path.stat().st_mtime_ns == timestamp


@pytest.mark.skipif(not os.environ.get('LIGTOOLS_PSI4'), reason='requires external Psi4 executable')
@pytest.mark.parametrize('threads', [1, 8])
def test_native_optimization_uses_new_default_tolerances(tmp_path, threads):
    options = QMGeometryOptimizationOptions(full_hess_every=-1)
    path = options.write_input(distorted_water(), tmp_path)
    Smi2Gmx._run_psi4_input(os.environ['LIGTOOLS_PSI4'], path, n_threads=threads)
    result = options.read_output(distorted_water(), tmp_path)
    assert result.success
    script = path.with_suffix('.opt.py')
    log = script.with_suffix('.out').read_text()
    import re
    row = re.findall(r'Convergence Criteria[^\n]+', log)[-1]
    assert row.count('1.00e-03') == 4
    assert re.search(rf'Threads:\s+{threads}\b', log)
    assert result.properties.optimization_iterations > 1


@pytest.mark.skipif(not os.environ.get('LIGTOOLS_PSI4'), reason='requires external Psi4 executable')
def test_native_single_point_uses_eight_threads(tmp_path):
    options = QMEnergyOptions()
    path = options.write_input(distorted_water(), tmp_path)
    Smi2Gmx._run_psi4_input(os.environ['LIGTOOLS_PSI4'], path, n_threads=8)
    data = qcel.util.deserialize(path.read_bytes(), 'msgpack')
    assert data['success']
    import re
    assert re.search(r'OpenMP threads:\s+8\b', data['stdout'])


@pytest.mark.skipif(not os.environ.get('LIGTOOLS_PSI4'), reason='requires external Psi4 executable')
def test_public_default_preserves_mmff_geometry_through_qm_and_resp(tmp_path, monkeypatch):
    initial = []
    generate = psiresp.Job.generate_conformers

    def capture(job):
        generate(job)
        initial.extend(conf.qcmol.geometry.copy() for mol in job.molecules for conf in mol.conformers)

    def unexpected(*args, **kwargs):
        raise AssertionError('The default pipeline must not call QM geometry optimization')

    monkeypatch.setattr(psiresp.Job, 'generate_conformers', capture)
    monkeypatch.setattr(psiresp.Job, 'optimize_geometries', unexpected)
    ligands = Smi2Gmx(smiles=['CO'], res_names=['MET'], work_dir=tmp_path,
                     n_conformer_pool=8, max_conf_n=[1], max_n_cores=8, n_conformer_threads=2)
    ligands.run_psi4(psi4_path=os.environ['LIGTOOLS_PSI4'], n_threads=8)
    assert initial
    assert not (tmp_path / 'optimization').exists()
    ligands.run_resp()
    for conf, before in zip(ligands.mols[0].conformers, initial):
        assert not conf.is_optimized
        np.testing.assert_array_equal(conf.qcmol.geometry, before)
        assert all(orientation.qc_wavefunction is not None for orientation in conf.orientations)
    assert np.isfinite(ligands.charges[0]).all()
    assert abs(sum(ligands.charges[0])) < 1e-5
