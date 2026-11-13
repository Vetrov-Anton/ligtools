"""Standalone helper executed inside the selected Psi4 executable's interpreter."""
import os
from pathlib import Path
import tempfile
import traceback


def run_optimization(input_path):
    import numpy as np
    import psi4
    import qcelemental as qcel

    path = Path(input_path)
    data = qcel.util.deserialize(path.read_bytes(), 'msgpack')
    try:
        specification = data['specification']['specification']
        molecule = psi4.core.Molecule.from_schema(data['initial_molecule'])
        keywords = dict(specification['keywords'])
        pcm_input = keywords.pop('pcm__input', None)
        keywords['basis'] = specification['model']['basis']
        psi4.set_options(keywords)
        if pcm_input:
            psi4.pcm_helper(pcm_input)
        energy, wavefunction, history = psi4.optimize(
            specification['model']['method'], molecule=molecule,
            return_wfn=True, return_history=True,
        )
        energies = [float(value) for value in history['energy']]
        gradients = [np.asarray(value) for value in history['gradient']]
        gradient = gradients[-1]
        final_molecule = dict(data['initial_molecule'])
        final_molecule['geometry'] = np.asarray(wavefunction.molecule().geometry())
        if not all(np.isfinite(value).all() for value in
                   [np.asarray(energies), gradient, final_molecule['geometry']]):
            raise ValueError('Optimizer returned non-finite results')
        result = {
            'schema_name': 'qcschema_optimization_result', 'schema_version': 2,
            'input_data': data, 'final_molecule': final_molecule,
            'trajectory_results': [],
            'trajectory_properties': [
                {'return_energy': e, 'return_gradient': g} for e, g in zip(energies, gradients)
            ],
            'properties': {
                'return_energy': float(energy), 'return_gradient': gradient,
                'optimization_iterations': len(energies),
                'final_max_force': float(np.abs(gradient).max()),
                'final_rms_force': float(np.sqrt(np.mean(gradient**2))),
            },
            'success': True,
            'provenance': {'creator': 'Psi4', 'version': psi4.__version__, 'routine': 'psi4.optimize'},
            'extras': {'ligtools_optimizer': 'psi4-optimize-v1', 'converged': True},
        }
    except Exception as error:
        result = {
            'success': False,
            'error': {'error_type': type(error).__name__, 'error_message': str(error)},
            'extras': {'ligtools_optimizer': 'psi4-optimize-v1', 'converged': False},
        }
        traceback.print_exc()
        _write_result(path, result)
        raise
    _write_result(path, result)


def _write_result(path, result):
    import qcelemental as qcel

    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=path.name, delete=False) as file:
        temporary = file.name
        file.write(qcel.util.serialize(result, 'msgpack-ext'))
    os.replace(temporary, path)
