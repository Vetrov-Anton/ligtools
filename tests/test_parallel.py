from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import os
import threading

import numpy as np
from threadpoolctl import threadpool_info

from ligtools._parallel import initialize_worker
from ligtools import psiresp


def _worker_state():
    return multiprocessing.get_start_method(), threadpool_info()


def test_spawn_with_active_parent_threads():
    before = dict(os.environ)
    method = multiprocessing.get_start_method()
    stop = threading.Event()
    thread = threading.Thread(target=stop.wait)
    thread.start()
    try:
        np.ones((128, 128)) @ np.ones((128, 128))
        with ProcessPoolExecutor(
            max_workers=1, mp_context=multiprocessing.get_context('spawn'),
            initializer=initialize_worker,
        ) as pool:
            worker_method, libraries = pool.submit(_worker_state).result(timeout=30)
            assert worker_method == 'spawn'
            assert libraries
            assert all(lib['num_threads'] == 1 for lib in libraries)
    finally:
        stop.set()
        thread.join()
    assert dict(os.environ) == before
    assert multiprocessing.get_start_method() == method


def test_no_pool_for_cached_esps(monkeypatch):
    from ligtools.psiresp import job as module

    def unexpected_pool(**kwargs):
        raise AssertionError('No pending ESPs: no pool should start')

    monkeypatch.setattr(module, 'ProcessPoolExecutor', unexpected_pool)
    molecule = psiresp.Molecule.from_smiles('O')
    molecule.generate_conformers()
    molecule.conformers[0].add_orientation_with_coordinates(molecule.coordinates)
    molecule.conformers[0].orientations[0].esp = np.ones(3)
    psiresp.Job(molecules=[molecule]).compute_esps()
