"""Exercise actual notebook execution with pre-existing native threads."""
import inspect
import json
import os
import sys

import pytest

from test_pipeline import test_four_ligand_pipeline as _pipeline


@pytest.mark.skipif(not os.environ.get('LIGTOOLS_FORCEFIELD'), reason='requires external chemistry tools')
def test_jupyter_without_thread_environment(tmp_path):
    nbformat = pytest.importorskip('nbformat')
    nbclient = pytest.importorskip('nbclient')
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpecManager

    kernels = tmp_path / 'kernels'
    spec = kernels / 'ligtools-test'
    spec.mkdir(parents=True)
    (spec / 'kernel.json').write_text(json.dumps({
        'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name': 'ligtools test', 'language': 'python',
    }))
    # Reuse the four-ligand assertions, executed inside the live notebook kernel.
    source = inspect.getsource(_pipeline)
    source = source[source.index('def test_four_ligand_pipeline'):]
    source = source.replace('parallel=False', 'parallel=True')
    warmup = '''import os, shutil, pytest
from pathlib import Path
import numpy as np
import parmed as pmd
import psi4
from ligtools import Smi2Gmx
assert get_ipython() is not None
limits = ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']
assert all(name not in os.environ for name in limits)
psi4.set_num_threads(2)
np.ones((256, 256)) @ np.ones((256, 256))
parent_threads = psi4.get_num_threads()
'''
    verification = '''assert all(name not in os.environ for name in limits)
assert psi4.get_num_threads() == parent_threads
print('NOTEBOOK PIPELINE PASSED: no environment settings, parent threads preserved')
'''
    notebook = nbformat.v4.new_notebook(cells=[
        nbformat.v4.new_code_cell(warmup),
        nbformat.v4.new_code_cell(source + '\ntest_four_ligand_pipeline(Path("calculation"))'),
        nbformat.v4.new_code_cell(verification),
    ])
    manager = KernelManager(kernel_name='ligtools-test', kernel_spec_manager=KernelSpecManager(
        kernel_dirs=[str(kernels)],
    ))
    environment = dict(os.environ)
    blocker = tmp_path / 'bin' / 'parallel'
    blocker.parent.mkdir()
    blocker.write_text(f'#!{sys.executable}\nfrom pathlib import Path\n'
                       f'Path({str(tmp_path / "parallel-called")!r}).touch()\nraise SystemExit(99)\n')
    blocker.chmod(0o755)
    environment['PATH'] = str(blocker.parent) + os.pathsep + environment.get('PATH', '')
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        environment.pop(key, None)
    try:
        nbclient.NotebookClient(notebook, km=manager, timeout=1200).execute(
            cwd=str(tmp_path), env=environment,
        )
    finally:
        nbformat.write(notebook, tmp_path / 'executed.ipynb')
    output = ''.join(o.get('text', '') for o in notebook.cells[-1].get('outputs', []))
    assert 'NOTEBOOK PIPELINE PASSED' in output
    assert not (tmp_path / 'parallel-called').exists()
