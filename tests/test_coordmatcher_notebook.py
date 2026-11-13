import importlib.util
import json
import os
from pathlib import Path
import sys

import pytest


@pytest.mark.skipif(not all(importlib.util.find_spec(m) for m in ['nbclient', 'ipykernel']),
                    reason='requires notebook test dependencies')
def test_all_coordmatcher_modes_in_jupyter(tmp_path):
    import nbclient
    import nbformat
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpecManager

    kernels = tmp_path / 'kernels'
    spec = kernels / 'coordmatcher-test'
    spec.mkdir(parents=True)
    (spec / 'kernel.json').write_text(json.dumps({
        'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name': 'CoordMatcher test', 'language': 'python',
    }))
    code = '''from pathlib import Path
import sys
import numpy as np
import parmed as pmd
from ligtools import CoordMatcher
from test_coordmatcher import make_pair, SMILES
assert get_ipython() is not None
probe, reference = make_pair(Path('inputs'))
for mode in range(5):
    matcher = CoordMatcher([probe], [reference], [SMILES])
    target = Path(f'aligned-{mode}.pdb')
    matcher.align(modes=[mode], n_confs=12, n_cores=2, outputs=[target], draw=True)
    assert target.exists()
    assert np.isfinite(matcher.pmd_objs[0].coordinates).all()
    matcher.get_pdbs([Path(f'export-{mode}.pdb')])
    complex_ = matcher.get_complex(output=Path(f'complex-{mode}.pdb'),
                                   scaffold=reference, pdb_names=['LIG'])
    assert len(complex_.atoms) == len(matcher.pmd_objs[0].atoms)
    print('PASS MODE', mode)
assert 'psi4' not in sys.modules
assert 'ligtools.psiresp' not in sys.modules
'''
    notebook = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(code)])
    manager = KernelManager(kernel_name='coordmatcher-test', kernel_spec_manager=KernelSpecManager(
        kernel_dirs=[str(kernels)],
    ))
    env = dict(os.environ)
    root = Path(__file__).resolve().parents[1]
    env['PYTHONPATH'] = os.pathsep.join([str(root / 'src'), str(root / 'tests')])
    try:
        nbclient.NotebookClient(notebook, km=manager, timeout=60).execute(cwd=str(tmp_path), env=env)
    finally:
        nbformat.write(notebook, tmp_path / 'executed.ipynb')
    outputs = notebook.cells[0].outputs
    text = ''.join(o.get('text', '') for o in outputs)
    assert all(f'PASS MODE {mode}' in text for mode in range(5))
    assert any('image/png' in o.get('data', {}) for o in outputs)
