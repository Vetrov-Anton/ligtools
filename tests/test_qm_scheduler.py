"""Run real child executables to check scheduling, arguments and failures."""
import json
import os
from pathlib import Path
import sys

import pytest

from ligtools import Smi2Gmx
from ligtools.psiresp.qm import QMTasksPending


@pytest.fixture
def fake_psi4(tmp_path):
    if os.name == 'nt':
        pytest.skip('POSIX executable fixture')
    executable = tmp_path / 'executable with spaces' / 'fake psi4'
    executable.parent.mkdir()
    executable.write_text(f'#!{sys.executable}\n' + '''import json, sys, time
from pathlib import Path
path = Path(sys.argv[2])
start = time.time()
mode = path.read_text()
if mode == 'fail':
    print('simulated Psi4 failure', file=sys.stderr)
    sys.exit(7)
if mode != 'noresult':
    time.sleep(0.2)
    path.with_suffix('.done').write_text(json.dumps({
        'start': start, 'end': time.time(), 'argv': sys.argv[1:],
    }))
''')
    executable.chmod(0o755)
    return executable


def pending_stage(paths):
    def action():
        remaining = [path for path in paths if not path.with_suffix('.done').exists()]
        if remaining:
            raise QMTasksPending(remaining, 'manual-script.sh')
    return action


@pytest.mark.parametrize('parallel, expected_overlap', [(False, 1), (True, 2)])
def test_public_scheduler_and_resume(tmp_path, monkeypatch, fake_psi4, parallel, expected_overlap):
    import ligtools.smi2gmx as module

    paths = [tmp_path / f'job {i} $(touch PWNED).msgpack' for i in range(5)]
    for path in paths:
        path.write_text('success')

    class Job:
        def __init__(self, **kwargs):
            assert kwargs['n_processes'] == 2

        def generate_conformers(self):
            pass

        def optimize_geometries(self):
            pending_stage(paths[:3])()

        def generate_orientations(self):
            assert all(path.with_suffix('.done').exists() for path in paths[:3])

        def compute_orientation_energies(self):
            pending_stage(paths[3:])()

    monkeypatch.setattr(module.psiresp, 'Job', Job)
    ligands = Smi2Gmx(smiles=['O'], res_names=['HOH'], max_n_cores=2)
    ligands.run_psi4(psi4_path=fake_psi4, parallel=parallel)
    events = []
    for path in paths:
        record = json.loads(path.with_suffix('.done').read_text())
        assert record['argv'] == ['--qcschema', str(path), '--nthread', '1']
        assert path.with_suffix('.psi4.log').exists()
        events.extend([(record['start'], 1), (record['end'], -1)])
    running = maximum = 0
    for _, change in sorted(events):
        running += change
        maximum = max(maximum, running)
    assert maximum == expected_overlap
    assert not (tmp_path / 'PWNED').exists()
    # Finished calculations are reused rather than executed a second time.
    before = [path.with_suffix('.done').read_bytes() for path in paths]
    ligands.run_psi4(psi4_path=fake_psi4, parallel=parallel)
    assert before == [path.with_suffix('.done').read_bytes() for path in paths]


@pytest.mark.parametrize('mode, expected', [('fail', 'code 7'), ('noresult', 'without producing')])
def test_external_failure(tmp_path, fake_psi4, mode, expected):
    path = tmp_path / 'broken.msgpack'
    path.write_text(mode)
    ligands = Smi2Gmx()
    with pytest.raises(RuntimeError, match=expected):
        ligands._run_qm_stage(pending_stage([path]), str(fake_psi4), workers=1)
    if mode == 'fail':
        assert 'simulated Psi4 failure' in path.with_suffix('.psi4.log').read_text()


@pytest.mark.parametrize('error', [ValueError('invalid QM options'), KeyboardInterrupt(), SystemExit('stop')])
def test_unrelated_errors_are_not_swallowed(error):
    def action():
        raise error

    with pytest.raises(type(error)):
        Smi2Gmx()._run_qm_stage(action, 'unused executable', workers=2)


def test_missing_executable():
    with pytest.raises(FileNotFoundError, match='Psi4 executable'):
        Smi2Gmx().run_psi4(psi4_path='/nonexistent/ligtools-test-psi4')
