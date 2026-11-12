# Local backend validation — 2026-11-13

## Python QM scheduler

GNU parallel is no longer a runtime dependency. Python launches Psi4 with
argument lists through subprocess, using a bounded ThreadPoolExecutor to wait
for independent external calculations. Each Psi4 process explicitly uses one
thread; `max_n_cores` bounds simultaneous calculations. `parallel=False` runs
QM commands sequentially. Shell scripts are not parsed or executed by this path.

32 unit/numerical tests passed on Python 3.11–3.14. New tests launch real stub
executables and check concurrency, sequential execution, staged ordering, reuse
of completed calculations, paths containing spaces and shell metacharacters,
nonzero exits, missing output, missing executables and propagation of unrelated
exceptions and interrupts.

The full four-ligand Jupyter test passed in 56.95 seconds without thread-limit
environment variables. A failing stand-in for the `parallel` executable was
placed first on PATH; it was never invoked. All pipeline output checks passed.
The run log is in `ligtools-audit/python-qm-notebook.log` beside this repository.

## Notebook-safe parallel execution

ESP and parallel ACPYPE now use explicit spawn contexts and importable library
targets. Native thread limits are configured automatically inside workers only;
the notebook's environment and Psi4 thread count are preserved. The standalone
example has a main guard; the notebook example requires no such guard.

24 unit/numerical tests passed on Python 3.11–3.14, including a live parent
thread, warmed-up BLAS and a real spawned worker whose native thread pools are
limited automatically. Cached ESPs do not create a pool.

`tests/test_notebook.py` passed in 57.25 seconds in a real ipykernel with
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS absent. The kernel
imported Psi4, set its parent thread count to two and warmed up NumPy before
executing all stages for all four example ligands (initial pool 20, RDKit threads
2), including **parallel** ACPYPE. The test checked atom counts, charge sums,
coordinates, output files and force-field insertion, then verified that the
parent environment and Psi4 thread count were unchanged.

The executed notebook was saved under
`/tmp/ligtools-notebook-validation/test_jupyter_without_thread_en0/executed.ipynb`;
the run log is in the adjacent `ligtools-audit/notebook-validation.log`.
Older thread-limited runs below document earlier validation, not a current
requirement to set environment variables.

## Conformer controls

After adding `n_conformer_pool` and `n_conformer_threads`, 22 tests passed on
Python 3.11–3.14 (the external integration test was skipped in these venvs).
The tests exercise real multithreaded RDKit embedding, scalar and per-molecule
pool sizes, the selection limit, UFF/MMFF thread limits and invalid inputs.
The synthetic charge regression now fixes its random seed to avoid dependence
on the ordering of earlier RDKit calls.

All 23 tests passed in 53.72 seconds with Psi4 1.11 and AmberTools 26, including
all four ligands with `n_conformer_pool=20` and `n_conformer_threads=2` supplied
through the public Smi2Gmx constructor. Coordinate files, topology charge sums
and force-field insertion were checked by that integration test.

## Python packages

Clean virtual environments were installed from the local project, without
upstream psiresp or acpype. Six unit/numerical tests passed on each interpreter:

| Python | Result |
| --- | --- |
| 3.11.15 | 6 passed, external integration skipped |
| 3.12.3 | 6 passed, external integration skipped |
| 3.13.14 | 6 passed, external integration skipped |
| 3.14.6 | 6 passed from the built wheel, external integration skipped |

`uv pip check` passed in all four virtual environments. The Python 3.13 venv
used Pydantic 2.14.3, QCElemental 0.51.4, NumPy 2.5.4 and RDKit 2026.09.2.
The GitHub Actions matrix is configured but was not run remotely.

Both sdist and wheel built successfully. The wheel contains both local backends
and their license notices, excludes upstream setup.py/Versioneer and Amber
binaries, and has no dependency on the upstream psiresp/acpype distributions.

## External chemistry tools

All seven tests, including `tests/test_pipeline.py`, passed in an isolated
Python 3.13.15 conda environment with Psi4 1.11 and AmberTools 26.0. That test
runs all four example ligands through QM, ESP/charge fitting, topology generation,
all output formats and force-field injection into a temporary copy.

For this integration test, the initial conformer pool was 20 rather than the
example's default 4000, with at most two simultaneous QM commands. It tests all
stages but does not establish the scientific convergence of conformer sampling.
Atom counts, finite coordinates/charges, net charge -1, consistency between fitted
charges and ITP charges, and copying of RTP files were verified.

The same stages also completed using the older Psi4 1.10.2 executable and the
new local backends with Pydantic 2.14.3 and QCElemental 0.51.4.

For an additional numerical regression, the maintained and original RESP
implementations fitted **identical saved ESP data** for PYR, API, MPO and MPI.
The maximum charge difference was 0 for every ligand. Comparing independent runs
with different RDKit releases is not equivalent: generated conformers and final
charges can differ (the earlier independent PYR runs differed by about 0.051 e).

## Known environment limitation

The conda AmberTools 26 distribution also installs Python metadata for auxiliary
programs that are not used by this pipeline. In the isolated conda environment,
`uv pip check` reported 11 issues involving edgembar, fetkutils, proprep,
packmol-memgen, ndfes and nlopt (NumPy/Python constraints, missing optional
dependencies, or incomplete wheel metadata). These were not repaired or removed
as part of maintaining the two included backends. Thus the successful pipeline
test is not a claim that every auxiliary AmberTools program is compatible.

Reproduce the checks using the commands in THIRD_PARTY.md. Local run logs are
under the adjacent `ligtools-audit/modern-psi411/` and
`ligtools-audit/modern-pipeline/` directories; these logs are not packaged.
