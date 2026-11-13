# Local backend validation — 2026-11-13

## Final validation before committing

The current checkout passed 126 tests on each of Python 3.11, 3.12, 3.13 and
3.14; ten optional integrations were skipped in those minimal environments.
The conda environment with Psi4, AmberTools and notebook dependencies passed
all 136 tests in 69.73 seconds, including native optimization/thread checks,
the four-ligand pipeline, its real Jupyter execution, parallel conformer
selection/pruning and all CoordMatcher modes. Both wheel and sdist built.

README, standalone and notebook examples now describe MMFF94 geometries by
default, optional QM optimization and explicit eight-thread Psi4 execution.
The Russian PDF report and its build script are under `reports/`. Generated
QM working files are excluded from version control. Logs for these checks are
in `ligtools-audit/local-optimization/commit-tests-*.log`,
`commit-conda-full.log` and `commit-build.log` beside the checkout.

## API, MPO and MPI full pipeline on MMFF94 geometries

A fresh run of the three GitHub-example ligands completed from SMILES through
MMFF94, eight-thread Psi4 single points, ESP/RESP, parallel ACPYPE and ITP/GRO/
PDB/RTP export, including RTP insertion into an isolated amber99sb-ildn.ff copy.
Each molecule used 4000 initial candidates and a maximum of five selected
conformers. MMFF94 parameters were available for all three inputs.

| Ligand | SMILES | Atoms | Selected conformers | RESP net charge | Result |
| --- | --- | ---: | ---: | ---: | --- |
| API | `CC(=O)[PH](=O)[O-]` | 10 | 1 | -1 | Passed |
| MPO | `CC(=O)P(OC)(=O)[O-]` | 14 | 4 | -1 | Passed |
| MPI | `CC(=O)P(C)(=O)[O-]` | 13 | 1 | -1 | Passed |

No QM geometry optimization ran. Conformer coordinates were exactly preserved
through QM/RESP. Exported GRO/PDB geometries matched MMFF94 within aligned RMSD
0.00472, 0.00525 and 0.00402 angstrom respectively, consistent with output
rounding. Checks also covered atom counts/names/order, finite coordinates and
charges, ITP/RTP charge agreement with RESP, charge sums and byte-identical RTP
copies in the force field.

Total wall time was 44.192 seconds: conformer preparation 25.965 s; preparation
plus six fresh QM calculations 36.291 s; ESP/RESP 2.203 s; ACPYPE 5.454 s;
export/forcefield insertion 0.069 s. These timings overlap only where indicated
(the preparation time is included in the QM-stage total).

Reproducible script, logs, result JSON and final force field:
`ligtools-audit/api-mpo-mpi-mmff/` beside the checkout. Files
`amber99sb-ildn.ff/API.rtp`, `MPO.rtp` and `MPI.rtp` are the validated residues.

## Completed TPP(-2) pipeline on MMFF94 geometries

The complete TPP(-2) pipeline passed from SMILES through MMFF94, eight-thread
Psi4 single points, ESP, RESP, ACPYPE/GAFF2 and RTP insertion into an isolated
amber99sb-ildn.ff copy. The result has 42 atoms and charge -2 in fitted charges,
ITP and RTP. Atom names/order and RTP copying were checked. No Psi4 optimization
directory was created; conformer coordinates remained exactly unchanged through
QM and RESP. Exported GRO geometry has aligned RMSD 0.00498 angstrom against
the QM/MMFF structure, within output rounding.

The initial successful QM calculation exposed a loading bug: QCWaveFunction
inherited BaseMolecule's RDKit graph reconstruction, which failed on the returned
TPP aromatic graph without mapped SMILES metadata. That hook is now bypassed
for wavefunctions; only coordinates, basis and orbital data are required.
The cached real TPP wavefunction (434 by 434 density matrix) was read successfully
without RDKit conversion. A regression test explicitly rejects any such call.
The full unit suite passed 126 tests with ten optional integrations skipped.

The successful repeat used the existing QM result: MMFF preparation 41.370 s,
cached QM read 0.005 s, ESP/RESP 9.264 s, ACPYPE 30.661 s, export/forcefield
insertion 0.034 s, total 81.529 s. The initial QM calculation including the failed
read took 38.367 s separately. These are measured attempts, not an 81.529 s cold
run. Timing instrumentation in the audit script was restored before spawning
ACPYPE workers to make its ligand object serializable.

ACPYPE currently falls back from PDB/user charges to an internal AM1-BCC job;
the exported ITP/RTP charges are subsequently replaced with and verified against
RESP charges. Exported geometry preservation was checked independently above.

Artifacts and reproducible runner are in
`ligtools-audit/tpp-minus2-mmff-threads8/` beside the checkout. The final RTP is
`amber99sb-ildn.ff/TPP.rtp`; logs and timing JSON record the completed run.
Previously stopped QM-optimization attempts below remain historical records.

## MMFF94 geometry as the default

At the user's request, Smi2Gmx now defaults to `optimize_geometry=False`.
Conformers are generated, pruned, minimized in MMFF94 and selected, then used
directly for Psi4 orientation wavefunctions, ESP and RESP. The QM optimization
stage is skipped entirely; it remains available with `optimize_geometry=True`.
MMFF conformers are not falsely marked as QM-optimized.

A native methanol regression check captured coordinates after MMFF94 selection,
ran Psi4 with eight threads and fitted RESP, and confirmed exact preservation
of conformer coordinates, populated wavefunctions, finite charges, neutral net
charge and absence of an optimization directory. An unexpected call to the QM
optimizer is an explicit test failure. All 14 local-optimization/native checks
passed. The full TPP run remains stopped; it was not restarted automatically.
Logs: `ligtools-audit/local-optimization/mmff-only-tests.log` and
`mmff-only-unit.log` beside the checkout.

## Native Psi4 thread control

`Smi2Gmx.run_psi4(n_threads=1)` now exposes the native thread count for each
external Psi4 process. The scheduler limits concurrent jobs to
`max_n_cores // n_threads`; sequential mode retains the requested native thread
count. Invalid/noninteger thread counts and requests exceeding the total budget
are rejected before starting calculations. The default remains one thread.
Thread limits for spawned ESP/ACPYPE workers do not constrain the external
Psi4 process. Regression tests check executable arguments, overlapping job
counts, resumed results and failure handling. Native tests inspect Psi4's
reported threads during both optimization and single-point calculations.

## Updated geometry convergence defaults

At the user's request, local optimization now defaults to explicit maximum/RMS
force and maximum/RMS displacement thresholds of 1e-3 in OptKing's internal
atomic units. A common `convergence_tolerance` option supplies these defaults;
explicit individual Psi4 keywords take precedence case-insensitively. Setting
it to None restores the named preset (gau_tight by default). The energy-change
criterion remains inactive; iteration count and Hessian settings are unchanged.
The changed serialized options also change the optimization cache key.

All 11 local-optimization checks passed with native Psi4, including an actual
optimization whose output printed all four thresholds as 1.00e-03, strict-preset
convergence/nonconvergence checks, and explicit-override validation. Historical
water and TPP measurements below used the earlier strict preset. The stopped
TPP run was not resumed with these new defaults.

## Parallel conformer selection and bottleneck audit

The TPP(-2) QM run was stopped at the user's request. Its conformer preparation
took 111.791 seconds; it retained one conformer, so only one Psi4 process ran.
No successful full TPP parameterization is claimed.

An independent controlled benchmark generated 4000 TPP(-2) candidates with four
RDKit threads and applied MMFF94 before comparing old and new selection on the
same pool. Embedding took 91.507 seconds, MMFF94 took 20.765 seconds, and the
original selection took 0.337 seconds. Updated selection took 0.118 seconds
with one worker and 0.086 seconds with four. The selected IDs matched exactly.
Energy scoring matched the original formula: 0.346 seconds originally versus
0.047 seconds using four workers. These figures describe this benchmark pool,
not a resumed full parameterization or a universal speedup.

Selection now avoids the full O(N^2) RMSD matrix and computes at most k-1 rows
for k retained conformers, with O(N) selection storage. Larger candidate sets
use spawned workers; small sets use RDKit's native threading. Energy scoring
uses bounded blocks rather than allocating distances for the entire pool.
Tests cover equivalence to full-matrix selection, symmetry-aware RMSD for
reordered conformer IDs, energy equivalence, pose preservation, thread-limit
propagation, and a real Jupyter kernel without environment-variable setup.
The full unit suite passed 87 tests with six optional tests skipped; the conda
environment separately passed all eight selection checks including Jupyter.

The earlier audit identified a limitation: a diagnostic run with embedding-time pruning disabled
took 16.857 seconds for 4000 candidates, identifying the built-in greedy RMSD
pruning as the main overhead in this example. Production pruning was preserved
to avoid silently changing the initial pool. This code has not eliminated that
serial section of RDKit at that point. This is superseded by the separate
parallel-pruning implementation described below. Molecules are prepared in sequence to share the worker
budget; QM and ESP run across conformers/orientations, ACPYPE across molecules.
The small coupled RESP fit and greedy selection decisions remain sequential.
With one available task, the relevant external-process stage uses one worker.

Scripts, binary pool, timing JSON and logs are in the adjacent directory
`ligtools-audit/conformer-selection/`. The Jupyter artifact is
`/tmp/ligtools-selection-notebook/test_parallel_selection_in_jup0/executed.ipynb`.

## Separate parallel embedding-time pruning

`_pruning.py` now replaces the serial RDKit pruning pass. Embedding runs with
pruning disabled, after which worker batches check the accepted prefix and
earlier candidates in the same batch. Acceptance is resolved in original order;
matches to rejected candidates cannot discard later candidates. Initial
conformers are preserved, and newly accepted IDs match native RDKit numbering.
The native point-alignment routine, strict SSD cutoff, heavy-atom mappings,
terminal-group symmetry and 1000-match cap are retained. A conservative
singular-value lower bound skips provably distant pairs before exact alignment.
Only bounded batch results are stored, without a full distance/adjacency matrix.

Exact comparisons against native embedding passed on flexible, aromatic,
carboxylate, nitro, phosphate, chiral and isotopically labelled molecules with
one and four workers, including existing non-contiguous conformer IDs and
multiple worker batches. All 21 pruning tests passed on Python 3.11–3.14.
The full unit suite passed 108 tests with six optional integrations skipped.
The conda environment separately passed 29 selection/pruning tests, including
parallel pruning and selection in a real Jupyter kernel. Wheel and sdist built.

TPP(-2), 4000 candidates, RMSD threshold 0.05 angstrom, four embedding threads:

| Operation | Wall time (seconds) |
| --- | ---: |
| Native embedding + pruning | 119.895 |
| Embedding without native pruning | 18.764 |
| External pruning, one worker | 0.969 |
| External pruning, four workers | 2.448 |
| New embedding + four-worker pruning | 21.211 |
| Subsequent MMFF94, four threads | 22.673 |
| Subsequent final selection, four workers | 0.231 |

All 4000 retained coordinates, IDs and their order matched native RDKit exactly;
the selected IDs also matched the prior controlled benchmark. These are single
wall-clock measurements on a shared workstation (some regression tests overlapped
the native baseline), not isolated performance guarantees. The earlier baseline
was 91.507 seconds. Process startup makes the four-worker pruning pass slower
than the single-worker pass on this particular pool after lower-bound screening;
the principal gain is avoiding unnecessary exact alignments. Expensive alignment
work is parallel, while ordered acceptance and batch boundaries remain serial.
No full TPP QM/RESP/ACPYPE run was resumed.

The reproducible benchmark and results are
`ligtools-audit/conformer-selection/compare_pruning.py` and
`ligtools-audit/conformer-selection/pruning_comparison.json` beside the checkout.
The executed notebook is under `/tmp/ligtools-pruning-notebook/`.

## Real local geometry optimization

The former local path evaluated one gradient and incorrectly marked the input
geometry as optimized. It now runs iterative `psi4.optimize`, reads the final
geometry only after convergence, and rejects single-gradient results. A changed
cache key excludes files produced by the old path. Existing in-memory objects
must be recreated after upgrading.

Native Psi4 1.11 checks use distorted water with HF/6-31G*. Optimization took
12 steps: energy decreased from -75.80612411545192 to -76.01072065750574 hartree,
and the final maximum gradient component was 7.58e-8 hartree/bohr. Coordinates
changed and the converged cache was reused. A one-iteration limit deliberately
failed and left the original geometry and optimization flag unchanged.

A real Jupyter kernel completed methanol through optimization, single points,
ESP, two-stage RESP, parallel ACPYPE and PDB/GRO/ITP/RTP exports, with atom-count
and charge-sum checks. This also verifies that `run_resp()` reuses the prepared
orientations instead of requesting a second QM calculation. The charge method
remains gas-phase HF/6-31G* RESP; no new charge protocols were added.

The updated unit suite passed 80 tests, with five optional integrations skipped.
With `PYTHONPATH=src` and `LIGTOOLS_PSI4` pointing to the conda executable,
all five local-optimization regression cases passed in 9.89 seconds, including
the native convergence and deliberate nonconvergence cases.
Logs, the notebook runner and executed notebook are in the adjacent
`ligtools-audit/local-optimization/` directory. The more expensive four-ligand
notebook audit was stopped; the earlier four-ligand results below predate this
fix and do not establish convergence of real geometry optimization.

## CoordMatcher audit and simplified interface

All five alignment modes were tested on transformed reference poses, both with
explicit PDB hydrogens and with heavy atoms only. Checks cover exact coordinate
transfer, rigid-body distance preservation, shape overlap, output writing in all
modes, atom names/order/counts, multiple probes, default modes and repeated calls.
Conformer-search tests include reordered atoms, aromatic, charged and chiral
molecules, PDBs without bond orders, and different probe/reference ligands.
Complex rebuilding checks replacement of the old ligand, retention of scaffold,
water and custom ions, their coordinates and the periodic box.

Fixed issues include missing output files in modes 0/1, ambiguous SMILES mapping
when valid chemical bond orders are available, and unclear errors for unmatched
structures, missing SMILES, failed embedding or mismatched argument lists.
MCS alignment uses heavy atoms with a bounded search time. PDB atom sets are
preserved throughout. Topological fallback can still give multiple equivalent
assignments for symmetric atoms when bond-order information is absent.

At the user's request, `addHs_to_probe`, `addHs_to_ref` and `opt_geom` were
removed from the API and implementation. There is no PDB hydrogen addition or
MMFF/UFF geometry optimization. Modes 2/4 still embed conformers using a temporary
SMILES molecule with internal explicit hydrogens, then map only the original
probe atoms back to the output. Mode 0 copies reference geometry, modes 1/3
preserve the probe geometry by rigid transforms, and modes 2/4 replace it with
the selected embedded conformer.

45 CoordMatcher regression tests passed on Python 3.11–3.14. The full unit suite
passed 77 tests on those interpreters (three optional integrations skipped).
An additional real-ipykernel test passed all five modes with drawing enabled,
PDB exports and complex rebuilding, without importing Psi4 or the RESP backend.
All five modes also completed on actual PYR, API, MPO and MPI PDB outputs of the
parameterization pipeline: 20 successful checks with finite coordinates and
preserved atom counts. These checks establish tested behavior, not global
optimality of shape/conformer searches or sampling convergence.

Audit logs and the separate pipeline-ligand script are in the adjacent
`ligtools-audit/coordmatcher/` directory. The executed notebook is under
`/tmp/ligtools-coordmatcher-notebook/test_all_coordmatcher_modes_in0/executed.ipynb`.

## Fresh conda installation audit

Both installation paths were executed online with conda's libmamba solver on
Linux x86_64. Existing environments were not modified. The README's named
environment was replaced only by a separate audit prefix to avoid name collisions.

| Installation path | Isolated prefix | Result |
| --- | --- | --- |
| README: conda create + pip install from GitHub | `/tmp/ligtools-install-readme` | Installed; 33 tests passed in 57.70 s |
| conda env create -f environment.yml | `/tmp/ligtools-install-yaml` | Installed; 33 tests passed in 60.09 s |

The Git installation resolved to commit
`4d4615decd86b4b7a1bc71d33d54c1ca2c4adc24`. The YAML installation built the local
checkout at that same commit. Imports were checked from outside the checkout;
`ligtools.__file__` pointed inside the tested environment's site-packages.
Psi4, antechamber, parmchk2, tleap and obabel also resolved inside each prefix.

Both environments selected Python 3.13.18, Psi4 1.11, AmberTools 26.0, NumPy
2.5.4, Pydantic 2.14.3, QCElemental 0.51.4, SciPy 1.18.1, RDKit 2026.09.2,
ParmEd 4.3.1 and MDAnalysis 2.11.0. Installed local backends were
`psiresp 0.4.2+ligtools.1` and `acpype 2023.10.27+ligtools.1`.

The initial README commands did not install ipykernel. Installing it with
`python -m pip install ipykernel` succeeded; kernel registration was tested with
an isolated `--prefix /tmp/ligtools-kernel-registration` rather than modifying the
user's global kernels. README now includes notebook kernel installation and
registration. The YAML already includes ipykernel.

Baseline import/executable/conformer checks ran before installing test tools.
For the test suite, pytest, nbclient and nbformat were subsequently installed
(and ipykernel in the README environment). Tests exercised the installed package,
not an editable checkout. Each 33-test run included the actual Jupyter
four-ligand integration test with pool size 20, RDKit threads 2, parallel QM,
parallel ACPYPE and all output/charge checks. The duplicate standalone integration
test was excluded because its same assertions run inside that notebook.

`pip check` is **not clean**: initially it reported nine issues for the README
environment and seven for the YAML environment. They concern bundled AmberTools
auxiliary programs: fetkutils, proprep, ndfes, edgembar and packmol-memgen, with
NumPy/Biopython constraints or missing netCDF4/pdb2pqr. No ligtools dependency
conflicts were reported. Successful pipeline execution does not validate those
other tools; no package metadata was removed to hide the conflicts.

| Auxiliary package | Requirement | Observed installation |
| --- | --- | --- |
| fetkutils, ndfes | NumPy <2; netCDF4 | NumPy 2.5.4; netCDF4 missing in the manual environment, present in the YAML environment |
| edgembar | NumPy <2 | NumPy 2.5.4 |
| proprep | NumPy <2; Biopython >=1.83,<1.86; pdb2pqr >=3.5 | NumPy 2.5.4, Biopython 1.88, pdb2pqr missing |
| packmol-memgen | pdb2pqr | Missing in both environments |

These requirements may prevent the auxiliary programs from working; they are
not treated as harmless warnings. The antechamber/parmchk2/tleap workflow used
by ligtools passed the integration test, while these other programs remain
unverified.

Commands, complete conda/pip logs, baseline checks, test logs and final package
snapshots are in the adjacent `ligtools-audit/installation/` directory.

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
