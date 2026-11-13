# Included backends

Imported and modified on 2026-11-13. These are maintained forks, packaged only
under `ligtools.*`. PsiRESP author: Lily Wang. ACPYPE author: Alan Wilter
Sousa da Silva and contributors. Original copyright notices and licenses remain in each module.
The original ligtools files retain their MIT license. The combined distribution
includes GPL-3.0 ACPYPE; its package metadata therefore declares GPL-3.0-only.

| Module | Upstream release | Source | Source archive SHA-256 | License |
| --- | --- | --- | --- | --- |
| `ligtools.psiresp` | 0.4.2 | https://pypi.org/project/psiresp/0.4.2/ ; https://github.com/lilyminium/psiresp | `7fb069351824df32a5064bc4061877d35cf8907c7ef6aa6d83a4af0797b72ae4` | LGPL-3.0, `src/ligtools/psiresp/LICENSE` |
| `ligtools.acpype` | 2023.10.27 | https://pypi.org/project/acpype/2023.10.27/ ; https://github.com/alanwilter/acpype | `2041ef01031015b6901aabce75cf39c17c6d7e1a034f8b0f4d168f0a6bd06a99` | GPL-3.0, `src/ligtools/acpype/LICENSE` |

Only runtime Python modules and licenses were copied from the PyPI source
archives. Upstream build systems, Versioneer, tests/data, and bundled Amber
binaries are excluded. Psi4, AmberTools and Open Babel remain external programs.

## Local changes

- ACPYPE absolute imports now resolve inside `ligtools.acpype`.
- AmberTools is found on PATH. There is no fallback to bundled executables and
  no replacement of LD_LIBRARY_PATH or DYLD_LIBRARY_PATH.
- PsiRESP models and validators use native Pydantic 2 APIs. QCElemental models
  explicitly use v2; v1 file conversion is isolated in the Psi4 single-point
  reader/writer. Local optimization uses QCSchema v2 inputs and results.
  Python 3.14 uses QCElemental's Pydantic 2 compatibility implementation of the
  v1 wire format because Pydantic 1 is unsupported on that interpreter.
- NumPy `in1d` is replaced by `isin`; importing ligtools no longer patches NumPy.
- Conformer generation exposes initial pool size and RDKit thread count through
  Smi2Gmx; the thread limit also applies to optional UFF/MMFF minimization.
- MMFF94 minimization is enabled by default for the generated psiresp conformer
  pool before selection. Explicit `minimize=None` disables it in backend options.
- Conformer energy scoring uses bounded, threaded blocks. Greedy selection
  computes only RMSD rows for selected candidates instead of the full matrix;
  larger sets use a spawn process pool and small sets use native RDKit threads.
  Native full-matrix calls also accept the thread limit. Alignment calculations
  preserve coordinates and symmetry handling.
- `_pruning.py` replaces RDKit's serial embedding-time pruning with bounded
  spawned-worker batches and ordered acceptance. The algorithm follows RDKit's
  `getMolSelfMatches` / `_isConfFarFromRest` in
  https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/DistGeomHelpers/Embedder.cpp
  (RDKit contributors, BSD-3-Clause). It uses the same native point-alignment
  routine, heavy-atom and terminal-group symmetry matching, 1000-match cap,
  strict cutoff, existing-conformer handling and output IDs. A singular-value
  lower bound avoids unnecessary alignments without changing the acceptance
  rule. Terminal-group matching follows the installed RDKit's default option.
- The sparse solver uses public SciPy `spsolve`, retaining the least-squares
  fallback for singular systems. NumPy block arrays replace deprecated matrices.
- Removed the ineffective assignment that replaced the RESP2 configuration class
  with a boolean; the existing pipeline still uses Job with its existing defaults.
- Replaced the local optimizer's single gradient evaluation with iterative
  `psi4.optimize` through `_optimization.py`. Only converged results supply
  final coordinates; failures propagate and old gradient caches are excluded.
  RESP fitting reuses the orientations and wavefunctions already prepared by
  `run_psi4`, avoiding a second orientation-generation and QM cycle.
- Geometry convergence defaults now explicitly set maximum/RMS force and
  displacement thresholds to 1e-3 in internal-coordinate atomic units.
  `convergence_tolerance=None` restores named preset thresholds; individual
  explicit Psi4 convergence keywords override the common tolerance.
- Smi2Gmx now leaves QM geometry optimization disabled by default and computes
  wavefunctions/ESP on selected MMFF94 geometries. `optimize_geometry=True` on
  construction restores the optional QM optimization stage. The backend's
  `is_optimized` flag is reserved for successful QM optimization.
- QCWaveFunction no longer runs the BaseMolecule RDKit initialization hook.
  Reading a QM wavefunction uses only QCSchema coordinates, basis and orbitals;
  no bond, formal-charge or aromaticity reconstruction is needed. The original
  ligand's chemical graph remains available independently from its input SMILES.
- ESP and ACPYPE workers use a local spawn context with worker-local native
  thread limits. Notebook users do not need BLAS environment variables.
- QM inputs are passed as structured paths to Python's subprocess scheduler,
  without shell parsing or GNU parallel. Only the specific pending-QM signal
  triggers external execution; other errors propagate. Process exit codes and
  returned QCSchema results are checked. Upstream shell scripts are retained
  only as an optional aid for manual execution.
- `run_psi4(n_threads=...)` sets native threads for each external Psi4 process
  in both optimization and single-point modes. Concurrent QM jobs are bounded
  by `max_n_cores // n_threads`; ESP/ACPYPE worker-local thread limits remain
  independent of the external Psi4 executable.
- The public `Smi2Gmx` and `CoordMatcher` imports are preserved. The geometry
  defaults described above are intentional changes; the RESP fitting method is
  unchanged. Optional legacy QCFractal integrations are not covered by our tests.

## Maintenance

Make changes in these directories directly and record them here. Compare future
upstream fixes against the source versions above rather than overwriting the
forks. Do not reintroduce top-level `psiresp`/`acpype` imports or their packaging.
Preserve license notices when redistributing modified source or wheels.

Run `python -m pytest tests` on Python 3.11–3.14 and build both wheel and sdist.
Tests include recovery of known charges from a synthetic Coulomb potential,
two-stage charge constraints, serialization and the Psi4 QCSchema boundary.
When changing scientific code, additionally run all stages in
`examples/parameterize.py` with external Psi4/AmberTools, verify generated atom
counts and charge sums, and compare fitted charges with a prior run.

External integration test (20 initial conformers per ligand, all four example
ligands and all pipeline stages; an isolated force-field copy is used):

```bash
LIGTOOLS_PSI4="$(command -v psi4)" \
LIGTOOLS_FORCEFIELD=/path/to/amber99sb-ildn.ff python -m pytest -q tests
```

The native optimization/thread checks require `LIGTOOLS_PSI4`; the full ligand
and Jupyter tests additionally require `LIGTOOLS_FORCEFIELD` and notebook test
dependencies. The tests copy the given force field into temporary directories.
