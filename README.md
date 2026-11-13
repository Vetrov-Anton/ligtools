# ligtools

Ligand preparation for **GROMACS**/**AMBER** molecular dynamics: build the topology from a
SMILES string, then put the molecule where it belongs in your system.

Two independent tools:

| Tool | What it does | Needs |
| --- | --- | --- |
| **`Smi2Gmx`** | SMILES -> GAFF2 topology with RESP charges (`.itp`, `.gro`, `.pdb`, `.rtp`) | Psi4, psiresp, AmberTools, acpype |
| **`CoordMatcher`** | Move a ligand onto a reference pose and rebuild the complex | RDKit, ParmEd, NumPy only |

```python
from ligtools import Smi2Gmx, CoordMatcher
```

The two are imported lazily, so `CoordMatcher` works in a minimal environment - importing it
does not pull in Psi4, psiresp or acpype.

## Smi2Gmx

`Smi2Gmx` chains together the standard tools of the parameterization workflow so that you do not
have to run them by hand:

```
SMILES -> RDKit/MMFF94 -> conformer selection -> Psi4 single points -> RESP -> acpype/GAFF2 -> .itp / .gro / .pdb / .rtp
```

The result is a set of GROMACS topology files, and optionally a force field directory with the
new residues already injected (`.rtp` entries, new atom types in `atomtypes.atp` and
`ffnonbonded.itp`).

---

## Features

- RDKit/MMFF94 conformers, Psi4 single points and optional QM geometry optimization via the maintained psiresp backend
- **RESP** charge fitting with automatic detection of symmetry-equivalent atoms
  (topologically equivalent atoms get identical charges)
- Net-charge constraints per molecule, or fully custom charge-equivalence groups
- GAFF2 atom types and bonded parameters via [acpype](https://github.com/alanwilter/acpype) / AmberTools
- Output as GROMACS `.itp`, `.gro`, `.pdb` and — unlike plain acpype — **`.rtp`** residue
  entries, so ligands can be built by `gmx pdb2gmx` like any standard residue
- Optional injection of the new residues and atom types into an existing `*.ff` directory
- Parallel QM execution with Python’s standard library and parallel ACPYPE workers

---

## Requirements

Two of the dependencies are **not** available on PyPI and must be installed with conda:

| Dependency | Purpose | Source |
| --- | --- | --- |
| Psi4 | QM optimization and single points | `conda-forge` |
| AmberTools (`antechamber`, `parmchk2`) | GAFF2 typing, used by acpype | `conda-forge` |

The maintained Python backends are included as `ligtools.psiresp` and
`ligtools.acpype`; do not install upstream `psiresp` or `acpype` separately.
Other Python dependencies are installed by pip. Open Babel is also needed for
PDB-to-MOL2 conversion in the parameterization workflow (install with conda).

**Python 3.11–3.14.** Unit and numerical regression tests cover this range.
Availability of Psi4/AmberTools binaries additionally depends on the platform.

### Maintained backends

The code originally released as psiresp 0.4.2 and acpype 2023.10.27 now lives in
`src/ligtools/psiresp/` and `src/ligtools/acpype/`. It is built and versioned with
ligtools, without upstream setup.py, Versioneer, or bundled Amber executables.

Smi2Gmx currently uses gas-phase HF/6-31G* and two-stage RESP. The RESP2
vacuum/solvent mixing protocol is not enabled by this interface.

The RESP backend uses Pydantic 2 and QCElemental's explicit v2 models, NumPy
`isin`, and SciPy's public sparse solver. QCSchema v1 conversion is kept at the
Psi4 single-point file boundary; local optimization uses QCSchema v2.
Pydantic 1 and the old `qcelemental<0.50` constraint are no
longer required. ACPYPE resolves AmberTools from the active environment.

By default, Smi2Gmx uses the selected RDKit/MMFF94 geometries directly for Psi4
single-point wavefunction and ESP calculations. It does not run QM geometry
optimization. To enable it explicitly, construct
`Smi2Gmx(..., optimize_geometry=True)`. In that mode, `psi4.optimize` must converge
before the final coordinates are passed on; failure stops the pipeline. Old
single-gradient optimization caches are excluded. Restart the notebook kernel
and create a fresh `Smi2Gmx` instance after upgrading.

When QM geometry optimization is enabled, it uses explicit maximum/RMS force and displacement
thresholds of `1e-3`, in OptKing's internal-coordinate atomic units. Length
coordinates use hartree/bohr and bohr; angular coordinates use hartree/radian and
radians. `QMGeometryOptimizationOptions(convergence_tolerance=None)` restores
the named `g_convergence` preset (default `gau_tight`). Individual explicit
`*_g_convergence` Psi4 keywords override the common tolerance. The iteration
limit and Hessian frequency are separate settings. The energy-change criterion
remains inactive by default.

See [THIRD_PARTY.md](THIRD_PARTY.md) for provenance, licenses, and maintenance
instructions. The charge-fitting method remains two-stage RESP; by default it
now uses MMFF94 geometries instead of QM-optimized geometries.

---

## Installation

### Recommended: conda environment + pip

```bash
# 1. Install the external quantum chemistry and topology tools
conda create -n ligtools -c conda-forge python=3.13 'psi4>=1.11,<1.12' \
    ambertools openbabel pip
conda activate ligtools

# 2. Install ligtools itself
python -m pip install git+https://github.com/Vetrov-Anton/ligtools.git

# 3. For Jupyter, install and register the environment's kernel
python -m pip install ipykernel
python -m ipykernel install --user --name ligtools --display-name "Python (ligtools)"
```

Select **Python (ligtools)** in your notebook. If you only use standalone Python
scripts, step 3 is optional. No BLAS environment settings are needed.

`mamba` (or Miniforge, which ships it) is strongly recommended over plain `conda` here:
the environment is ~3.6 GB and 245 packages, and the classic solver is slow on it.

### One-liner with the bundled environment file

```bash
git clone https://github.com/Vetrov-Anton/ligtools.git
cd ligtools
conda env create -f environment.yml     # installs the dependencies *and* ligtools
conda activate ligtools
```

This environment file already includes `ipykernel`. If Jupyter runs from another
environment, register this kernel using the `python -m ipykernel install` command
above and select **Python (ligtools)**.

Both installation paths were tested from fresh environments on Linux with
Python 3.13.18, Psi4 1.11 and AmberTools 26, including the four-ligand notebook
pipeline. The tested AmberTools distribution also includes auxiliary Python
tools with incompatible or missing dependencies reported by `pip check` (for
example, `numpy<2` requirements). Those auxiliary tools are not used by ligtools;
the successful pipeline check does not establish that every AmberTools utility works. Details
are recorded in [VALIDATION.md](VALIDATION.md).

### Development install

```bash
git clone https://github.com/Vetrov-Anton/ligtools.git
cd ligtools
conda activate ligtools
python -m pip install -e .
```

### Check the installation

```bash
python -c "import ligtools; print(ligtools.__version__)"
which psi4 antechamber parmchk2
```

To check the included backends:

```bash
python -c "from ligtools import Smi2Gmx, CoordMatcher; from ligtools import psiresp, acpype"
```

If `psi4` is not on your `$PATH`, either pass its location explicitly
(`ligands.run_psi4(psi4_path='/path/to/psi4')`) or export `PSI4_EXE=/path/to/psi4`.

---

## Jupyter and parallel execution

Run the API directly in notebook cells; no `OMP_NUM_THREADS`, BLAS settings or
multiprocessing configuration is required. ESP and parallel ACPYPE workers start
fresh Python interpreters with the local `spawn` context. Their BLAS/OpenMP thread
limits are applied inside the workers, without changing the notebook environment
or the thread count used by RDKit conformer generation.

`max_n_cores` limits ESP/ACPYPE workers and the total QM thread budget.
`run_psi4(n_threads=8)` uses eight threads per Psi4 process; set
`max_n_cores=8` to run one such process at a time. With `parallel=True`, at most
`max_n_cores // n_threads` QM processes run concurrently. `parallel=False`
retains the requested thread count and runs one process at a time. The default
is one thread per Psi4 process. Existing ESP
results are reused without starting a worker pool. A notebook example is in
[`examples/parameterize.ipynb`](examples/parameterize.ipynb).

For standalone Python scripts, put the calls inside a function and invoke it
under `if __name__ == "__main__":`, as in `examples/parameterize.py`. This guard
is needed for spawn in scripts; notebook cells do not need it.

## Quick start

```python
from ligtools import Smi2Gmx

ligands = Smi2Gmx(
    smiles=['CC(=O)C(=O)[O-]', 'CC(=O)P(C)(=O)[O-]'],
    res_names=['PYR', 'MPI'],
    work_dir='./output_param/',
    n_conformer_threads=8,          # RDKit/MMFF94/pruning worker limit
    max_n_cores=8,                  # Total QM thread budget
    optimize_geometry=False,       # Default: keep the MMFF94 geometries
)

ligands.run_psi4(n_threads=8)       # 1. Eight threads per QM single point
ligands.run_resp()                  # 2. RESP charge fit
ligands.run_acpype(parallel=False)  # 3. GAFF2 parameters

ligands.get_itps()                  # 4. GROMACS topologies with the fitted charges
ligands.get_gros()
ligands.get_pdbs()
ligands.get_rtps()

# 5. optional: add residues to a COPY of your force field (modifies that directory)
ligands.get_ff_with_rtps(ff_path='./amber99sb-ildn.ff')
```

A complete, commented script is in [`examples/parameterize.py`](examples/parameterize.py).

The MMFF94 workflow was checked through force-field RTP export for API, MPO and
MPI (44.2 seconds total on the audit workstation) and for TPP with net charge -2.
Runtime depends on molecule size and selected conformer count. Optional QM
geometry optimization can take substantially longer; see [VALIDATION.md](VALIDATION.md).

---

## How it works

| Step | Method | What happens |
| --- | --- | --- |
| 1 | `run_psi4()` | Generate and select MMFF94 conformers, then compute Psi4 orientation wavefunctions at fixed geometry. QM optimization runs only with `optimize_geometry=True`. Python schedules up to `max_n_cores // n_threads` processes at once. Each calculation has a `.psi4.log`; failed exits/results raise errors. |
| 2 | `run_resp()` | Symmetry-equivalent atoms are found with `Chem.CanonicalRankAtoms` on a "flattened" copy of the molecule (all bonds single, no formal charges), so e.g. the three methyl hydrogens or the two carboxylate oxygens share one charge. Net charge is constrained per molecule. |
| 3 | `run_acpype()` | Each molecule is written as a PDB with the requested residue name and passed to acpype for GAFF2 parameters. With PDB input, ACPYPE falls back from `chargeType='user'` to an internal AM1-BCC calculation. Final exported ITP/RTP charges are replaced with RESP charges. |
| 4 | `get_*()` | Results are collected; `get_itps()` overwrites acpype's charges with the RESP ones. |
| 5 | `get_ff_with_rtps()` | New atom types are appended to `atomtypes.atp` and `ffnonbonded.itp`, and one `.rtp` file per residue is copied into the force field directory. |

### Layout of `work_dir`

```
work_dir/
├── optimization/     # Present only when QM geometry optimization is enabled
├── single_point/     # psiresp/Psi4 single points
├── temporary_pdbs/   # input PDBs handed to acpype
├── acpype/           # <RES>.acpype/ directories produced by acpype
├── itps/             # <RES>_GMX.itp  (RESP charges)
├── gros/             # <RES>_GMX.gro
├── pdbs/             # <RES>.pdb
└── rtps/             # <RES>.rtp
```

Rerunning with the same options and `work_dir` regenerates conformers and reuses
matching completed QM files. ESP/RESP and exports may run again. An interrupted
geometry optimization restarts from its input; there is no automatic continuation
from its last geometry step. Use a new work directory for an independent run.

---

## API

### `Smi2Gmx(...)`

| Argument | Default | Meaning |
| --- | --- | --- |
| `smiles` | `[]` | List of SMILES strings to parameterize. |
| `res_names` | `[]` | GROMACS residue name per molecule (required). |
| `max_conf_n` | `5` each | Maximum number of conformers per molecule. |
| `n_conformer_pool` | `4000` | Initial conformer candidates per molecule, before filtering; a positive integer or a list with one positive integer per SMILES. |
| `n_conformer_threads` | `1` | Positive worker limit for embedding, force-field minimization, energy scoring and RMSD selection. |
| `optimize_geometry` | `False` | Enable optional Psi4 geometry optimization after MMFF94 selection. The default uses MMFF94 coordinates directly for QM single points and ESP. |
| `chiral` | `False` each | Pass acpype's chirality flag for that molecule. |
| `net_charges` | from SMILES | Total charge per molecule; overrides the value RDKit infers. |
| `ff_path` | `''` | Force field directory used by `get_ff_with_rtps()` when it is called without arguments. |
| `atomtypes_filename` | `''` | Stored for reference (see note below). |
| `ffnonbonded_filename` | `''` | Stored for reference (see note below). |
| `work_dir` | `./tmp_Smi2Gmx` | Working directory for all intermediate and output files. |
| `max_n_cores` | `32` | Total thread budget for QM; upper bound on parallel workers for ESP and ACPYPE. |

All path arguments accept `str` or `pathlib.Path`, relative or absolute, with `~` expansion;
they are resolved to absolute paths on construction.

For example, generate up to 200 initial candidates for each molecule on four
RDKit threads, then keep at most five per molecule for QM:

```python
ligands = Smi2Gmx(
    smiles=['CC(=O)C(=O)[O-]', 'CC(=O)P(C)(=O)[O-]'],
    res_names=['PYR', 'MPI'],
    n_conformer_pool=200,       # Or [200, 500] for separate pool sizes
    n_conformer_threads=4,
    max_conf_n=[5, 5],
    max_n_cores=2,             # Separate limit for QM/acpype workers
)
```

Embedding can produce fewer candidates because similar geometries are pruned
or RDKit cannot embed a requested conformation. Generated candidates are minimized
with MMFF94 by default before energy-window and RMSD selection. The selected
coordinates are used directly for Psi4 single points; subsequent QM geometry
optimization is optional (`optimize_geometry=True`). MMFF94 minimization applies
to the psiresp pipeline; CoordMatcher does not perform force-field minimization.
Molecules are processed in sequence;
`n_conformer_threads` controls RDKit threads within each molecule and requires a
build of RDKit with thread support. Energy scoring uses bounded blocks and
threads. RMSD selection computes only distances to already selected conformers,
using spawned processes for larger candidate sets and native RDKit threads for
small sets. This avoids a full pairwise distance matrix. Notebook users do not
need thread environment variables. Embedding-time RMSD pruning is performed
separately in bounded parallel batches, using RDKit's atom-symmetry mappings,
strict threshold and original acceptance order. A rotation/permutation-invariant
lower bound skips provably distant pairs before exact alignment. The greedy
acceptance decisions remain sequential; a worker limit does not guarantee
continuous use of every core. With only one retained conformer, only one Psi4
process can run; `run_psi4(n_threads=8)` lets that process use eight threads
when `max_n_cores` is at least eight.

### Methods

| Method | Description |
| --- | --- |
| `draw_smiles(charge=True, index=True, molsPerRow=3, subImgSize=(500, 500))` | Draw the molecules in a Jupyter notebook; after `run_resp()` the atoms are annotated with their partial charges and indices. |
| `run_psi4(psi4_path=DEFAULT_PSI4, parallel=True, *, n_threads=1)` | Run the QM stage with a positive integer thread count per Psi4 process, no greater than `max_n_cores`. `psi4_path` may be a bare executable name or a full path. |
| `run_resp(constraints=[])` | Fit RESP charges. Per molecule pass `'auto'` (symmetry detected automatically, the default), `'no'` (no equivalence constraints) or a list of index tuples for custom groups. |
| `run_acpype(parallel=True)` | Run acpype for every molecule. |
| `get_itps(dir='')` | Write `.itp` files with the RESP charges substituted in. |
| `get_gros(dir='')` | Copy acpype's `.gro` files. |
| `get_pdbs(dir='')` | Convert the `.gro` files to `.pdb`. |
| `get_rtps(dir='')` | Write GROMACS `.rtp` residue entries (bonds, angles, dihedrals, impropers converted to GROMACS units). |
| `get_ff_with_rtps(ff_path='', atomtypes_filename='atomtypes.atp', ffnonbonded_filename='ffnonbonded.itp')` | Append the new atom types to the force field and copy the `.rtp` files into it. Calls `get_rtps()` first if needed. |

Each `get_*` method writes into `work_dir/<name>/` unless another directory is passed.

---

## CoordMatcher

A freshly parameterized ligand comes out of acpype in an arbitrary orientation. `CoordMatcher`
moves it onto a reference pose - a crystal ligand, a docked pose, the same ligand from an earlier
simulation - and can rebuild the whole complex around it.

```python
from ligtools import CoordMatcher

cm = CoordMatcher(
    xyz_probes=['PYR.pdb'],            # the ligand to move
    xyz_refs=['site_ligand.pdb'],      # where it should end up
    smiles=['CC(=O)C(=O)[O-]'],        # only needed for modes 2 and 4
)

cm.align(modes=[1], outputs=['PYR_aligned.pdb'])

# Drop the old ligand copy from the complex and insert the aligned one,
# keeping water and ions
cm.get_complex(output='new_start.pdb', scaffold='complex.pdb', pdb_names=['PYR'])
```

### Alignment modes

| mode | name | Strategy | Needs SMILES |
| --- | --- | --- | --- |
| 0 | `EXACT` | Copy reference coordinates atom by atom. The probe must be a substructure of the reference, and the result is exact. | no |
| 1 | `ONE_CONF` | Rigid-body fit over the maximum common substructure (`rdFMCS`). Keeps the probe conformer. | no |
| 2 | `CONFS` | Generate conformers from SMILES, keep the one with the lowest MCS RMSD to the reference. | yes |
| 3 | `ONE_SHAPE` | Rigid-body fit maximizing shape overlap (`rdShapeAlign`). Useful when probe and reference are different molecules. | no |
| 4 | `SHAPES` | Generate conformers from SMILES, keep the best shape + color overlap. | yes |

Defaults to mode 3 when no SMILES were given, mode 4 otherwise. Mode 0 copies the reference
geometry; modes 1 and 3 rigidly move the probe; modes 2 and 4 replace it with a generated conformer, so expect
the bond lengths and angles to come from RDKit's conformer embedding rather than from your topology.

All modes preserve the original probe's atom count, order, names and residues in
the output. Existing hydrogens are retained; there is no option to add hydrogens
to the input PDB or to optimize its geometry. The former `addHs_to_probe`,
`addHs_to_ref` and `opt_geom` parameters have been removed. Modes 2 and 4 use
explicit hydrogens only in the temporary SMILES molecule used by RDKit; only
atoms mapped to the original probe are written back.

MCS fits use heavy atoms. SMILES mapping first respects chemical bond orders
and falls back to topology for PDBs without bond-order information. Symmetric
atoms can have multiple equivalent mappings. Invalid inputs or missing matches
raise an error instead of silently producing an incomplete structure. `outputs`
is supported in all five modes.

### Methods

| Method | Description |
| --- | --- |
| `align(modes=None, n_confs=1000, n_cores=16, opt_param=0.5, seed=80, outputs=None, draw=False)` | Align every probe onto its reference. Results are collected in `self.pmd_objs`; pass `outputs` to also write them. Drawing is disabled by default; pass `draw=True` to show RDKit depictions. |
| `get_pdbs(filenames)` | Write the aligned structures to the given file names. |
| `get_pmd()` | Return a list of ParmEd structures: original probes before alignment, aligned results afterwards. Originals are loaded on first access. Structure objects are returned by reference; the list itself is separate. |
| `get_complex(output='new_start.pdb', scaffold='', pdb_names=[], water=[], ions=[])` | Rebuild a complex: scaffold minus the residues named in `pdb_names`, plus the aligned ligands, plus ions and water. Extra residue names can be declared with `water=` and `ions=`; `HOH` and `NA/CL/CA/MG` are always included. |

## Notes and caveats

- **`get_ff_with_rtps()` modifies the force field in place.** It appends to `atomtypes.atp` and
  rewrites `ffnonbonded.itp`. Copy your `*.ff` directory before calling it, and do not call it
  twice with the same residues unless you want duplicate entries.
- The `atomtypes_filename` / `ffnonbonded_filename` constructor arguments are stored but not read
  by `get_ff_with_rtps()`; pass the file names to the method itself if they differ from the
  GROMACS defaults.
- `run_acpype()` temporarily changes the process working directory (acpype writes next to the
  current directory); it is restored afterwards, including on failure.
- `psiresp.grid.GridOptions.vdw_radii` is set at import time to add a van der Waals radius for
  bromine. Edit the top of `src/ligtools/smi2gmx.py` if you need radii for other elements.
- Always sanity-check the generated parameters (charges, missing dihedrals reported by
  `parmchk2`, and the geometry after a short minimization) before running production MD.

---

## Acknowledgements

`ligtools` is a wrapper. The scientific work is done by
[Psi4](https://psicode.org), [psiresp](https://psiresp.readthedocs.io),
[RDKit](https://www.rdkit.org), [AmberTools](https://ambermd.org/AmberTools.php),
[acpype](https://github.com/alanwilter/acpype), [ParmEd](https://parmed.github.io/ParmEd/)
and [MDAnalysis](https://www.mdanalysis.org). Please cite them if you use this tool in
published work.

## License

Original ligtools code: MIT — see [LICENSE](LICENSE).
Included psiresp code: LGPL-3.0; included acpype code: GPL-3.0.
The combined distribution is GPL-3.0; see [THIRD_PARTY.md](THIRD_PARTY.md).
