# ligtools

Ligand preparation for **GROMACS**/**AMBER** molecular dynamics: build the topology from a
SMILES string, then put the molecule where it belongs in your system.

Two independent tools:

| Tool | What it does | Needs |
| --- | --- | --- |
| **`Smi2Gmx`** | SMILES -> GAFF2 topology with RESP2 charges (`.itp`, `.gro`, `.pdb`, `.rtp`) | Psi4, psiresp, AmberTools, acpype |
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
SMILES -> RDKit/psiresp conformers -> Psi4 (QM) -> RESP2 charges -> acpype/GAFF2 -> .itp / .gro / .pdb / .rtp
```

The result is a set of GROMACS topology files, and optionally a force field directory with the
new residues already injected (`.rtp` entries, new atom types in `atomtypes.atp` and
`ffnonbonded.itp`).

---

## Features

- Conformer generation and geometry optimization driven by [psiresp](https://psiresp.readthedocs.io) + [Psi4](https://psicode.org)
- **RESP2** charge fitting with automatic detection of symmetry-equivalent atoms
  (topologically equivalent atoms get identical charges)
- Net-charge constraints per molecule, or fully custom charge-equivalence groups
- GAFF2 atom types and bonded parameters via [acpype](https://github.com/alanwilter/acpype) / AmberTools
- Output as GROMACS `.itp`, `.gro`, `.pdb` and — unlike plain acpype — **`.rtp`** residue
  entries, so ligands can be built by `gmx pdb2gmx` like any standard residue
- Optional injection of the new residues and atom types into an existing `*.ff` directory
- Parallel execution of the QM step (GNU `parallel`) and of the acpype step (`multiprocessing`)

---

## Requirements

Two of the dependencies are **not** available on PyPI and must be installed with conda:

| Dependency | Purpose | Source |
| --- | --- | --- |
| Psi4 | QM optimization and single points | `conda-forge` |
| AmberTools (`antechamber`, `parmchk2`) | GAFF2 typing, used by acpype | `conda-forge` |
| GNU `parallel` | only for `run_psi4(parallel=True)` | `conda-forge` / system package manager |

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

The RESP backend uses Pydantic 2 and QCElemental's explicit v2 models, NumPy
`isin`, and SciPy's public sparse solver. QCSchema v1 conversion is kept at the
Psi4 file boundary. Pydantic 1 and the old `qcelemental<0.50` constraint are no
longer required. ACPYPE resolves AmberTools from the active environment.

See [THIRD_PARTY.md](THIRD_PARTY.md) for provenance, licenses, and maintenance
instructions. These changes preserve the existing scientific configuration;
they do not change the charge-fitting method or its defaults.

---

## Installation

### Recommended: conda environment + pip

```bash
# 1. Install the external quantum chemistry and topology tools
conda create -n ligtools -c conda-forge python=3.13 'psi4>=1.11,<1.12' \
    ambertools openbabel parallel pip
conda activate ligtools

# 2. Install ligtools itself
pip install git+https://github.com/Vetrov-Anton/ligtools.git
```

`mamba` (or Miniforge, which ships it) is strongly recommended over plain `conda` here:
the environment is ~3.6 GB and 245 packages, and the classic solver is slow on it.

### One-liner with the bundled environment file

```bash
git clone https://github.com/Vetrov-Anton/ligtools.git
cd ligtools
conda env create -f environment.yml     # installs the dependencies *and* ligtools
conda activate ligtools
```

### Development install

```bash
git clone https://github.com/Vetrov-Anton/ligtools.git
cd ligtools
conda activate ligtools
pip install -e .
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

## Required: single-threaded BLAS

```bash
export OMP_NUM_THREADS=1
```

**Without this, ligtools deadlocks.** Psi4/MKL start OpenMP threads when they are imported;
psiresp then forks worker processes, and each child inherits a mutex that was held at fork time
by a thread that does not exist in the child. Every worker blocks on `futex` forever, at 0% CPU,
with no error message.

This costs nothing in practice: the QM work runs in separate `psi4` processes, and ligtools
parallelizes across conformers with `run_psi4(parallel=True)` rather than across threads.

## Quick start

```python
from ligtools import Smi2Gmx

ligands = Smi2Gmx(
    smiles=['CC(=O)C(=O)[O-]', 'CC(=O)P(C)(=O)[O-]'],
    res_names=['PYR', 'MPI'],
    work_dir='./output_param/',
)

ligands.run_psi4()                  # 1. QM optimization + single points
ligands.run_resp()                  # 2. RESP2 charge fit
ligands.run_acpype(parallel=False)  # 3. GAFF2 parameters

ligands.get_itps()                  # 4. GROMACS topologies with the fitted charges
ligands.get_gros()
ligands.get_pdbs()
ligands.get_rtps()

# 5. optional: add the residues to an existing force field
ligands.get_ff_with_rtps(ff_path='./amber99sb-ildn.ff')
```

A complete, commented script is in [`examples/parameterize.py`](examples/parameterize.py).

> The QM step is the expensive one. Four small ligands with 5 conformers each take on the order
> of hours on a workstation; plan accordingly and prefer `parallel=True` on a multi-core machine.

---

## How it works

| Step | Method | What happens |
| --- | --- | --- |
| 1 | `run_psi4()` | psiresp generates conformers, Psi4 optimizes them and computes the ESP grid. Because psiresp stops and writes shell scripts instead of blocking, `ligtools` rewrites `run_optimization.sh` / `run_single_point.sh` with absolute paths and executes them (through GNU `parallel` when `parallel=True`). |
| 2 | `run_resp()` | Symmetry-equivalent atoms are found with `Chem.CanonicalRankAtoms` on a "flattened" copy of the molecule (all bonds single, no formal charges), so e.g. the three methyl hydrogens or the two carboxylate oxygens share one charge. Net charge is constrained per molecule. |
| 3 | `run_acpype()` | Each molecule is written as a PDB with the requested residue name and passed to acpype (`chargeType='user'`), which derives GAFF2 atom types and bonded parameters. |
| 4 | `get_*()` | Results are collected; `get_itps()` overwrites acpype's charges with the RESP2 ones. |
| 5 | `get_ff_with_rtps()` | New atom types are appended to `atomtypes.atp` and `ffnonbonded.itp`, and one `.rtp` file per residue is copied into the force field directory. |

### Layout of `work_dir`

```
work_dir/
├── optimization/     # psiresp/Psi4 geometry optimization
├── single_point/     # psiresp/Psi4 single points
├── temporary_pdbs/   # input PDBs handed to acpype
├── acpype/           # <RES>.acpype/ directories produced by acpype
├── itps/             # <RES>_GMX.itp  (RESP2 charges)
├── gros/             # <RES>_GMX.gro
├── pdbs/             # <RES>.pdb
└── rtps/             # <RES>.rtp
```

Because every step writes to disk, an interrupted run can be resumed by re-running the script:
finished stages are picked up from `work_dir` instead of being recomputed.

---

## API

### `Smi2Gmx(...)`

| Argument | Default | Meaning |
| --- | --- | --- |
| `smiles` | `[]` | List of SMILES strings to parameterize. |
| `res_names` | `[]` | GROMACS residue name per molecule (required). |
| `max_conf_n` | `5` each | Maximum number of conformers per molecule. |
| `n_conformer_pool` | `4000` | Initial conformer candidates per molecule, before filtering; a positive integer or a list with one positive integer per SMILES. |
| `n_conformer_threads` | `1` | Positive number of RDKit threads for conformer embedding and optional force-field minimization. |
| `chiral` | `False` each | Pass acpype's chirality flag for that molecule. |
| `net_charges` | from SMILES | Total charge per molecule; overrides the value RDKit infers. |
| `ff_path` | `''` | Force field directory used by `get_ff_with_rtps()` when it is called without arguments. |
| `atomtypes_filename` | `''` | Stored for reference (see note below). |
| `ffnonbonded_filename` | `''` | Stored for reference (see note below). |
| `work_dir` | `./tmp_Smi2Gmx` | Working directory for all intermediate and output files. |
| `max_n_cores` | `32` | Upper bound on parallel workers for the QM and acpype steps. |

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

Embedding can produce fewer candidates because RDKit prunes similar geometries
or cannot embed a requested conformation. Molecules are processed in sequence;
`n_conformer_threads` controls RDKit threads within each molecule and requires a
build of RDKit with thread support. Energy and RMSD selection still run serially,
so increasing this setting does not accelerate every part of conformer selection.

### Methods

| Method | Description |
| --- | --- |
| `draw_smiles(charge=True, index=True, molsPerRow=3, subImgSize=(500, 500))` | Draw the molecules in a Jupyter notebook; after `run_resp()` the atoms are annotated with their partial charges and indices. |
| `run_psi4(psi4_path=DEFAULT_PSI4, parallel=True)` | Run the QM stage. `psi4_path` may be a bare executable name (looked up on `$PATH`) or a full path. |
| `run_resp(constraints=[])` | Fit RESP2 charges. Per molecule pass `'auto'` (symmetry detected automatically, the default), `'no'` (no equivalence constraints) or a list of index tuples for custom groups. |
| `run_acpype(parallel=True)` | Run acpype for every molecule. |
| `get_itps(dir='')` | Write `.itp` files with the RESP2 charges substituted in. |
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

Defaults to mode 3 when no SMILES were given, mode 4 otherwise. Modes 0, 1 and 3 preserve or
rigidly move the probe geometry; modes 2 and 4 replace it with a generated conformer, so expect
the bond lengths and angles to come from RDKit's force field rather than from your topology.

### Methods

| Method | Description |
| --- | --- |
| `align(modes=[], n_confs=1000, n_cores=16, opt_param=0.5, opt_geom=False, seed=80, addHs_to_probe=[], addHs_to_ref=[], outputs=None, draw=True)` | Align every probe onto its reference. Results are collected in `self.pmd_objs`; pass `outputs` to also write them. `draw=False` skips the RDKit depictions (use it outside Jupyter). |
| `get_pdbs(filenames)` | Write the aligned structures to the given file names. |
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
