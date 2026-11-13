# Local optimization report

[Russian PDF: previous and corrected psiresp optimization](psiresp_local_optimization_changes_ru.pdf)
describes the single-gradient bug, the `psi4.optimize` implementation, optional
QM geometry optimization, current convergence defaults and validation results.
The default parameterization workflow uses MMFF94 geometries and Psi4 single
points. Historical stopped TPP optimization attempts are distinguished from the
successful MMFF94-based parameterization.

To rebuild the PDF on Linux with DejaVu fonts installed:

```bash
python -m pip install reportlab
python reports/build_optimization_report.py
```

ReportLab is needed only to build this report, not to run ligtools.
