"""ligtools: prepare small-molecule ligands for GROMACS/AMBER simulations.

Two independent tools that are usually used one after the other:

- :class:`Smi2Gmx` turns SMILES into GROMACS topologies (`.itp`, `.rtp`) with
  RESP charges and GAFF2 parameters, via Psi4, psiresp and acpype.
- :class:`CoordMatcher` moves the resulting ligand into the coordinate frame of
  a reference structure, and rebuilds the complex around it.

::

    from ligtools import Smi2Gmx, CoordMatcher

The two are imported lazily, so using CoordMatcher does not require Psi4,
the parameterization backends - it only needs RDKit, ParmEd, NumPy and IPython.
"""

__version__ = '0.2.0'

__all__ = ['Smi2Gmx', 'CoordMatcher', 'DEFAULT_PSI4', '__version__']

_LAZY = {
    'Smi2Gmx': ('.smi2gmx', 'Smi2Gmx'),
    'DEFAULT_PSI4': ('.smi2gmx', 'DEFAULT_PSI4'),
    'CoordMatcher': ('.coordmatcher', 'CoordMatcher'),
}


def __getattr__(name):
    """Import Smi2Gmx / CoordMatcher on first use (PEP 562)."""
    try:
        module_name, attr = _LAZY[name]
    except KeyError:
        raise AttributeError(f'module {__name__!r} has no attribute {name!r}') from None

    from importlib import import_module

    value = getattr(import_module(module_name, __name__), attr)
    globals()[name] = value
    return value


def __dir__():
    return sorted(__all__)
