"""Place a parameterized ligand into a binding site and rebuild the complex.

Run this after `parameterize.py` has produced the ligand PDB.
Remember to `export OMP_NUM_THREADS=1` (see README).
"""

from ligtools import CoordMatcher

cm = CoordMatcher(
    # The ligand as it comes out of Smi2Gmx.get_pdbs(): correct topology,
    # arbitrary orientation.
    xyz_probes=['output_param/pdbs/PYR.pdb'],
    # The pose to reproduce: a crystal ligand, a docked pose, ...
    xyz_refs=['site_ligand.pdb'],
    # Only needed for the conformer-searching modes (2 and 4).
    smiles=['CC(=O)C(=O)[O-]'],
)

# mode 1 = rigid fit over the maximum common substructure, keeps the geometry
# that came with the topology. Use draw=False outside Jupyter.
cm.align(modes=[1], outputs=['PYR_aligned.pdb'], draw=False)

# Replace the old ligand copy in the complex with the aligned one, keeping
# water and ions.
cm.get_complex(
    output='new_start.pdb',
    scaffold='complex.pdb',
    pdb_names=['PYR'],
)
