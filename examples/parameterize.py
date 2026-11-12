"""Parameterize four small carboxylate/phosphinate ligands for GROMACS.

Run it from a directory that contains the force field you want to extend::

    python example.py

Every step writes into `work_dir`, so an interrupted run can be resumed by
re-running the script.
"""

from ligtools import Smi2Gmx

# SMILES of the molecules to parameterize
list_smiles = [
    'CC(=O)C(=O)[O-]',
    'CC(=O)[PH](=O)[O-]',
    'CC(=O)P(OC)(=O)[O-]',
    'CC(=O)P(C)(=O)[O-]',
]

# GROMACS residue names, one per SMILES (3-4 characters is the usual convention)
list_names = ['PYR', 'API', 'MPO', 'MPI']

ligands = Smi2Gmx(
    smiles=list_smiles,
    res_names=list_names,
    work_dir='./output_param/',
    n_conformer_pool=4000,     # Initial candidates per molecule, before filtering
    n_conformer_threads=4,     # RDKit embedding threads
)

# Draw the input molecules (Jupyter/IPython only)
ligands.draw_smiles()

# 1. QM: geometry optimization + single points with Psi4
ligands.run_psi4()

# 2. Fit RESP2 charges under symmetry and net-charge constraints
ligands.run_resp()

# 3. Assign GAFF2 parameters with acpype/antechamber
ligands.run_acpype(parallel=False)

# Draw again, now annotated with the fitted partial charges
ligands.draw_smiles()

# 4. Collect the results
ligands.get_itps()
ligands.get_gros()
ligands.get_pdbs()
ligands.get_rtps()

# 5. Add the new residues to an existing GROMACS force field directory
ligands.get_ff_with_rtps(ff_path='./amber-99sb-ildn.ff')
