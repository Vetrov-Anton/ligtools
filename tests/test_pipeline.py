"""Opt-in external-tool integration: LIGTOOLS_FORCEFIELD=/path/to/amber99sb-ildn.ff."""
import os
from pathlib import Path
import shutil

import numpy as np
import parmed as pmd
import pytest

from ligtools import Smi2Gmx


@pytest.mark.skipif(not os.environ.get('LIGTOOLS_FORCEFIELD'), reason='requires Psi4, AmberTools and force field')
def test_four_ligand_pipeline(tmp_path):
    forcefield = tmp_path / 'amber99sb-ildn.ff'
    shutil.copytree(Path(os.environ['LIGTOOLS_FORCEFIELD']), forcefield)
    names = ['PYR', 'API', 'MPO', 'MPI']
    ligands = Smi2Gmx(
        smiles=['CC(=O)C(=O)[O-]', 'CC(=O)[PH](=O)[O-]',
                'CC(=O)P(OC)(=O)[O-]', 'CC(=O)P(C)(=O)[O-]'],
        res_names=names, work_dir=tmp_path / 'output', max_n_cores=2,
        n_conformer_pool=20, n_conformer_threads=2,
    )
    ligands.run_psi4()
    ligands.run_resp()
    ligands.run_acpype(parallel=False)
    ligands.get_itps()
    ligands.get_gros()
    ligands.get_pdbs()
    ligands.get_rtps()
    ligands.get_ff_with_rtps(ff_path=forcefield)
    for name, molecule in zip(names, ligands.mols):
        root = tmp_path / 'output'
        top = pmd.load_file(str(root / 'itps' / f'{name}_GMX.itp'), parametrize=False)
        charges = np.array([atom.charge for atom in top.atoms])
        assert len(top.atoms) == molecule.n_atoms
        assert np.isfinite(charges).all()
        assert charges.sum() == pytest.approx(-1, abs=1e-5)
        np.testing.assert_allclose(charges, molecule.charges, atol=1e-6)
        for kind, filename in [('gros', f'{name}_GMX.gro'), ('pdbs', f'{name}.pdb')]:
            structure = pmd.load_file(str(root / kind / filename))
            assert len(structure.atoms) == molecule.n_atoms
            assert np.isfinite(structure.coordinates).all()
        assert (root / 'rtps' / f'{name}.rtp').read_bytes() == (forcefield / f'{name}.rtp').read_bytes()
