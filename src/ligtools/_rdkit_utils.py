"""Shared RDKit helpers used by more than one ligtools module."""

from rdkit import Chem


def simplify_molecule(mol, add_hs=True):
    """Return a purely topological copy of `mol`.

    Every bond is made single, aromaticity is cleared and formal charges are
    dropped, so that atoms differing only in bond order or charge (a
    carboxylate =O and [O-], for instance) become indistinguishable. This is
    what makes them rank as equivalent in `CanonicalRankAtoms` and match each
    other in `GetSubstructMatch`.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to flatten. It is not modified.
    add_hs : bool
        Add explicit hydrogens first. Use True when the caller works with
        hydrogens (charge fitting), False when matching a structure that
        already carries the hydrogens it needs.
    """
    if add_hs:
        mol = Chem.AddHs(mol)

    new_mol = Chem.RWMol(mol)

    for bond in list(new_mol.GetBonds()):
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if bond.GetBondType() != Chem.BondType.SINGLE:
            new_mol.RemoveBond(begin, end)
            new_mol.AddBond(begin, end, Chem.BondType.SINGLE)

    for atom in new_mol.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
        # Hydrogens are explicit here, so no implicit H must be inferred.
        # Without this, dropping the formal charge of e.g. a carboxylate [O-]
        # leaves it with 0 implicit H while its =O partner gains 1 (it becomes
        # an -OH), so the two oxygens are no longer interchangeable.
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)

    for bond in new_mol.GetBonds():
        bond.SetIsAromatic(False)

    result = new_mol.GetMol()
    # RDKit >= 2023 requires implicit valences and ring info to be computed
    # before CanonicalRankAtoms or substructure matching. A full SanitizeMol is
    # not an option: the graph is deliberately "wrong" chemically, so it would
    # raise on valence errors. strict=False computes what is needed without
    # validating valences.
    result.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(result)
    return result
