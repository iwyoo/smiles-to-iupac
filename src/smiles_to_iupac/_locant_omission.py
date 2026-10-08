"""Substituent locants that carry no information (P-14.3.4.3, P-14.3.4.5, P-14.3.4.6): a sole substituent on a
parent with one kind of substitutable hydrogen, all positions bearing the same substituent, or all hydrogens on one
carbon. The decision depends only on the parent structure, never on its numbering."""

from rdkit import Chem


def _is_substitutable(mol, atom, suffix_atoms):
    """The hydrogen of a formyl group or the carbon of a cyano group is not a substitutable position."""
    if atom.GetAtomicNum() != 6:
        return True
    for bond in atom.GetBonds():
        other = bond.GetOtherAtomIdx(atom.GetIdx())
        if other in suffix_atoms and (bond.GetBondTypeAsDouble() == 3.0 or (bond.GetBondTypeAsDouble() == 2.0 and atom.GetTotalNumHs())):
            return False
    return True


def _branches(mol, parent, suffix_atoms, free_atoms=frozenset()):
    return [
        (p, n.GetIdx())
        for p in parent
        for n in mol.GetAtomWithIdx(p).GetNeighbors()
        if n.GetIdx() not in parent and n.GetIdx() not in suffix_atoms and n.GetIdx() not in free_atoms
    ]


def _hydride_ranks(mol, parent, suffix_atoms):
    """Canonical symmetry classes of the parent hydride with its suffix groups, every branch replaced by hydrogen."""
    keep = set(parent) | set(suffix_atoms)
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        idx = atom.GetIdx()
        if idx not in keep:
            continue
        lost = sum(
            b.GetBondTypeAsDouble()
            for b in mol.GetAtomWithIdx(idx).GetBonds()
            if b.GetOtherAtomIdx(idx) not in keep
        )
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(mol.GetAtomWithIdx(idx).GetTotalNumHs() + int(lost))
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for idx in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    for bond in editable.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    fragment = editable.GetMol()
    try:
        Chem.SanitizeMol(fragment)
    except Exception:
        return None
    order = sorted(keep)
    ranks = Chem.CanonicalRankAtoms(fragment, breakTies=False)
    return {atom: ranks[i] for i, atom in enumerate(order)}, fragment, order


def omits_all_locants(mol, parent, grouped, suffix_atoms=(), single_kind=True, free_atoms=frozenset()):
    """True when the substituents of `grouped` are cited without locants.

    `parent`: atom indices of the parent skeleton; `suffix_atoms`: atoms of the
    principal characteristic groups on it; `free_atoms`: atoms the parent is attached to;
    `single_kind=False` limits the rule to full substitution (a ring carrying a suffix keeps its locants).
    """
    if mol is None or not grouped:
        return False
    parent = set(parent)
    suffix_atoms = set(suffix_atoms)
    branches = _branches(mol, parent, suffix_atoms, set(free_atoms))
    if not branches:
        return False
    if sum(len(info["locants"]) for info in grouped.values()) != len(branches):
        return False
    if len(grouped) == 1 and not any(
        mol.GetAtomWithIdx(p).GetTotalNumHs() and _is_substitutable(mol, mol.GetAtomWithIdx(p), suffix_atoms)
        for p in parent
    ):
        return True
    if not single_kind:
        return False
    resolved = _hydride_ranks(mol, parent, suffix_atoms)
    if resolved is None:
        return False
    ranks, fragment, order = resolved
    slots = [
        p
        for i, p in enumerate(order)
        if p in parent
        and fragment.GetAtomWithIdx(i).GetTotalNumHs()
        and _is_substitutable(mol, mol.GetAtomWithIdx(p), suffix_atoms)
    ]
    if len(slots) == 1:
        return mol.GetAtomWithIdx(slots[0]).GetAtomicNum() == 6
    return len(branches) == 1 and len({ranks[p] for p in slots}) == 1
