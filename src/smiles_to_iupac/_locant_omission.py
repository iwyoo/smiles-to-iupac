"""Substituent locants that carry no information (P-14.3.4.3 to P-14.3.4.6): a sole substituent on a parent with one
kind of substitutable hydrogen, all positions bearing the same substituent, all hydrogens on one carbon, or groups
that no other arrangement on the parent hydride could place differently. The decision depends only on the parent
structure, never on its numbering."""

from contextvars import ContextVar
from itertools import combinations

from rdkit import Chem

CITE_SUBSTITUENT_LOCANTS = ContextVar("cite_substituent_locants", default=False)


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
    if mol is None or not grouped or CITE_SUBSTITUENT_LOCANTS.get():
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
    slots, ranks = _hydrogen_slots(mol, parent, suffix_atoms)
    if len(slots) == 1:
        return mol.GetAtomWithIdx(slots[0]).GetAtomicNum() == 6
    return len(branches) == 1 and len({ranks[p] for p in slots}) == 1


def omits_group_locants(mol, parent, count):
    """True when the `count` identical principal groups, the only groups on `parent`, need no locant: a sole group on a
    parent with one kind of substitutable hydrogen (P-14.3.4.3) or groups on every substitutable position (P-14.3.4.5)."""
    if mol is None or CITE_SUBSTITUENT_LOCANTS.get():
        return False
    parent = set(parent)
    branches = _branches(mol, parent, frozenset())
    if len(branches) != count or any(mol.GetBondBetweenAtoms(*b).GetBondTypeAsDouble() != 1.0 for b in branches):
        return False
    if any(b.GetBondTypeAsDouble() > 1.0 and not b.GetIsAromatic() for b in mol.GetBonds() if {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= parent):
        return False
    if count > 1:
        return not any(mol.GetAtomWithIdx(p).GetTotalNumHs() for p in parent)
    slots, ranks = _hydrogen_slots(mol, parent, frozenset())
    return len({ranks[p] for p in slots}) == 1 and (len(slots) > 1 or mol.GetAtomWithIdx(slots[0]).GetAtomicNum() == 6)


def _hydrogen_slots(mol, parent, suffix_atoms):
    resolved = _hydride_ranks(mol, parent, suffix_atoms)
    if resolved is None:
        return [], {}
    ranks, fragment, order = resolved
    slots = [
        p
        for i, p in enumerate(order)
        if p in parent
        and fragment.GetAtomWithIdx(i).GetTotalNumHs()
        and _is_substitutable(mol, mol.GetAtomWithIdx(p), suffix_atoms)
    ]
    return slots, ranks


_ARRANGEMENT_NODES = 4000
_ARRANGEMENT_LEAVES = 300


def _outside(mol, root, parent):
    seen, stack = {root}, [root]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() not in parent and n.GetIdx() not in seen:
                seen.add(n.GetIdx())
                stack.append(n.GetIdx())
    return seen


def has_single_arrangement(mol, parent, mobile_charges=frozenset()):
    """True when the groups on the parent hydride `parent` (atom indices) can be placed on it in one way only, up to
    the symmetry of the parent: no isomer arises by moving a prefix or suffix to another position or by interchanging
    two of them, so no locant is needed (P-14.3.4.4). The `mobile_charges` atoms are cationic centres whose added
    hydron may sit on any atom of the parent, so the centre itself is one of the groups that can move."""
    if mol is None or CITE_SUBSTITUENT_LOCANTS.get():
        return False
    try:
        return _single_arrangement(mol, set(parent), set(mobile_charges))
    except Exception:
        return False


def _single_arrangement(mol, parent, mobile=frozenset()):
    groups, positions = [], sorted(parent)
    for p in positions:
        for n in mol.GetAtomWithIdx(p).GetNeighbors():
            root = n.GetIdx()
            if root in parent:
                continue
            inside = _outside(mol, root, parent)
            if any(m.GetIdx() in parent and (m.GetIdx(), a) != (p, root) for a in inside for m in mol.GetAtomWithIdx(a).GetNeighbors()):
                return False
            order = mol.GetBondBetweenAtoms(p, root).GetBondTypeAsDouble()
            if order not in (1.0, 2.0, 3.0):
                return False
            signature = (order, Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(inside), rootedAtAtom=root, isomericSmiles=False))
            groups.append((signature, p, root, int(order)))
    if not groups:
        return False
    free = {}
    for p in positions:
        atom = mol.GetAtomWithIdx(p)
        total = atom.GetTotalNumHs() + sum(order for _, q, _, order in groups if q == p)
        formyl = atom.GetTotalNumHs() and any(
            order == 2 and mol.GetAtomWithIdx(r).GetAtomicNum() == 8 for _, q, r, order in groups if q == p
        )
        total -= p in mobile
        free[p] = (total, total - (atom.GetTotalNumHs() if formyl else 0))
    groups.sort(key=lambda g: g[0])
    base = Chem.RWMol(mol)
    for _, p, root, _ in groups:
        base.RemoveBond(p, root)
    seen, budget = set(), [_ARRANGEMENT_NODES]

    def build(placement, charged):
        editable = Chem.RWMol(base)
        used = {p: 0 for p in positions}
        for (_, _, root, order), target in zip(groups, placement):
            editable.AddBond(target, root, Chem.BondType.values[order])
            used[target] += order
        for p in positions:
            atom = editable.GetAtomWithIdx(p)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(free[p][0] + (p in charged) - used[p])
            if mobile:
                atom.SetFormalCharge(int(p in charged))
        try:
            Chem.SanitizeMol(editable)
        except Exception:
            budget[0] = -1
            return None
        return Chem.MolToSmiles(editable, isomericSmiles=False)

    def place(index, placement, remaining, charged):
        budget[0] -= 1
        if budget[0] < 0 or len(seen) > 1:
            return
        if index == len(groups):
            found = build(placement, charged)
            if found is not None:
                seen.add(found)
                if len(seen) > _ARRANGEMENT_LEAVES:
                    budget[0] = -1
            return
        signature, order = groups[index][0], groups[index][3]
        floor = placement[-1] if index and groups[index - 1][0] == signature else None
        for p in positions:
            if floor is not None and p < floor:
                continue
            if remaining[p] >= order:
                remaining[p] -= order
                place(index + 1, placement + [p], remaining, charged)
                remaining[p] += order

    for charged in combinations(positions, len(mobile)) if mobile else [()]:
        place(0, [], {p: free[p][1] + (p in charged) for p in positions}, set(charged))
    return budget[0] >= 0 and len(seen) == 1
