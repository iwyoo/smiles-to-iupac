"""Choice among equally senior parent candidates of one class (P-15.3.1.1, P-45.2.1, P-45.5)."""

from rdkit import Chem

from ._common import UnsupportedStructure, alpha_sort_key


def _group_atoms(mol, centre):
    return {centre, *(n.GetIdx() for n in mol.GetAtomWithIdx(centre).GetNeighbors() if n.GetAtomicNum() != 6)}


def _reach(mol, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() not in seen and n.GetIdx() not in blocked:
                seen.add(n.GetIdx())
                stack.append(n.GetIdx())
    return seen


def _unit(mol, centre, centres):
    """Per atom of the group, its branches leading to the other centres and its own substituents."""
    group = _group_atoms(mol, centre)
    others = set(centres) - {centre}
    unit = []
    for atom in group:
        links, own = 0, []
        for n in mol.GetAtomWithIdx(atom).GetNeighbors():
            if n.GetIdx() in group:
                continue
            branch = _reach(mol, n.GetIdx(), group)
            if branch & others:
                links += 1
            else:
                own.append(Chem.MolFragmentToSmiles(mol, sorted(branch), rootedAtAtom=n.GetIdx()))
        if links or own:
            unit.append((links, tuple(sorted(own))))
    return tuple(sorted(unit))


def senior_parent_name(mol, centres, build):
    """The name of the candidate with the most substituent prefixes, then the alphanumerically earlier name; `build(centre)`
    gives (prefix count, name). Identically substituted groups are identical parents, named multiplicatively (P-15.3.1.1)."""
    groups = [_group_atoms(mol, c) for c in centres]
    if any(a & b for i, a in enumerate(groups) for b in groups[i + 1 :]):
        raise UnsupportedStructure("acid groups sharing an atom are not separate parent candidates")
    units = [_unit(mol, c, centres) for c in centres]
    if len(set(units)) != len(units):
        raise UnsupportedStructure("identical parent structures are named multiplicatively")
    return min((build(c) for c in centres), key=lambda built: (-built[0], alpha_sort_key(built[1])))[1]
