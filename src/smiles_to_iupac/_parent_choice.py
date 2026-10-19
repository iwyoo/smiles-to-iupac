"""Choice among equally senior parent candidates of one class (P-15.3.1.1, P-45.2.1, P-45.5)."""

import contextlib

from rdkit import Chem

from ._common import UnsupportedStructure, alpha_sort_key
from ._pin import mark


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


def senior_parent_name(mol, centres, build, prefix_context=contextlib.nullcontext):
    """The name of the candidate with the most substituent prefixes, then the alphanumerically earlier name; `build(centre)`
    gives (prefix count, name). Identically substituted groups are identical parents, named multiplicatively (P-15.3.1.1)."""
    groups = [_group_atoms(mol, c) for c in centres]
    if any(a & b for i, a in enumerate(groups) for b in groups[i + 1 :]):
        raise UnsupportedStructure("acid groups sharing an atom are not separate parent candidates")
    units = [_unit(mol, c, centres) for c in centres]
    if len(set(units)) != len(units):
        return name_one_of_several_groups(mol, centres, lambda centre: build(centre)[1], prefix_context)
    return min((build(c) for c in centres), key=lambda built: (-built[0], alpha_sort_key(built[1])))[1]


def name_one_of_several_groups(mol, centres, name_at, prefix_context):
    """The multiplicative name of identical acid groups (P-15.3.3.2.1) when a linking group allows it (P-15.3.4.1), else the
    alphanumerically first of the substitutive names with one group as the parent (P-45.5); `name_at(centre)` gives one."""
    from ._chain_multiplicative import has_chain_multiplicative_shape
    from ._polyfunctional import name_polyfunctional

    if has_chain_multiplicative_shape(mol):
        with contextlib.suppress(UnsupportedStructure), prefix_context():
            name = name_polyfunctional(mol)
            if name is not None:
                return name
    names = [name_at(centre) for centre in centres]
    if len(set(names)) == 1:
        return mark(names[0], "identical parent structures joined by a symmetrical group are named multiplicatively (P-15.3.3.2.1)")
    return min(names, key=lambda name: (alpha_sort_key(name), name))
