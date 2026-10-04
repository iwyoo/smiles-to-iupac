"""Glycosyl groups of O-glycosides (P-102.4, P-107.4.2, P-107.4.3.3).

The sugar part is cut at its anomeric C-O bond, named as the free aldopyranose, aldofuranose or ketopyranose,
and its '-ose' ending becomes '-osyl'.
"""

from rdkit import Chem

from ._carbohydrate import (
    has_cyclic_aldofuranose_shape,
    has_cyclic_aldopyranose_shape,
    has_cyclic_ketohexopyranose_shape,
    name_cyclic_aldofuranose,
    name_cyclic_aldopyranose,
    name_cyclic_ketohexopyranose,
)

_NAMERS = (
    (has_cyclic_aldopyranose_shape, name_cyclic_aldopyranose),
    (has_cyclic_aldofuranose_shape, name_cyclic_aldofuranose),
    (has_cyclic_ketohexopyranose_shape, name_cyclic_ketohexopyranose),
)


def _sugar_atoms(graph, root, glycosidic_oxygen):
    seen = {root}
    stack = [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != glycosidic_oxygen and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def glycosyl_group(mol, graph, root, glycosidic_oxygen):
    """(name, atoms) of the glycosyl group rooted at the anomeric `root` and bonded to the aglycone through
    `glycosidic_oxygen`, or None when the group is not a plain unsubstituted aldose/ketose ring."""
    if mol.GetAtomWithIdx(glycosidic_oxygen).GetAtomicNum() != 8 or not mol.GetAtomWithIdx(root).IsInRing():
        return None
    atoms = _sugar_atoms(graph, root, glycosidic_oxygen)
    if glycosidic_oxygen in atoms or len(atoms) > 20:
        return None
    editable = Chem.RWMol(mol)
    keep = atoms | {glycosidic_oxygen}
    for idx in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    sugar = editable.GetMol()
    try:
        Chem.SanitizeMol(sugar)
    except Exception:
        return None
    for has_shape, namer in _NAMERS:
        if has_shape(sugar):
            name = namer(sugar)
            return name[: -len("e")] + "yl", atoms
    return None


def glycosyl_branch(mol, graph, root, coming_from):
    """(name, True) of a glycosyl substituent, marking its stereocentres as cited in the enclosing name."""
    from ._substituents import BRANCH_STEREO

    found = glycosyl_group(mol, graph, root, coming_from)
    if found is None:
        return None
    name, atoms = found
    context = BRANCH_STEREO.get()
    if context:
        context["used"].update(("atom", a) for a in atoms)
    return name, True
