"""Glycosyl groups of O-glycosides (P-102.4, P-107.4.2, P-107.4.3.3).

The sugar part is cut at its anomeric C-O bond, named as the free aldopyranose, aldofuranose or ketopyranose,
and its '-ose' ending becomes '-osyl'.
"""

from rdkit import Chem, rdBase

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


def _anomeric_hydroxy(mol, root, link):
    """The terminal OH on an anomeric carbon that also carries the ring oxygen, a ring carbon and the non-oxygen `link`."""
    atom = mol.GetAtomWithIdx(root)
    if atom.GetDegree() != 4 or mol.GetAtomWithIdx(link).GetAtomicNum() == 8:
        return None
    hydroxy = [
        n.GetIdx()
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1 and not n.GetFormalCharge()
    ]
    return hydroxy[0] if len(hydroxy) == 1 else None


def glycosyl_group(mol, graph, root, glycosidic_oxygen):
    """(name, atoms) of the glycosyl group rooted at the anomeric `root` and bonded to the aglycone through
    `glycosidic_oxygen`, or None when the group is not a plain unsubstituted aldose/ketose ring. An anomeric hydroxy
    group kept beside the aglycone is the prefix 1-hydroxy, and alpha/beta then refer to the free valence (P-102.6.1.1.2)."""
    if not mol.GetAtomWithIdx(root).IsInRing():
        return None
    if not any(n.GetAtomicNum() == 8 and n.IsInRing() for n in mol.GetAtomWithIdx(root).GetNeighbors()):
        return None
    atoms = _sugar_atoms(graph, root, glycosidic_oxygen)
    if glycosidic_oxygen in atoms or len(atoms) > 20:
        return None
    hydroxy = _anomeric_hydroxy(mol, root, glycosidic_oxygen)
    editable = Chem.RWMol(mol)
    # The anomeric O takes the place of the aglycone atom by changing that atom's element, which keeps the bond order
    # around the anomeric carbon and so its configuration.
    swapped = [(glycosidic_oxygen, 8)] if mol.GetAtomWithIdx(glycosidic_oxygen).GetAtomicNum() != 8 else []
    if hydroxy is not None:
        swapped.append((hydroxy, 1))
    for index, atomic_num in swapped:
        atom = editable.GetAtomWithIdx(index)
        atom.SetAtomicNum(atomic_num)
        atom.SetIsAromatic(False)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(atomic_num == 1)
    keep = atoms | {glycosidic_oxygen}
    for idx in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    sugar = editable.GetMol()
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(sugar)
            if hydroxy is not None:
                sugar = Chem.RemoveHs(sugar)
    except Exception:
        return None
    for has_shape, namer in _NAMERS:
        if has_shape(sugar):
            name = namer(sugar)
            return ("1-hydroxy-" if hydroxy is not None else "") + name[: -len("e")] + "yl", atoms
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
