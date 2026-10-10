"""Three or more contiguous identical chalcogen atoms as a parent hydride (P-68.4.1.1): methyltrisulfane,
dimethyltrisulfane, methyl(phenyl)triselane, tetramethyltetraselane. The end groups are carbon groups or hydrogen; a
chain with a characteristic group of higher seniority is a prefix of that group's parent instead."""

from rdkit import Chem

from ._common import adjacency, halogen_substituents
from ._hetero_prefixes import is_functional_carbon
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, name_branch

_STEMS = {8: "oxidane", 16: "sulfane", 34: "selane", 52: "tellane"}
_CARBON_AND_HALOGEN = {6, 9, 17, 35, 53}


def _arm_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n != blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _run(mol, graph, start):
    """The straight chain of identical chalcogen atoms containing `start`, as an ordered list, else None."""
    element = mol.GetAtomWithIdx(start).GetAtomicNum()
    members = {start}
    stack = [start]
    while stack:
        for n in graph[stack.pop()]:
            if n not in members and mol.GetAtomWithIdx(n).GetAtomicNum() == element:
                members.add(n)
                stack.append(n)
    ends = [a for a in members if sum(n in members for n in graph[a]) <= 1]
    if len(members) < 3 or len(ends) != 2 or any(sum(n in members for n in graph[a]) > 2 for a in members):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(members):
        onward = [n for n in graph[order[-1]] if n in members and n != previous]
        if len(onward) != 1:
            return None
        previous = order[-1]
        order.append(onward[0])
    return order


def name_polychalcogen_hydride(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    ):
        return None
    graph = adjacency(mol)
    order = None
    for element in _STEMS:
        members = [a for a in mol.GetAtoms() if a.GetAtomicNum() == element]
        if any(a.IsInRing() for a in members):
            return None
        found = _run(mol, graph, members[0].GetIdx()) if len(members) >= 3 else None
        if found is not None and len(found) == len(members):
            order = found
            break
    if order is None:
        return None
    run = set(order)
    for atom in order:
        a = mol.GetAtomWithIdx(atom)
        if any(b.GetBondTypeAsDouble() != 1.0 for b in a.GetBonds()):
            return None
        outside = [n for n in graph[atom] if n not in run]
        if atom not in (order[0], order[-1]) and outside:
            return None
        if len(outside) > 1 or (outside and a.GetTotalNumHs()):
            return None
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries, covered = [], set(run)
    element = mol.GetAtomWithIdx(order[0]).GetAtomicNum()
    others = set(_STEMS) - {element}
    for end in (order[0], order[-1]):
        for root in (n for n in graph[end] if n not in run):
            if mol.GetAtomWithIdx(root).GetAtomicNum() not in {6, *others}:
                return None
            arm = _arm_atoms(graph, root, end)
            if end in arm or covered & arm or any(
                mol.GetAtomWithIdx(i).GetAtomicNum() not in _CARBON_AND_HALOGEN | others for i in arm
            ):
                return None
            if any(
                b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(mol.GetAtomWithIdx(i)).GetAtomicNum() == 8
                for i in arm
                for b in mol.GetAtomWithIdx(i).GetBonds()
            ) or any(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 and is_functional_carbon(mol, i) for i in arm):
                return None
            covered |= arm
            entries.append(name_branch(graph, root, end, halogens, aromatic, mol=mol, unsaturated=True))
    if len(covered) != mol.GetNumAtoms():
        return None
    parent = multiplying_prefix(len(order)) + _STEMS[mol.GetAtomWithIdx(order[0]).GetAtomicNum()]
    if not entries:
        return parent
    if len(entries) == 1:
        (name, compound), = entries
        return (f"({name})" if compound and name[0].isdigit() else name) + parent
    return format_mononuclear_prefixes(entries) + parent
