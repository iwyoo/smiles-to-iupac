"""Directly bonded Group 13-15 metal pairs (P-69.5.3, `tmp/bluebook/P6a.txt`
lines 8971-8985): the senior metal (As > Sb > Bi > Ge > Sn > Pb > Al > Ga >
In > Tl) is the parent hydride, other metals are '-yl' groups
('germylbismuthane'). Carbon chains linking metals and the multiplicative
'plumbanetetrayl' form are out of scope.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds, plain_phenyl_substituent_atoms
from ._group13_hydride import GROUP_13_STEMS, GROUP_14_STEMS, GROUP_15_STEMS
from ._substituents import format_mononuclear_prefixes, name_branch

_STEMS = {**GROUP_13_STEMS, **GROUP_14_STEMS, **GROUP_15_STEMS}
_SENIORITY = [33, 51, 83, 32, 50, 82, 13, 31, 49, 81]
_MAX_VALENCE = {**{n: 3 for n in GROUP_13_STEMS}, **{n: 4 for n in GROUP_14_STEMS}, **{n: 3 for n in GROUP_15_STEMS}}


def has_metal_pair_shape(mol) -> bool:
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS) >= 2


def _group_stem(atomic_num: int) -> str:
    stem = _STEMS[atomic_num]
    if atomic_num in GROUP_14_STEMS:
        return stem[:-3] + "yl"
    return stem[:-1] + "yl"


def _substituents(mol, graph, phenyl_atoms, idx, came_from):
    entries = []
    for n in graph[idx]:
        if n == came_from:
            continue
        atom = mol.GetAtomWithIdx(n)
        if atom.GetAtomicNum() in _STEMS:
            group = _metal_group_name(mol, graph, phenyl_atoms, n, idx)
            entries.append((group, "\x01" in group))
        elif atom.GetAtomicNum() == 6:
            if n in phenyl_atoms:
                entries.append(("phenyl", False))
            else:
                entries.append(name_branch(graph, n, idx, {}, mol=mol))
        else:
            raise UnsupportedStructure("only metal, carbon and hydrogen substituents are supported here")
    return entries


def _format(entries) -> str:
    return format_mononuclear_prefixes(entries) if entries else ""


def _brackets(name: str) -> str:
    return name.replace("(\x01", "[").replace("\x02)", "]").replace("\x01", "[").replace("\x02", "]")


def _metal_group_name(mol, graph, phenyl_atoms, idx, came_from) -> str:
    atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
    entries = _substituents(mol, graph, phenyl_atoms, idx, came_from)
    name = _format(entries) + _group_stem(atomic_num)
    return f"\x01{name}\x02" if "(" in name else name


def name_metal_pair(mol) -> str:
    metals = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS]
    parent_num = next(n for n in _SENIORITY if any(a.GetAtomicNum() == n for a in metals))
    parents = [a for a in metals if a.GetAtomicNum() == parent_num]
    if len(parents) != 1:
        raise UnsupportedStructure("more than one atom of the senior metal is not supported yet")
    (parent,) = parents
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetFormalCharge() != 0 or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if any(a.GetAtomicNum() not in (6, *_STEMS) for a in mol.GetAtoms()):
        raise UnsupportedStructure("only carbon and metal atoms are supported here")

    graph = adjacency(mol)
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, {n for m in metals for n in graph[m.GetIdx()]})
    if any(b[0] not in phenyl_atoms or b[1] not in phenyl_atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated substituent is out of scope here")
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if ring_atoms - phenyl_atoms:
        raise UnsupportedStructure("a ring other than plain phenyl is out of scope here")
    for m in metals:
        if m.GetDegree() > _MAX_VALENCE[m.GetAtomicNum()]:
            raise UnsupportedStructure("a metal atom exceeds its valence")
        for n in m.GetNeighbors():
            if n.GetAtomicNum() != 6 or n.GetIdx() in phenyl_atoms:
                continue
            seen, stack = {m.GetIdx()}, [n.GetIdx()]
            while stack:
                cur = stack.pop()
                if cur in seen:
                    continue
                seen.add(cur)
                if mol.GetAtomWithIdx(cur).GetAtomicNum() in _STEMS:
                    raise UnsupportedStructure("a carbon chain linking two metals is out of scope here")
                stack.extend(graph[cur])

    entries = _substituents(mol, graph, phenyl_atoms, parent.GetIdx(), -1)
    return _brackets(_format(entries) + _STEMS[parent_num])
