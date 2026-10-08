"""An acyclic cation with further onium centres of other elements cited as prefixes (P-73.6, P-73.7).

The parent is the cationic centre of the senior element (N > P > As > Sb > Bi > O > S > Se > Te); every junior centre with
its own substituents becomes an 'azaniumyl'-type prefix named from the isolated cation (method 1 of P-73.6). The junior
group is swapped for a placeholder ether while the parent is named, and its prefix is supplied in place of the ether name.
"""

from rdkit import Chem

from ._ammonium import name_ammonium
from ._common import UnsupportedStructure, adjacency
from ._oxonium import name_oxonium
from ._phosphonium import name_phosphonium
from ._substituents import FORCED_BRANCH_NAMES
from ._sulfonium import name_sulfonium
from ._zwitterion import _ONIUM_STEMS, _ammonium_prefix, _chain_neighbor

_SENIORITY = (7, 15, 33, 51, 83, 8, 16, 34, 52)
_PARENT_NAMERS = {7: name_ammonium, 15: name_phosphonium, 8: name_oxonium, 16: name_sulfonium}


def _is_onium(mol, atom):
    if atom.GetFormalCharge() != 1 or atom.IsInRing() or atom.GetIsotope():
        return False
    valence = 4 if atom.GetAtomicNum() == 7 else _ONIUM_STEMS.get(atom.GetAtomicNum(), (None, None))[1]
    orders = [mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in atom.GetNeighbors()]
    return (
        valence is not None
        and sum(orders) + atom.GetTotalNumHs() == valence
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
        and sorted(orders)[:-1] == [1.0] * (len(orders) - 1)
        and orders
        and max(orders) <= (2.0 if atom.GetAtomicNum() in (8, 16) else 1.0)
    )


def _ylidene(mol, atom):
    return any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0 for n in atom.GetNeighbors())


def _centres(mol):
    """(senior, juniors) of the cationic centres, or None when the molecule is not such a cation."""
    if len(Chem.GetMolFrags(mol)) > 1 or any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()):
        return None
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge() > 0]
    if len(cations) < 2 or not all(_is_onium(mol, a) for a in cations):
        return None
    rank = {z: i for i, z in enumerate(_SENIORITY)}
    best = min(rank[a.GetAtomicNum()] for a in cations)
    seniors = [a for a in cations if rank[a.GetAtomicNum()] == best]
    if len(seniors) > 1 and len([a for a in seniors if _ylidene(mol, a)]) == 1:
        seniors = [a for a in seniors if _ylidene(mol, a)]
    if len(seniors) != 1 or seniors[0].GetAtomicNum() not in _PARENT_NAMERS:
        return None
    return seniors[0], [a for a in cations if a is not seniors[0]]


def has_mixed_onium_shape(mol) -> bool:
    return _centres(mol) is not None


def name_mixed_onium(mol) -> str:
    senior, juniors = _centres(mol)
    graph = adjacency(mol)
    prefixes, groups = {}, {}
    for junior in juniors:
        anchor = _chain_neighbor(mol, junior.GetIdx(), senior.GetIdx())
        if anchor is None:
            raise UnsupportedStructure("an onium centre bonded to the parent cation")
        text = _ammonium_prefix(mol, junior.GetIdx(), anchor)
        simple = text in ("azaniumyl", "oxidaniumyl", "sulfaniumyl", "phosphaniumyl")
        prefixes[junior.GetIdx()] = (text if simple else text[1:-1], not simple)
        group, stack = {junior.GetIdx()}, [n for n in graph[junior.GetIdx()] if n != anchor]
        while stack:
            current = stack.pop()
            if current not in group:
                group.add(current)
                stack.extend(n for n in graph[current] if n != junior.GetIdx())
        groups[junior.GetIdx()] = (group, anchor)
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    placeholders = {}
    for junior, (group, anchor) in groups.items():
        oxygen = editable.AddAtom(Chem.Atom(8))
        carbon = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(anchor, oxygen, Chem.BondType.SINGLE)
        editable.AddBond(oxygen, carbon, Chem.BondType.SINGLE)
        editable.GetAtomWithIdx(oxygen).SetIntProp("_orig", -1 - junior)
        editable.GetAtomWithIdx(carbon).SetIntProp("_orig", -1000)
        placeholders[-1 - junior] = junior
    for group, _ in groups.values():
        for index in sorted(group, reverse=True):
            editable.RemoveAtom(index)
    reduced = editable.GetMol()
    Chem.SanitizeMol(reduced)
    forced = {a.GetIdx(): prefixes[placeholders[a.GetIntProp("_orig")]] for a in reduced.GetAtoms() if a.GetIntProp("_orig") in placeholders}
    token = FORCED_BRANCH_NAMES.set((reduced.GetNumAtoms(), forced))
    try:
        if _ylidene(mol, senior):
            from ._hydride_ylium import name_hydride_onium

            return name_hydride_onium(reduced)
        return _PARENT_NAMERS[senior.GetAtomicNum()](reduced)
    finally:
        FORCED_BRANCH_NAMES.reset(token)
