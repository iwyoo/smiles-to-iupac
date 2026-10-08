"""Unsaturated acyclic parent hydrides of Groups 14 and 15 (P-31.1.2.2, P-14.3.4.2(d)): unbranched homogeneous chains
(disilyne, triazene, pentaaz-2-ene, hexasil-2-ene) and chains of alternating atoms (tristannaphospha-1,3-diene) take
'ene'/'yne' endings by the general method of P-31.1.1; carbon groups and halogens are cited as prefixes."""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    adjacency,
    group_substituents,
    halogen_substituents,
    substituent_locant_set_and_citation,
    unsaturation_suffix,
)
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_ROOTS = {
    7: "az", 14: "sil", 15: "phosph", 32: "germ", 33: "ars", 50: "stann", 51: "stib", 82: "plumb", 83: "bismuth",
}
_VALENCE = {7: 3, 14: 4, 15: 3, 32: 4, 33: 3, 50: 4, 51: 3, 82: 4, 83: 3}
# P-21.2.3.1: seniority order of the elements of an alternating chain; nitrogen is left to amine names
_ALTERNATING = (15, 33, 51, 83, 14, 32, 50, 82)
_A_TERMS = {
    15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma", 14: "sila", 32: "germa", 50: "stanna", 82: "plumba",
}


def _chain(mol):
    """(chain atoms in path order, inner element or None) of the unbranched chain every Group 14/15 atom of `mol` forms."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    members = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _ROOTS]
    elements = {a.GetAtomicNum() for a in members}
    if not members or len(elements) > 2 or any(a.IsInRing() or a.GetFormalCharge() or a.GetIsotope() for a in members):
        return None
    inner = None
    if len(elements) == 2:
        if 7 in elements:
            return None
        inner, terminal = sorted(elements, key=_ALTERNATING.index)
        counts = [sum(a.GetAtomicNum() == z for a in members) for z in (terminal, inner)]
        if counts[0] != counts[1] + 1:
            return None
    idxs = {a.GetIdx() for a in members}
    degree = {i: sum(n.GetIdx() in idxs for n in mol.GetAtomWithIdx(i).GetNeighbors()) for i in idxs}
    ends = [i for i in idxs if degree[i] <= 1]
    if len(members) < 2 or len(ends) != 2 or any(d > 2 for d in degree.values()):
        return None
    chain, previous = [ends[0]], None
    while True:
        following = [n.GetIdx() for n in mol.GetAtomWithIdx(chain[-1]).GetNeighbors() if n.GetIdx() in idxs and n.GetIdx() != previous]
        if not following:
            break
        previous = chain[-1]
        chain.append(following[0])
    if len(chain) != len(members):
        return None
    if inner is not None and any(mol.GetAtomWithIdx(a).GetAtomicNum() != (inner if i % 2 else terminal) for i, a in enumerate(chain)):
        return None
    return chain, inner


def _shape(mol):
    found = _chain(mol)
    if found is None:
        return None
    chain, inner = found
    orders = [mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() for a, b in zip(chain, chain[1:])]
    if not any(o > 1.0 for o in orders) or (len(chain) == 2 and inner is None and mol.GetAtomWithIdx(chain[0]).GetAtomicNum() == 7):
        return None
    chain_set = set(chain)
    for i in chain:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetTotalValence() != _VALENCE[atom.GetAtomicNum()] or atom.GetNumRadicalElectrons():
            return None
        for n in atom.GetNeighbors():
            if n.GetIdx() in chain_set:
                continue
            if mol.GetBondBetweenAtoms(i, n.GetIdx()).GetBondTypeAsDouble() != 1.0 or n.GetFormalCharge() or n.GetIsotope():
                return None
    if any(
        bond.GetStereo() != Chem.BondStereo.STEREONONE or bond.GetStereo() == Chem.BondStereo.STEREOANY
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() in chain_set and bond.GetEndAtomIdx() in chain_set
    ):
        return None
    outside = [a for a in mol.GetAtoms() if a.GetIdx() not in chain_set]
    if any(a.GetAtomicNum() not in (6, *HALOGEN_PREFIXES) for a in outside):
        return None
    return chain, inner, orders


def has_unsaturated_hydride_chain_shape(mol) -> bool:
    return _shape(mol) is not None


def name_unsaturated_hydride_chain(mol) -> str:
    chain, inner, orders = _shape(mol)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    count = len(chain)
    best = None
    for path, bond_orders in ((chain, orders), (chain[::-1], orders[::-1])):
        multiple = [(i + 1, order) for i, order in enumerate(bond_orders) if order > 1.0]
        enes = [p for p, order in multiple if order == 2.0]
        ynes = [p for p, order in multiple if order == 3.0]
        substituents = {}
        for position, atom in enumerate(path, start=1):
            for n in graph[atom]:
                if n not in set(path):
                    substituents.setdefault(position, []).append(
                        name_branch(graph, n, atom, halogens, aromatic, mol=mol, unsaturated=True)
                    )
        grouped = group_substituents(substituents)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (sorted(p for p, _ in multiple), enes, locant_set, citation)
        if best is None or key < best[0]:
            best = (key, enes, ynes, grouped)
    _, enes, ynes, grouped = best
    body, needs_a = unsaturation_suffix(enes, ynes)
    single = len(enes) + len(ynes) == 1
    total_substituents = sum(len(info["locants"]) for info in grouped.values())
    unlocanted = count == 2 or (count == 3 and not grouped and single)
    if unlocanted:
        body = body.split("-", 1)[1] if "-" in body else body
    if inner is None:
        stem = multiplying_prefix(count) + _ROOTS[mol.GetAtomWithIdx(chain[0]).GetAtomicNum()]
    else:
        terminal = mol.GetAtomWithIdx(chain[0]).GetAtomicNum()
        joined = _A_TERMS[terminal][:-1] if _A_TERMS[inner][0] in "aeiou" else _A_TERMS[terminal]
        stem = multiplying_prefix((count + 1) // 2) + joined + _A_TERMS[inner][:-1]
    stem += "a" if needs_a else ""
    parent = f"{stem}{'' if unlocanted else '-'}{body}"
    omit = (count == 2 and total_substituents == 1) or (len(grouped) == 1 and all(mol.GetAtomWithIdx(a).GetTotalNumHs() == 0 for a in chain))
    prefixes = format_substituent_prefixes(grouped, omit_locants=omit)
    return prefixes + parent
