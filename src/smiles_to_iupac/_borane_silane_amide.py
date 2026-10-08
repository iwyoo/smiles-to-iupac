"""Amides and hydrazides of the boron acids and silicic acid (P-67.1.2.6.2): substitutive names on the preselected parent
hydrides borane and silane, 'boranetriamine', 'silanetetramine', '1,1',1''-boranetriyltrihydrazine'."""

from rdkit import Chem

from ._cited_group import cited_group
from ._common import UnsupportedStructure, adjacency
from ._common import multiplied_word
from ._numerals import numerical_term
from ._substituents import format_mononuclear_prefixes

_PARENTS = {5: ("borane", 3, "boran"), 14: ("silane", 4, "silan")}
_PRIMES = ("", "′", "″", "‴")


def _nitrogen_group(mol, center, neighbor):
    """'amine' for -NH2, 'hydrazine' for -NH-NH2, else None."""
    if neighbor.GetAtomicNum() != 7 or neighbor.GetFormalCharge():
        return None
    others = [n for n in neighbor.GetNeighbors() if n.GetIdx() != center.GetIdx()]
    if not others and neighbor.GetTotalNumHs() == 2:
        return "amine", {neighbor.GetIdx()}
    if len(others) == 1 and neighbor.GetTotalNumHs() == 1:
        far = others[0]
        if far.GetAtomicNum() == 7 and far.GetDegree() == 1 and far.GetTotalNumHs() == 2 and not far.GetFormalCharge():
            return "hydrazine", {neighbor.GetIdx(), far.GetIdx()}
    return None


def borane_silane_amide_name(mol):
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _PARENTS]
    if len(centers) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    center = centers[0]
    parent, valence, stem = _PARENTS[center.GetAtomicNum()]
    if center.GetFormalCharge() or center.GetIsotope() or center.IsInRing():
        return None
    kinds, carbons, covered = [], [], {center.GetIdx()}
    for neighbor in center.GetNeighbors():
        if mol.GetBondBetweenAtoms(center.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
        group = _nitrogen_group(mol, center, neighbor)
        if group is not None:
            kinds.append(group[0])
            covered |= group[1]
        elif neighbor.GetAtomicNum() == 6:
            carbons.append(neighbor.GetIdx())
        else:
            return None
    if not kinds or carbons or len(set(kinds)) != 1 or len(kinds) + center.GetTotalNumHs() != valence:
        return None
    graph = adjacency(mol)
    try:
        names = [cited_group(mol, graph, c, center.GetIdx()) for c in carbons]
    except UnsupportedStructure:
        return None
    from ._cited_group import subtree

    for c in carbons:
        covered |= subtree(graph, c, center.GetIdx())
    if covered != set(range(mol.GetNumAtoms())):
        return None
    count = len(kinds)
    prefixes = format_mononuclear_prefixes(names) if names else ""
    if kinds[0] == "amine":
        multiplied = f"{parent}{multiplied_word(count, 'amine')}" if count > 1 else f"{stem}amine"
        return prefixes + multiplied
    if carbons or center.GetAtomicNum() not in _PARENTS:
        return None
    if count == 1:
        return f"{stem}ylhydrazine"
    locants = ",".join(f"1{chr(0x2032) * i}" for i in range(count))
    linker = f"{stem}e{numerical_term(count)}yl"
    return f"{locants}-{linker}{numerical_term(count)}hydrazine"
