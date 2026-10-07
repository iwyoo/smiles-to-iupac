"""Phosphane, arsane and stibane chalcogenides ('triphenyl-λ5-phosphanone', 'trimethyl-λ5-arsanone',
'triphenyl-λ5-phosphanethione'), per the IUPAC 2013 Recommendations (P-68.3.2.3.1, P-74):

A group 15 atom carrying one double-bonded terminal chalcogen is named substitutively on its hydride with the suffix
'-one', '-thione', '-selone' or '-tellone' (the final 'e' of the hydride is elided before the vowel, P-16.3.3), which
outranks functional class nomenclature ('triphenylphosphane oxide'): 'phenylphosphanone (PIN)',
'triphenyl-λ5-phosphanone (PIN)'. Normal valence 3 needs no λ label; two or more other substituents bring the valence
above 3, flagged as 'λ5' (P-14.1). The other substituents are named by the shared substituent namer, so any group
that namer supports (alkyl, aryl, ring, halogenated) is accepted; a group that must be cited as a suffix instead is
left to the other namers.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_mononuclear_prefixes, name_branch

_CENTERS = {15: "phosphane", 33: "arsane", 51: "stibane", 83: "bismuthane"}
_SUFFIX = {8: "one", 16: "thione", 34: "selone", 52: "tellone"}
_SENIOR = [
    Chem.MolFromSmarts(smarts)
    for smarts in ("[CX3](=O)[OX2]", "[CX3](=O)[NX3]", "[CX2]#[NX1]", "[CX3H1](=O)[#6]", "[#16X4](=O)(=O)[OX2]")
]


def _find_center(mol):
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _CENTERS]
    if len(centers) != 1 or centers[0].IsInRing():
        return None
    center = centers[0]
    chalcogens = [
        b.GetOtherAtom(center)
        for b in center.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(center).GetAtomicNum() in _SUFFIX
    ]
    if len(chalcogens) != 1 or chalcogens[0].GetDegree() != 1 or chalcogens[0].GetFormalCharge():
        return None
    if any(n.GetAtomicNum() != 6 for n in center.GetNeighbors() if n.GetIdx() != chalcogens[0].GetIdx()):
        return None
    return center, chalcogens[0]


def has_phosphanone_shape(mol) -> bool:
    return _find_center(mol) is not None


def name_phosphanone(mol) -> str:
    found = _find_center(mol)
    if found is None:
        raise UnsupportedStructure("no group 15 chalcogenide (R3E=X) shape found")
    center, chalcogen = found
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if center.GetFormalCharge() or any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR):
        raise UnsupportedStructure("an acid, ester, amide, nitrile or aldehyde group outranks the chalcogenide")
    graph = adjacency(mol)
    roots = [n for n in graph[center.GetIdx()] if n != chalcogen.GetIdx()]
    if not roots:
        raise UnsupportedStructure("an unsubstituted group 15 chalcogenide is not supported yet")
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [name_branch(graph, r, center.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for r in roots]
    parent = _CENTERS[center.GetAtomicNum()]
    suffix = _SUFFIX[chalcogen.GetAtomicNum()]
    stem = parent[:-1] if suffix[0] in "aeiouy" else parent
    lambda_label = "-λ5-" if len(entries) >= 2 else ""
    return f"{format_mononuclear_prefixes(entries)}{lambda_label}{stem}{suffix}"
