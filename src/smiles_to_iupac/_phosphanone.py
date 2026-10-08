"""Group 14 and 15 heterones ('dimethylsilanone', 'triphenyl-λ5-phosphanone', 'methyl-λ5-phosphanedione',
'[(2-methylpropyl)amino]arsanone'; P-61.6, P-64.1.2.2, P-64.4.1, P-74.2.1.4): a mononuclear hydride atom carrying
terminal =O/=S/=Se/=Te takes the suffix 'one'/'thione'/..., multiplied for several, and a λ label when its total
valence differs from the standard one. Carbon and amino groups are cited as prefixes; acid amides are left out.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents, heteroatom_stereo_prefix, multiplied_word
from ._substituents import format_mononuclear_prefixes, name_branch

_CENTERS = {
    14: "silane", 32: "germane", 50: "stannane", 82: "plumbane",
    15: "phosphane", 33: "arsane", 51: "stibane", 83: "bismuthane",
    16: "sulfane", 34: "selane", 52: "tellane",
}
_STANDARD_VALENCE = {14: 4, 32: 4, 50: 4, 82: 4, 15: 3, 33: 3, 51: 3, 83: 3, 16: 2, 34: 2, 52: 2}
_YLIDENE_ONLY = {16, 34, 52}
_SUFFIX = {8: "one", 16: "thione", 34: "selone", 52: "tellone"}
_SENIOR = [
    Chem.MolFromSmarts(smarts)
    for smarts in ("[CX3](=O)[OX2]", "[CX3](=O)[NX3]", "[CX2]#[NX1]", "[CX3H1](=O)[#6]", "[#16X4](=O)(=O)[OX2]")
]


def _find_center(mol):
    centers = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _CENTERS and (a.GetDegree() > 1 or a.GetAtomicNum() not in _YLIDENE_ONLY)
    ]
    if len(centers) != 1 or centers[0].IsInRing():
        return None
    center = centers[0]
    chalcogens = [
        b.GetOtherAtom(center)
        for b in center.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(center).GetAtomicNum() in _SUFFIX
    ]
    if not chalcogens or any(c.GetDegree() != 1 or c.GetFormalCharge() for c in chalcogens):
        return None
    if len({c.GetAtomicNum() for c in chalcogens}) != 1:
        return None
    others = [n for n in center.GetNeighbors() if n.GetIdx() not in {c.GetIdx() for c in chalcogens}]
    amino = [n for n in others if n.GetAtomicNum() == 7]
    if len(others) != sum(n.GetAtomicNum() == 6 for n in others) + len(amino):
        return None
    if amino and len(chalcogens) == 1 and center.GetTotalValence() == 5:
        return None
    if amino and center.GetAtomicNum() in _YLIDENE_ONLY:
        return None
    ylidene = any(
        mol.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0 for n in others
    )
    if center.GetAtomicNum() in _YLIDENE_ONLY and not ylidene:
        return None
    return center, chalcogens


def has_phosphanone_shape(mol) -> bool:
    return _find_center(mol) is not None


def name_phosphanone(mol) -> str:
    found = _find_center(mol)
    if found is None:
        raise UnsupportedStructure("no group 15 chalcogenide (R3E=X) shape found")
    center, chalcogens = found
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if center.GetFormalCharge() or any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR):
        raise UnsupportedStructure("an acid, ester, amide, nitrile or aldehyde group outranks the chalcogenide")
    graph = adjacency(mol)
    roots = [n for n in graph[center.GetIdx()] if n not in {c.GetIdx() for c in chalcogens}]
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [name_branch(graph, r, center.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for r in roots]
    parent = _CENTERS[center.GetAtomicNum()]
    suffix = multiplied_word(len(chalcogens), _SUFFIX[chalcogens[0].GetAtomicNum()])
    stem = parent[:-1] if suffix[0] in "aeiouy" else parent
    valence = center.GetTotalValence()
    lambda_label = f"λ{valence}-" if valence != _STANDARD_VALENCE[center.GetAtomicNum()] else ""
    stereo = heteroatom_stereo_prefix(mol, center.GetIdx()) or ""
    prefixes = format_mononuclear_prefixes(entries) if entries else ""
    joint = "-" if prefixes and lambda_label else ""
    return f"{stereo}{prefixes}{joint}{lambda_label}{stem}{suffix}"


_IMINE_CENTERS = {
    14: "silanimine", 32: "germanimine", 50: "stannanimine", 82: "plumbanimine",
    15: "phosphanimine", 33: "arsanimine", 51: "stibanimine", 83: "bismuthanimine",
}


def _find_imine_center(mol):
    """(center, imino nitrogen) of a Group 15 atom carrying one =N-H or =N-R and otherwise carbon groups, else None."""
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _IMINE_CENTERS and a.GetDegree() > 1]
    if len(centers) != 1 or centers[0].IsInRing() or len(Chem.GetMolFrags(mol)) != 1:
        return None
    center = centers[0]
    imines = [
        b.GetOtherAtom(center)
        for b in center.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(center).GetAtomicNum() == 7
    ]
    if len(imines) != 1 or imines[0].GetFormalCharge() or imines[0].IsInRing() or imines[0].GetDegree() > 2:
        return None
    others = [n for n in center.GetNeighbors() if n.GetIdx() != imines[0].GetIdx()]
    if not others or any(n.GetAtomicNum() != 6 or n.GetFormalCharge() for n in others):
        return None
    if any(n.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(imines[0].GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
           for n in imines[0].GetNeighbors() if n.GetIdx() != center.GetIdx()):
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    if any(mol.HasSubstructMatch(query) for query in _SENIOR):
        return None
    if any(
        b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() in (7, 8, 16)
        for n in imines[0].GetNeighbors()
        if n.GetIdx() != center.GetIdx()
        for b in n.GetBonds()
    ):
        return None
    return center, imines[0]


def has_phosphanimine_shape(mol) -> bool:
    return _find_imine_center(mol) is not None


def name_phosphanimine(mol) -> str:
    """P-62.3.1.3: =NH or =NR on a hydride atom of standard bonding number is the imine suffix (1-methylphosphanimine);
    P-74.2.1.5: R3P=NR' is a N,P-substituted lambda5-phosphanimine (N-ethyl-P,P,P-triphenyl-λ5-phosphanimine)."""
    from ._substituents import format_substituent_prefixes

    found = _find_imine_center(mol)
    if found is None:
        raise UnsupportedStructure("no Group 15 imide (R3E=NR) shape found")
    center, nitrogen = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    grouped = {}
    for locant, atom, root in (
        *(("P", center.GetIdx(), n) for n in graph[center.GetIdx()] if n != nitrogen.GetIdx()),
        *(("N", nitrogen.GetIdx(), n) for n in graph[nitrogen.GetIdx()] if n != center.GetIdx()),
    ):
        name, compound = name_branch(graph, root, atom, halogens, aromatic, mol=mol, unsaturated=True)
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    valence = center.GetTotalValence()
    parent = _IMINE_CENTERS[center.GetAtomicNum()]
    if valence == _STANDARD_VALENCE[center.GetAtomicNum()]:
        for info in grouped.values():
            info["locants"] = ["1" if locant == "P" else locant for locant in info["locants"]]
        return f"{format_substituent_prefixes(grouped)}{parent}"
    if not any("N" in info["locants"] for info in grouped.values()):
        entries = [(name, info["compound"]) for name, info in grouped.items() for _ in info["locants"]]
        return f"{format_mononuclear_prefixes(entries)}-λ{valence}-{parent}"
    return f"{format_substituent_prefixes(grouped)}-λ{valence}-{parent}"
