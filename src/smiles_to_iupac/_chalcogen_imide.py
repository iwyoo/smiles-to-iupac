"""Sulfimides, sulfinylamines, sulfonylamines, sulfoximides, sulfonediimines and sulfur di- and triimides, with their
selenium and tellurium analogues (P-68.4.3.3 to P-68.4.3.8): a chalcogen centre of bonding number 4 or 6 with a doubly
bonded nitrogen is named on the lambda-sulfane (selane, tellane) parent. A doubly bonded chalcogen is the suffix 'one'
and the =NR groups are 'imino' prefixes; otherwise the =NH groups are the suffix 'imine', 'diimine' or 'triimine'.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, group_substituents, halogen_substituents, multiplied_word
from ._multiplicative_text import enclose
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_PARENT = {16: "sulfane", 34: "selane", 52: "tellane"}
_ONE = {8: "one", 16: "thione", 34: "selone", 52: "tellone"}
_N_LOCANTS = ("N", "N'", "N''")


def _terminal_chalcogen(bond, center):
    atom = bond.GetOtherAtom(center)
    return atom.GetAtomicNum() in _ONE and atom.GetDegree() == 1 and not atom.GetFormalCharge() and bond.GetBondTypeAsDouble() == 2.0


def _imino_nitrogen(bond, center):
    atom = bond.GetOtherAtom(center)
    return (
        atom.GetAtomicNum() == 7
        and bond.GetBondTypeAsDouble() == 2.0
        and not atom.GetFormalCharge()
        and not atom.IsInRing()
        and atom.GetDegree() <= 2
        and all(b.GetBondTypeAsDouble() == 1.0 and b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds() if b.GetOtherAtom(atom).GetIdx() != center.GetIdx())
    )


def _find(mol):
    """(centre, imino nitrogens, doubly bonded chalcogens, carbon neighbours) when `mol` is only a chalcogen centre
    with nonstandard bonding, its =N, =O (=S, =Se, =Te) and carbon groups, else None."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    centers = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _PARENT and not a.IsInRing() and any(_imino_nitrogen(b, a) for b in a.GetBonds())
    ]
    if len(centers) != 1:
        return None
    center = centers[0]
    if center.GetFormalCharge() or center.GetIsotope() or center.GetNumRadicalElectrons() or center.GetIsAromatic():
        return None
    imines, chalcogens, carbons = [], [], []
    for bond in center.GetBonds():
        other = bond.GetOtherAtom(center)
        if _imino_nitrogen(bond, center):
            imines.append(other)
        elif _terminal_chalcogen(bond, center):
            chalcogens.append(other)
        elif other.GetAtomicNum() == 6 and bond.GetBondTypeAsDouble() == 1.0:
            carbons.append(other)
        else:
            return None
    valence = center.GetTotalValence()
    if not imines or valence not in (4, 6) or len({a.GetAtomicNum() for a in chalcogens}) > 1:
        return None
    allowed = {center.GetIdx(), *(a.GetIdx() for a in imines), *(a.GetIdx() for a in chalcogens)}
    if any(
        a.GetIdx() not in allowed and (a.GetAtomicNum() != 6 and a.GetAtomicNum() not in HALOGEN_PREFIXES or a.GetFormalCharge() or a.GetIsotope())
        for a in mol.GetAtoms()
    ):
        return None
    return center, imines, chalcogens, carbons


def has_chalcogen_imide_shape(mol) -> bool:
    return _find(mol) is not None


def _distinct_arrangements(sulfur_slots, nitrogen_slots, names):
    """The number of different compounds obtained by placing the substituents on the sulfur hydrogens or on the imino
    nitrogens (one each), the nitrogens being interchangeable (P-14.3.4.4)."""
    found = set()

    def place(index, on_sulfur, on_nitrogen):
        if index == len(names):
            found.add((tuple(sorted(on_sulfur)), tuple(sorted(on_nitrogen))))
            return
        if len(on_sulfur) < sulfur_slots:
            place(index + 1, on_sulfur + [names[index]], on_nitrogen)
        if len(on_nitrogen) < nitrogen_slots:
            place(index + 1, on_sulfur, on_nitrogen + [names[index]])

    place(0, [], [])
    return len(found)


def name_chalcogen_imide(mol) -> str:
    found = _find(mol)
    if found is None:
        raise UnsupportedStructure("no chalcogen imide shape found")
    center, imines, chalcogens, carbons = found
    if any(e.specified == Chem.StereoSpecified.Specified for e in Chem.FindPotentialStereo(mol)):
        raise UnsupportedStructure("the configuration of a chalcogen imide is not cited yet")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())

    def named(root, parent):
        return name_branch(graph, root, parent, halogens, aromatic, mol=mol, unsaturated=True)

    sulfur_groups = [named(c.GetIdx(), center.GetIdx()) for c in carbons]
    nitrogen_groups = []
    for nitrogen in imines:
        roots = [n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != center.GetIdx()]
        nitrogen_groups.append(named(roots[0], nitrogen.GetIdx()) if roots else None)
    valence = center.GetTotalValence()
    parent = _PARENT[center.GetAtomicNum()]

    if chalcogens:
        entries = list(sulfur_groups)
        for group in nitrogen_groups:
            entries.append(("imino", False) if group is None else ((enclose(group[0]) if group[1] else group[0]) + "imino", True))
        suffix = multiplied_word(len(chalcogens), _ONE[chalcogens[0].GetAtomicNum()])
        stem = parent[:-1] if suffix[0] in "aeiouy" else parent
        return f"{format_mononuclear_prefixes(entries)}-λ{valence}-{stem}{suffix}"

    suffix = multiplied_word(len(imines), "imine")
    stem = parent[:-1] if suffix[0] in "aeiouy" else parent
    substituents = [g for g in nitrogen_groups if g is not None]
    keys = [g for g in (*sulfur_groups, *substituents)]
    sulfur_slots = valence - 2 * len(imines)
    cite_locants = bool(substituents) and _distinct_arrangements(sulfur_slots, len(imines), keys) > 1
    if not cite_locants:
        prefixes = format_mononuclear_prefixes(keys) if keys else ""
        return f"{prefixes}{'-' if prefixes else ''}λ{valence}-{stem}{suffix}"
    positions = {}
    for group in sulfur_groups:
        positions.setdefault(center.GetSymbol(), []).append(group)
    ordered = sorted(
        (g for g in nitrogen_groups if g is not None),
        key=lambda g: g[0],
    )
    for locant, group in zip(_N_LOCANTS, ordered):
        positions.setdefault(locant, []).append(group)
    grouped = group_substituents(positions)
    return f"{format_substituent_prefixes(grouped)}-λ{valence}-{stem}{suffix}"
