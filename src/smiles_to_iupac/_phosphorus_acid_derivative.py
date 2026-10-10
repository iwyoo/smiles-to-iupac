"""Acid halides and amides of the mononuclear Group 15 oxoacids (P-67.1.2.5.1, P-67.1.2.6.1): the centre (P, As, Sb)
bears k = 0, 1 or 2 organyl groups, one double-bonded oxygen or sulfur, and 3 - k identical halogens or amino groups.
Halides take the class name added to the acid name, 'phenylphosphonic dichloride', with 'phosphoryl trichloride' for
the halides of phosphoric acid itself; amides replace 'acid' by 'amide', with N locants for the amino groups and the
element symbol as locant of the organyl groups ('N,N'-dimethyl-P-phenylphosphonic diamide')."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._carbonic_family import pseudohalide_at
from ._oxoacid_acyl import is_senior_centre
from ._phosphonic_acid import _SENIOR_ACIDS, CENTER_STEMS
from ._substituents import (
    alpha_sort_key,
    cited_stereo_around,
    format_mononuclear_prefixes,
    format_substituent_prefixes,
    name_branch,
)

_HALIDE = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_MULTIPLIER = {1: "", 2: "di", 3: "tri"}
_PSEUDO = {
    "N3": "azide",
    "CN": "cyanide",
    "NC": "isocyanide",
    "NCO": "isocyanate",
    "NCS": "isothiocyanate",
    "NCSe": "isoselenocyanate",
    "NCTe": "isotellurocyanate",
}
_PRIMES = ["N", "N'", "N''"]
_SYMBOL = {15: "P", 33: "As", 51: "Sb"}


def _find(mol):
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in CENTER_STEMS]
    senior = [a for a in centers if is_senior_centre(mol, a)]
    if len(senior) != 1:
        return None
    center = senior[0]
    if center.GetDegree() != 4 or center.GetFormalCharge() or center.IsInRing():
        return None
    double = [
        n
        for n in center.GetNeighbors()
        if n.GetAtomicNum() in (8, 16)
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if len(double) != 1:
        return None
    others = [n for n in center.GetNeighbors() if n.GetIdx() != double[0].GetIdx()]
    pseudo = {n.GetIdx(): pseudohalide_at(mol, n.GetIdx(), center.GetIdx()) for n in others}
    carbons = [n for n in others if n.GetAtomicNum() == 6 and not pseudo[n.GetIdx()]]
    rest = [n for n in others if n not in carbons]
    if not rest or len(carbons) > 2:
        return None
    if any(mol.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in others):
        return None
    words = {_halide_word(mol, n, center) for n in rest}
    if None not in words and len(words) == 1:
        return center, double[0], carbons, rest, "halide"
    if all(_is_amino(n, center) for n in rest):
        return center, double[0], carbons, rest, "amide"
    amino = [n for n in rest if _is_amino(n, center)]
    halides = [n for n in rest if n.GetAtomicNum() in HALOGEN_PREFIXES]
    if len(amino) == 1 and len(halides) == len(rest) - 1 and len({n.GetAtomicNum() for n in halides}) == 1:
        return center, double[0], carbons, rest, "amidic halide"
    return None


def _halide_word(mol, atom, center):
    if atom.GetAtomicNum() in HALOGEN_PREFIXES:
        return _HALIDE[atom.GetAtomicNum()]
    pseudo = pseudohalide_at(mol, atom.GetIdx(), center.GetIdx())
    return _PSEUDO.get(pseudo) if pseudo else None


def _is_amino(nitrogen, center):
    if nitrogen.GetAtomicNum() != 7 or nitrogen.GetFormalCharge() or nitrogen.IsInRing() or nitrogen.GetIsAromatic():
        return False
    return all(n.GetAtomicNum() == 6 for n in nitrogen.GetNeighbors() if n.GetIdx() != center.GetIdx())


def has_phosphorus_acid_derivative_shape(mol) -> bool:
    return _find(mol) is not None


def _acid_word(stem, carbons, sulfur):
    thio = "thio" if sulfur else ""
    if carbons == 0:
        return f"{stem}or{thio}ic" if sulfur else f"{stem}oric"
    return f"{stem}{'on' if carbons == 1 else 'in'}o{thio}ic" if sulfur else f"{stem}{'on' if carbons == 1 else 'in'}ic"


def name_phosphorus_acid_derivative(mol) -> str:
    found = _find(mol)
    if found is None:
        raise UnsupportedStructure("no Group 15 acid halide or amide shape")
    with cited_stereo_around(mol, found[0].GetIdx()):
        return _name_derivative(mol, found)


def _name_derivative(mol, found) -> str:
    center, chalcogen, carbons, rest, kind = found
    if len(Chem.GetMolFrags(mol)) > 1 or any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("multi-fragment or isotopically modified structures are not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the phosphorus derivative")
    sulfur = chalcogen.GetAtomicNum() == 16
    if sulfur and not carbons:
        raise UnsupportedStructure("thio acid derivatives of the acid with no organyl group are not supported yet")
    stem = CENTER_STEMS[center.GetAtomicNum()]
    acid = _acid_word(stem, len(carbons), sulfur)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    count = len(rest)
    if kind == "halide":
        word = _halide_word(mol, rest[0], center)
        head = "" if carbons else f"{stem}oryl"
        prefixes = [name_branch(graph, c.GetIdx(), center.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for c in carbons]
        organyl = format_mononuclear_prefixes(prefixes) if prefixes else ""
        if not carbons:
            if sulfur:
                raise UnsupportedStructure("thio acid halides of the acid with no organyl group are not supported yet")
            return f"{head} {_MULTIPLIER[count]}{word}"
        return f"{organyl}{acid} {_MULTIPLIER[count]}{word}"
    if kind == "amidic halide":
        return _amidic_halide(mol, graph, halogens, aromatic, center, carbons, rest, stem, sulfur)
    amino = []
    for nitrogen in rest:
        subs = [
            name_branch(graph, n.GetIdx(), nitrogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
            for n in nitrogen.GetNeighbors()
            if n.GetIdx() != center.GetIdx()
        ]
        amino.append(subs)
    amino.sort(key=lambda subs: (-len(subs), [alpha_sort_key(name) for name, _ in subs]))
    positions = {}
    for locant, subs in zip(_PRIMES, amino):
        if subs:
            positions[locant] = subs
    organyl = [
        name_branch(graph, c.GetIdx(), center.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for c in carbons
    ]
    if organyl:
        positions[_SYMBOL[center.GetAtomicNum()]] = organyl
    complete = not carbons and all(len(s) == 2 for s in amino) and len({tuple(s) for s in amino}) == 1
    grouped = group_substituents(positions)
    prefix = format_substituent_prefixes(grouped, omit_locants=False)
    if complete:
        prefix = _strip_locants(prefix)
    return f"{prefix}{acid} {_MULTIPLIER[count]}amide"


def _strip_locants(prefix):
    import re

    return re.sub(r"^(?:N'*,)*N'*-", "", prefix)


def _amidic_halide(mol, graph, halogens, aromatic, center, carbons, rest, stem, sulfur):
    """One amino group beside identical halogens: 'N,N-dimethylphosphoramidic dichloride' (P-67.1.2.5.1)."""
    if sulfur:
        raise UnsupportedStructure("thio amidic halides are not supported yet")
    nitrogen = next(n for n in rest if n.GetAtomicNum() == 7)
    halides = [n for n in rest if n is not nitrogen]
    positions = {}
    subs = [
        name_branch(graph, n.GetIdx(), nitrogen.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
        for n in nitrogen.GetNeighbors()
        if n.GetIdx() != center.GetIdx()
    ]
    if subs:
        positions["N"] = subs
    organyl = [
        name_branch(graph, c.GetIdx(), center.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True) for c in carbons
    ]
    if organyl:
        positions[_SYMBOL[center.GetAtomicNum()]] = organyl
    prefix = format_substituent_prefixes(group_substituents(positions)) if positions else ""
    acid = f"{stem}{'or' if not carbons else 'on' if len(carbons) == 1 else 'in'}amidic"
    return f"{prefix}{acid} {_MULTIPLIER[len(halides)]}{_HALIDE[halides[0].GetAtomicNum()]}"
