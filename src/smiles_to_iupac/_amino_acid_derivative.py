"""Esters (P-103.2.6) and singly ionized forms (P-103.2.4.2) of amino acids named on the retained name.

The ester or ion is reduced to its neutral acid, named by the amino acid module, and the ending is changed:
'-ine' to '-inate' / '-inium', '-ic acid' to '-ate'; the alcohol groups precede the anion name.
"""

from rdkit import Chem

from ._amino_acid_retained import has_retained_amino_acid_shape, name_retained_amino_acid
from ._amino_acid import SYSTEMATIC_ACID_PROBE, _match, has_amino_acid_shape as _has_plain, name_amino_acid as _name_plain
from ._histidine import has_histidine_shape, name_histidine
from ._cited_group import cited_group, subtree
from ._common import UnsupportedStructure, adjacency

_TWO_AMINO_GROUPS = {"lysine", "ornithine", "arginine", "histidine"}
_DIACID_SIDE_LOCANT = {"aspartic acid": "4", "glutamic acid": "5"}
_ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2;!R]([#6])")


def _anion_stem(plain):
    for ending, replacement in (("aspartic acid", "aspartate"), ("glutamic acid", "glutamate"), ("ine", "inate"), ("tryptophan", "tryptophanate")):
        if plain.endswith(ending):
            return plain[: -len(ending)] + replacement
    return None


def _cation_stem(plain):
    for ending, replacement in (("ine", "inium"), ("tryptophan", "tryptophanium")):
        if plain.endswith(ending):
            return plain[: -len(ending)] + replacement
    return None


def _neutralize(mol):
    """(neutral copy, carboxylate count, ammonium count); None if another charge or a fragment is present."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    editable = Chem.RWMol(mol)
    carboxylates = ammoniums = 0
    for atom in editable.GetAtoms():
        charge = atom.GetFormalCharge()
        if not charge:
            continue
        if charge == -1 and atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:
            carboxylates += 1
        elif charge == 1 and atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() >= 1:
            ammoniums += 1
            atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
            atom.SetNoImplicit(True)
            atom.SetFormalCharge(0)
            continue
        else:
            return None
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    neutral = editable.GetMol()
    Chem.SanitizeMol(neutral)
    return neutral, carboxylates, ammoniums


def _acid_and_esters(mol):
    """(acid copy, [(carbonyl carbon, alcohol group)]) with each ester oxygen turned into a hydroxy group."""
    graph = adjacency(mol)
    sites = []
    for carbon, _, oxygen, root in mol.GetSubstructMatches(_ESTER):
        if mol.GetAtomWithIdx(carbon).GetAtomicNum() == 6 and root != carbon:
            sites.append((carbon, oxygen, root))
    removed = set()
    groups = []
    for carbon, oxygen, root in sites:
        inside = subtree(graph, root, oxygen)
        if carbon in inside:
            raise UnsupportedStructure("a cyclic ester is not an amino acid ester")
        removed |= inside
        groups.append((carbon, cited_group(mol, graph, root, oxygen)))
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    for _, oxygen, _ in sites:
        atom = editable.GetAtomWithIdx(oxygen)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    return acid, groups


def _ester_words(groups, locants, anion):
    if locants is None:
        return f"{groups[0][1][0]} {anion}"
    if len(groups) == 2 and groups[0][1] == groups[1][1]:
        name, compound = groups[0][1]
        word = f"bis({name})" if compound else f"di{name}"
        return f"{word} {anion}"
    parts = [f"{locants[carbon]}-{name}" for carbon, (name, _) in sorted(groups, key=lambda item: locants[item[0]])]
    return f"{' '.join(parts)} {anion}"


def _derivative_name(mol):
    neutral = _neutralize(mol)
    if neutral is None:
        return None
    neutral, carboxylates, ammoniums = neutral
    if carboxylates and ammoniums:
        return None
    try:
        acid, groups = _acid_and_esters(neutral)
    except UnsupportedStructure:
        return None
    if not (groups or carboxylates or ammoniums) or (groups and (carboxylates or ammoniums)):
        return None
    if has_histidine_shape(acid):
        if not ammoniums or groups:
            return None
        return f"{_cation_stem(name_histidine(acid))}({ammoniums}+)"
    if not _has_plain(acid):
        return None
    found = _match(acid)
    base = found[0]
    plain = _name_plain(acid)
    diacid = base in _DIACID_SIDE_LOCANT
    if carboxylates:
        if diacid and carboxylates == 1:
            stem = _anion_stem(plain)
            return stem and f"{stem}(1–)"
        if diacid != (carboxylates == 2):
            return None
        return _anion_stem(plain)
    if ammoniums:
        stem = _cation_stem(plain)
        if stem and base in _TWO_AMINO_GROUPS:
            return f"{stem}({ammoniums}+)"
        return stem if ammoniums == 1 else None
    anion = _anion_stem(plain)
    if anion is None:
        return None
    if not diacid:
        return _ester_words(groups, None, anion) if len(groups) == 1 else None
    alpha = next(a.GetIntProp("_orig") for a in acid.GetAtoms() if a.GetIdx() == found[1])
    graph = adjacency(neutral)
    locants = {}
    for carbon, _ in groups:
        locants[carbon] = "1" if alpha in graph[carbon] else _DIACID_SIDE_LOCANT[base]
    if len(groups) > 2 or len(set(locants.values())) != len(groups):
        return None
    return _ester_words(groups, locants, anion)


def has_amino_acid_shape(mol) -> bool:
    if SYSTEMATIC_ACID_PROBE.get():
        return False
    if has_retained_amino_acid_shape(mol):
        return True
    try:
        if _has_plain(mol):
            return True
    except UnsupportedStructure:
        pass
    return _derivative_name(mol) is not None


def name_amino_acid(mol) -> str:
    if has_retained_amino_acid_shape(mol):
        return name_retained_amino_acid(mol)
    try:
        if _has_plain(mol):
            return _name_plain(mol)
    except UnsupportedStructure:
        pass
    name = _derivative_name(mol)
    if name is None:
        raise UnsupportedStructure("this amino acid has no retained-name derivative form")
    return name
