"""Esters (P-103.2.6) and singly ionized forms (P-103.2.4.2) of amino acids named on the retained name.

The ester or ion is reduced to its neutral acid, named by the amino acid module, and the ending is changed:
'-ine' to '-inate' / '-inium', '-ic acid' to '-ate'; the alcohol groups precede the anion name.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._amino_acid import SYSTEMATIC_ACID_PROBE, _match, has_amino_acid_shape as _has_plain, name_amino_acid as _name_plain
from ._common import UnsupportedStructure, adjacency, halogen_substituents

_DIACID_SIDE_LOCANT = {"aspartic acid": "4", "glutamic acid": "5"}
_ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2;!R]([#6])")


def _subtree(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _alcohol_group(mol, graph, oxygen, root):
    from ._substituents import BRANCH_STEREO, name_branch

    inside = _subtree(graph, root, oxygen)
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in inside and a.HasProp("_CIPCode")}
    bonds = {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
        for b in probe.GetBonds()
        if b.HasProp("_CIPCode") and b.GetBeginAtomIdx() in inside and b.GetEndAtomIdx() in inside
    }
    context = {"atoms": atoms, "bonds": bonds, "used": set()}
    aromatic = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    token = BRANCH_STEREO.set(context)
    try:
        name, compound = name_branch(graph, root, oxygen, halogen_substituents(mol), aromatic, mol)
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in atoms) or any(("bond", b) not in context["used"] for b in bonds):
        raise UnsupportedStructure("a stereo element of the ester group is not cited by any supported name")
    return name, compound


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
        elif charge == 1 and atom.GetAtomicNum() == 7 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3:
            ammoniums += 1
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
        inside = _subtree(graph, root, oxygen)
        if carbon in inside:
            raise UnsupportedStructure("a cyclic ester is not an amino acid ester")
        removed |= inside
        groups.append((carbon, _alcohol_group(mol, graph, oxygen, root)))
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
    if not _has_plain(acid):
        return None
    found = _match(acid)
    base = found[0]
    plain = _name_plain(acid)
    diacid = base in _DIACID_SIDE_LOCANT
    if carboxylates:
        if diacid != (carboxylates == 2):
            return None
        return _anion_stem(plain)
    if ammoniums:
        return _cation_stem(plain)
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
    return _has_plain(mol) or _derivative_name(mol) is not None


def name_amino_acid(mol) -> str:
    if _has_plain(mol):
        return _name_plain(mol)
    return _derivative_name(mol)
