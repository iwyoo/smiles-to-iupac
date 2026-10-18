"""Phosphorus, arsenic and antimony acids R-E(=O)(OH)2, R2E(=O)OH and their -ous analogues (P-67.1.1.1, P-67.1.1.2): the
senior acid centre (P-42.3) is the parent, the organyl groups on it are cited as prefixes and every other acid centre is
a phosphono, arsono, stibono ... prefix (P-41)."""

import contextlib

from rdkit import Chem

from ._hetero_prefixes import EXTENDED_PREFIXES
from ._multiplicative_text import enclose
from ._parent_choice import name_one_of_several_groups, senior_parent_name
from ._common import UnsupportedStructure, adjacency, halogen_substituents, nitro_group_atoms
from ._substituents import centre_nuclides, cited_stereo_around, format_mononuclear_prefixes, name_branch

CENTER_STEMS = {15: "phosph", 33: "ars", 51: "stib"}


_SENIOR_ACIDS = [
    Chem.MolFromSmarts(smarts)
    for smarts in ("[CX3](=O)[OX2H1]", "[#16,#34,#52;X3,X4](=O)([#6])[OX2H1]")
]


@contextlib.contextmanager
def acid_prefixes():
    """The acid-family prefix set (sulfonyl esters, boranyl groups ...) that substituents of a principal acid may use."""
    token = EXTENDED_PREFIXES.set(True)
    try:
        yield
    finally:
        EXTENDED_PREFIXES.reset(token)


_ELEMENT_RANK = {15: 0, 33: 1, 51: 2}
_ACID_CLASS_RANK = {(1, 2): 0, (0, 2): 1, (1, 1): 2, (0, 1): 3}
_ACID_ENDINGS = {(1, 2): "onic acid", (0, 2): "onous acid", (1, 1): "inic acid", (0, 1): "inous acid"}
_NOT_AN_ACID = (1,)
_ACYL_PARTNERS = {7, 8, 9, 16, 17, 34, 35, 52, 53}


def acid_group_parts(mol, atom):
    """(oxo, hydroxyl, organyl) counts of a phosphonic, phosphonous, phosphinic or phosphinous acid centre; None for any
    other P, As or Sb atom (acid halide, amide, ester, chalcogenide, ring member, organyl derivative)."""
    if atom.GetAtomicNum() not in CENTER_STEMS or atom.GetFormalCharge() or atom.GetDegree() > 4 or atom.IsInRing():
        return None
    oxo = hydroxyl = organyl = 0
    for n in atom.GetNeighbors():
        double = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        if n.GetAtomicNum() == 8 and double and n.GetDegree() == 1:
            oxo += 1
        elif n.GetAtomicNum() == 8 and not double and n.GetDegree() == 1 and n.GetTotalNumHs() == 1:
            hydroxyl += 1
        elif n.GetAtomicNum() not in _ACYL_PARTNERS and not double:
            organyl += 1
        else:
            return None
    if hydroxyl + organyl + atom.GetTotalNumHs() != 3 or (oxo, hydroxyl) not in _ACID_CLASS_RANK:
        return None
    return oxo, hydroxyl, organyl


def acid_centre_rank(mol, atom):
    """Seniority key of a P, As or Sb atom (P-41, P-42.3): an acid centre ranks by element, then acid type; any other
    centre is a class below every acid and shares one junior key."""
    parts = acid_group_parts(mol, atom)
    if parts is None:
        return _NOT_AN_ACID
    return 0, _ELEMENT_RANK[atom.GetAtomicNum()], _ACID_CLASS_RANK[parts[:2]]


def senior_acid_centres(mol, candidates):
    """The candidates whose acid is as senior as every other acid centre of the molecule (several when they tie)."""
    ranks = {a.GetIdx(): acid_centre_rank(mol, a) for a in mol.GetAtoms() if a.GetAtomicNum() in CENTER_STEMS}
    ranked = [c for c in candidates if ranks.get(c.GetIdx(), _NOT_AN_ACID) != _NOT_AN_ACID]
    best = min((ranks[c.GetIdx()] for c in ranked), default=None)
    if best is None or best != min(ranks.values()):
        return []
    return [c for c in ranked if ranks[c.GetIdx()] == best]


def require_phosphorus_acid_scope(mol):
    """Carboxylic and the sulfur-group acids are senior to the phosphorus acids, which become phosphono prefixes (P-41)."""
    in_nitro_group = set().union(*nitro_group_atoms(mol))
    for atom in mol.GetAtoms():
        if (atom.GetFormalCharge() != 0 and atom.GetIdx() not in in_nitro_group) or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_ACIDS):
        raise UnsupportedStructure("a carboxylic or sulfur-group acid outranks the phosphorus acid")


def _phosphonic_acid_phosphorus_atoms(mol):
    """Phosphorus atoms shaped like a phosphonic acid group: bonded to
    exactly one non-oxygen substituent, one double-bonded (terminal)
    oxygen, and two single-bonded hydroxyl oxygens (each terminal, one
    H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in CENTER_STEMS or atom.GetDegree() != 4 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 3 or len(non_oxygens) != 1:
            continue
        double_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(hydroxyl_os) != 2:
            continue
        (double_o,) = double_os
        if double_o.GetDegree() != 1:
            continue
        if any(o.GetDegree() != 1 or o.GetTotalNumHs() != 1 for o in hydroxyl_os):
            continue
        matches.append(atom)
    return matches


def has_phosphorus_acid_shape(mol) -> bool:
    return any((parts := acid_group_parts(mol, atom)) is not None and parts[2] > 0 for atom in mol.GetAtoms())


def name_phosphorus_acid(mol) -> str:
    """The name with the senior acid centre (P-42.3) as the parent; identical centres on different parts are named
    through the choice among several parents (P-45.5)."""
    require_phosphorus_acid_scope(mol)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    centres = senior_acid_centres(mol, [a for a in mol.GetAtoms() if a.GetAtomicNum() in CENTER_STEMS])
    if len(centres) == 1:
        return name_acid_group(mol, centres[0])[1]
    if any(n.GetAtomicNum() in CENTER_STEMS for c in centres for n in c.GetNeighbors()):
        raise UnsupportedStructure("acid centres bonded to each other are named as one polynuclear acid")
    if acid_group_parts(mol, centres[0])[1] == 2:
        return name_one_of_several_groups(mol, centres, lambda phosphorus: name_acid_group(mol, phosphorus)[1], acid_prefixes)
    return senior_parent_name(mol, [c.GetIdx() for c in centres], lambda idx: name_acid_group(mol, mol.GetAtomWithIdx(idx)))


def name_acid_group(mol, phosphorus):
    """(number of substituent prefixes, name) with the acid group at `phosphorus` as the parent."""
    oxo, hydroxyl, _ = acid_group_parts(mol, phosphorus)
    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}
    graph = adjacency(mol)
    roots = [n for n in graph[phosphorus.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with acid_prefixes(), cited_stereo_around(mol, phosphorus.GetIdx()):
        entries = [name_branch(graph, root, phosphorus.GetIdx(), halogens, aromatic_atoms, mol=mol) for root in roots]
    centre = centre_nuclides(phosphorus.GetIdx(), oxygens=True)
    stem = f"{centre}{CENTER_STEMS[phosphorus.GetAtomicNum()]}{_ACID_ENDINGS[(oxo, hydroxyl)]}"
    if hydroxyl == 2:
        ((name, is_compound),) = entries
        return 1, f"{enclose(name) if is_compound else name}{stem}"
    return len(entries), f"{format_mononuclear_prefixes(entries)}{stem}"
