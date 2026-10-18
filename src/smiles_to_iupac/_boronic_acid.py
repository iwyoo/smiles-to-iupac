"""Naming of boronic acids R-B(OH)2 (P-68.1.4.1): the organyl group R is cited before the retained name.

Any group of a class junior to the acid (P-41, P-42.3) is a prefix on R; a carboxylic, sulfur-group or
phosphorus acid outranks the boronic acid, which then becomes the borono prefix. Of several boronic acid groups one is the
parent and the rest are borono prefixes (P-44.1.1).
"""


from ._borinic_acid import _borinic_acid_boron_atoms, boron_acid_scope_problem, boronic_type_boron_atoms, name_borinic_acid
from ._multiplicative_text import enclose
from ._parent_choice import name_one_of_several_groups
from ._hetero_prefixes import acid_outranks_boron, boron_acid_parent
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import centre_nuclides, cited_stereo_around, name_branch

_BORON = 5


def _boronic_acid_boron_atoms(mol):
    """Boron atoms shaped like a boronic acid group: bonded to exactly one
    non-oxygen substituent and two single-bonded hydroxyl oxygens (each
    terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 2 or len(non_oxygens) != 1:
            continue
        if non_oxygens[0].GetAtomicNum() in HALOGEN_PREFIXES:
            continue
        if any(
            mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() != 1.0 for o in oxygens
        ):
            continue
        if any(o.GetDegree() != 1 or o.GetTotalNumHs() != 1 for o in oxygens):
            continue
        matches.append(atom)
    return matches


def has_boronic_acid_shape(mol) -> bool:
    return bool(_boronic_acid_boron_atoms(mol))


def name_boronic_acid(mol) -> str:
    boron_atoms = _boronic_acid_boron_atoms(mol)
    free = {boron.GetIdx() for boron in boron_atoms}
    if any(boron.GetIdx() not in free for boron in boronic_type_boron_atoms(mol)):
        raise UnsupportedStructure("several boronic acid groups or esters are named as one polybasic acid (P-65.6.3.2.3)")

    problem = boron_acid_scope_problem(mol)
    if problem:
        raise UnsupportedStructure(problem)

    if len(boron_atoms) == 1:
        return _name_at(mol, boron_atoms[0])
    return name_one_of_several_groups(mol, boron_atoms, lambda boron: _name_at(mol, boron), boron_acid_parent)


def _name_at(mol, boron):
    """The name with the boronic acid group at `boron` as the parent; every other one is a borono prefix (P-44.1.1)."""
    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8}
    graph = adjacency(mol)
    (root,) = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with boron_acid_parent(), cited_stereo_around(mol, boron.GetIdx()):
        name, is_compound = name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol)
    prefix = enclose(name) if is_compound else name
    return f"{prefix}{centre_nuclides(boron.GetIdx(), oxygens=True)}boronic acid"


def boron_acid_is_parent(mol) -> bool:
    """A boronic or borinic acid that no senior acid outranks."""
    return bool(_boronic_acid_boron_atoms(mol) or _borinic_acid_boron_atoms(mol)) and not acid_outranks_boron(mol)


def name_boron_acid(mol) -> str:
    return name_boronic_acid(mol) if has_boronic_acid_shape(mol) else name_borinic_acid(mol)
