"""Naming of boronic acids R-B(OH)2 (P-68.1.4.1): the organyl group R is cited before the retained name.

Any group of a class junior to the acid (P-41, P-42.3) is a prefix on R; a carboxylic, sulfur-group or
phosphorus acid outranks the boronic acid, which then becomes the borono prefix.
"""

from rdkit import Chem

from ._borinic_acid import _borinic_acid_boron_atoms, name_borinic_acid
from ._multiplicative_text import enclose
from ._hetero_prefixes import acid_outranks_boron, boron_acid_parent
from ._common import UnsupportedStructure, adjacency, halogen_substituents, nitro_group_atoms
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
    if len(boron_atoms) != 1:
        raise UnsupportedStructure("more than one boronic acid group is not supported yet")
    (boron,) = boron_atoms
    junior_boron = {a.GetIdx() for a in _borinic_acid_boron_atoms(mol)}

    in_nitro_group = set().union(*nitro_group_atoms(mol))
    for atom in mol.GetAtoms():
        if (atom.GetFormalCharge() != 0 and atom.GetIdx() not in in_nitro_group) or atom.GetIsotope() != 0 or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("charged, radical or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _BORON and atom.GetIdx() != boron.GetIdx() and atom.GetIdx() not in junior_boron:
            raise UnsupportedStructure("more than one boron atom is not supported yet")

    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8}
    if acid_outranks_boron(mol):
        raise UnsupportedStructure("a carboxylic, sulfur-group or phosphorus acid outranks the boronic acid")

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    (root,) = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with boron_acid_parent(), cited_stereo_around(mol, boron.GetIdx()):
        name, is_compound = name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol)
    prefix = enclose(name) if is_compound else name
    return f"{prefix}{centre_nuclides(boron.GetIdx(), oxygens=True)}boronic acid"


def boron_acid_is_parent(mol) -> bool:
    """A boronic or borinic acid whose boron carries only carbon groups (a boron bonded to N, S or Si is an
    amide, thioester or silyl derivative named elsewhere) and that no senior acid outranks."""
    groups = _boronic_acid_boron_atoms(mol) or _borinic_acid_boron_atoms(mol)
    return (
        bool(groups)
        and all(n.GetAtomicNum() in (6, 8) for boron in groups for n in boron.GetNeighbors())
        and not acid_outranks_boron(mol)
    )


def name_boron_acid(mol) -> str:
    return name_boronic_acid(mol) if has_boronic_acid_shape(mol) else name_borinic_acid(mol)
