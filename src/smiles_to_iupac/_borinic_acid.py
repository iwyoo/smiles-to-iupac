"""Naming of borinic acids R2B-OH (P-68.1.4.1): both organyl groups are cited as prefixes on the retained
name, alphanumerically and with multiplying prefixes (P-16.5.1.3.1).

As for boronic acids, junior groups on R are prefixes and a senior acid turns the borinic acid into a prefix.
"""

from rdkit import Chem

from ._hetero_prefixes import acid_outranks_boron, boron_acid_parent
from ._parent_choice import senior_parent_name
from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import cited_stereo_around, format_mononuclear_prefixes, name_branch

_BORON = 5


def _borinic_acid_boron_atoms(mol):
    """Boron atoms shaped like a borinic acid group: bonded to exactly two
    non-oxygen substituents and one single-bonded hydroxyl oxygen
    (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0 or atom.IsInRing():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        non_oxygens = [n for n in neighbors if n.GetAtomicNum() != 8]
        if len(oxygens) != 1 or len(non_oxygens) != 2:
            continue
        (hydroxyl_o,) = oxygens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), hydroxyl_o.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_borinic_acid_shape(mol) -> bool:
    return bool(_borinic_acid_boron_atoms(mol))


def name_borinic_acid(mol) -> str:
    from ._boronic_acid import _boronic_acid_boron_atoms

    boron_atoms = _borinic_acid_boron_atoms(mol)
    if not boron_atoms:
        raise UnsupportedStructure("no borinic acid group")
    if _boronic_acid_boron_atoms(mol):
        raise UnsupportedStructure("a boronic acid group outranks the borinic acid (P-42.3)")
    borons = {a.GetIdx() for a in boron_atoms}

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("charged, radical or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _BORON and atom.GetIdx() not in borons:
            raise UnsupportedStructure("more than one boron atom is not supported yet")

    if acid_outranks_boron(mol):
        raise UnsupportedStructure("a carboxylic, sulfur-group or phosphorus acid outranks the borinic acid")

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    if len(borons) == 1:
        return _name_at(mol, boron_atoms[0])[1]
    return senior_parent_name(mol, sorted(borons), lambda idx: _name_at(mol, mol.GetAtomWithIdx(idx)))


def _name_at(mol, boron):
    """(number of substituent prefixes, name) with the borinic acid group at `boron` as the parent."""
    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8}
    graph = adjacency(mol)
    roots = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with boron_acid_parent(), cited_stereo_around(mol, boron.GetIdx()):
        entries = [name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol) for root in roots]
    return len(entries), f"{format_mononuclear_prefixes(entries)}borinic acid"
