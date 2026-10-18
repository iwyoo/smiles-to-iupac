"""Naming of borinic acids R2B-OH (P-68.1.4.1): both organyl groups are cited as prefixes on the retained
name, alphanumerically and with multiplying prefixes (P-16.5.1.3.1).

As for boronic acids, junior groups on R are prefixes and a senior acid turns the borinic acid into a prefix.
"""

from rdkit import Chem

from ._hetero_prefixes import acid_outranks_boron, boron_acid_parent
from ._parent_choice import senior_parent_name
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents, nitro_group_atoms
from ._substituents import cited_stereo_around, format_mononuclear_prefixes, name_branch

_BORON = 5


def _aminooxy_oxygen(atom):
    """An oxygen between the boron and an acyclic amine nitrogen: that nitrogen gives no traditional ester, so the group is
    the prefix 'aminooxy' (P-65.6.3.4.1)."""
    nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
    return (
        atom.GetAtomicNum() == 8
        and atom.GetDegree() == 2
        and len(nitrogens) == 1
        and not nitrogens[0].IsInRing()
        and not nitrogens[0].GetFormalCharge()
        and all(b.GetBondTypeAsDouble() == 1.0 for b in nitrogens[0].GetBonds())
    )


def _borinic_acid_boron_atoms(mol):
    """Boron atoms shaped like a borinic acid group: one single-bonded hydroxyl
    oxygen (terminal, one H) and two substituents other than oxygen, of which
    one may be the hydrogen atom of H-B(OH)-R."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetFormalCharge() != 0 or atom.IsInRing():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8 and not _aminooxy_oxygen(n)]
        non_oxygens = [n for n in neighbors if n not in oxygens]
        if len(oxygens) != 1 or len(non_oxygens) + atom.GetTotalNumHs() != 2:
            continue
        if atom.GetTotalNumHs() and (atom.GetTotalNumHs() > 1 or non_oxygens[0].GetAtomicNum() not in (6, 8)):
            continue
        if any(n.GetAtomicNum() in HALOGEN_PREFIXES for n in non_oxygens):
            continue
        (hydroxyl_o,) = oxygens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), hydroxyl_o.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def boronic_type_boron_atoms(mol):
    """Acyclic boron atoms R-B(OX)2 (X = H or carbon): a boronic acid, a partial ester or an ester of it."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _BORON or atom.GetDegree() != 3 or atom.GetFormalCharge() != 0 or atom.IsInRing():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8 and not _aminooxy_oxygen(n)]
        others = [n for n in neighbors if n not in oxygens]
        if len(oxygens) == 2 and others[0].GetAtomicNum() == 6 and all(
            mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0 for o in oxygens
        ):
            matches.append(atom)
    return matches


def boronic_type_group_outranks(mol, parent_boron, boronic_family):
    """Another boronic-type group ranks with or above the acid at `parent_boron`: every one for a boronic acid (an acid and
    its esters are named as one polybasic acid, P-65.6.3.2.3), one still carrying a hydroxy group for a borinic acid (P-42.3)."""
    return any(
        b.GetIdx() != parent_boron
        and (boronic_family or any(o.GetTotalNumHs() for o in b.GetNeighbors() if o.GetAtomicNum() == 8))
        for b in boronic_type_boron_atoms(mol)
    )


def boron_acid_scope_problem(mol, charged_ligands=frozenset()):
    """Why a boron acid group of `mol` is not named as the parent here: an ion, a radical, an isotope, a senior acid or a
    second fragment; None when none applies. Atoms of `charged_ligands` (the azido group) may carry a charge."""
    in_nitro_group = set().union(*nitro_group_atoms(mol)) | set(charged_ligands)
    for atom in mol.GetAtoms():
        if (atom.GetFormalCharge() != 0 and atom.GetIdx() not in in_nitro_group) or atom.GetIsotope() != 0 or atom.GetNumRadicalElectrons():
            return "charged, radical or isotopically modified atoms are not supported yet"
    if acid_outranks_boron(mol):
        return "a carboxylic, sulfur-group or phosphorus acid outranks the boron acid"
    if len(Chem.GetMolFrags(mol)) > 1:
        return "multi-fragment structures are not supported yet"
    return None


def has_borinic_acid_shape(mol) -> bool:
    return bool(_borinic_acid_boron_atoms(mol))


def name_borinic_acid(mol) -> str:
    boron_atoms = _borinic_acid_boron_atoms(mol)
    if not boron_atoms:
        raise UnsupportedStructure("no borinic acid group")
    if boronic_type_group_outranks(mol, -1, False):
        raise UnsupportedStructure("a boronic acid group outranks the borinic acid (P-42.3)")
    borons = {a.GetIdx() for a in boron_atoms}

    problem = boron_acid_scope_problem(mol)
    if problem:
        raise UnsupportedStructure(problem)

    if len(borons) == 1:
        return _name_at(mol, boron_atoms[0])[1]
    return senior_parent_name(mol, sorted(borons), lambda idx: _name_at(mol, mol.GetAtomWithIdx(idx)))


def _name_at(mol, boron):
    """(number of substituent prefixes, name) with the borinic acid group at `boron` as the parent."""
    group_oxygens = {n.GetIdx() for n in boron.GetNeighbors() if n.GetAtomicNum() == 8 and not _aminooxy_oxygen(n)}
    graph = adjacency(mol)
    roots = [n for n in graph[boron.GetIdx()] if n not in group_oxygens]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    with boron_acid_parent(), cited_stereo_around(mol, boron.GetIdx()):
        entries = [name_branch(graph, root, boron.GetIdx(), halogens, aromatic_atoms, mol=mol) for root in roots]
    return len(entries), f"{format_mononuclear_prefixes(entries)}borinic acid"
