"""Substituents on a mononuclear onium centre ('sulfanium', 'oxidanium', 'phosphanium'; P-73.1.1.2): a cation outranks
every neutral class (P-41), so each group of the substituents is an ordinary prefix named through `name_branch`."""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    is_nitro_nitrogen,
    plain_phenyl_substituent_atoms,
    specified_stereo_elements,
)
from ._substituents import cited_branch_stereo, format_mononuclear_prefixes, name_branch


def _nitro_part(mol, atom):
    if is_nitro_nitrogen(mol, atom.GetIdx()):
        return True
    return atom.GetAtomicNum() == 8 and any(is_nitro_nitrogen(mol, n.GetIdx()) for n in atom.GetNeighbors())


def only_nitro_charges_besides(mol, center):
    """No charged atom other than `center` and the charge-separated nitro groups: an anionic group elsewhere makes the
    molecule a zwitterion named on the anion."""
    return all(
        atom.GetIdx() == center.GetIdx() or not atom.GetFormalCharge() or _nitro_part(mol, atom) for atom in mol.GetAtoms()
    )


def centre_is_specified_stereocentre(mol, center):
    """Whether the onium centre itself is a specified stereocentre; the descriptors of the groups on it are cited with
    their prefixes (P-93.4.2.1, P-91.3)."""
    return any(
        element.centeredOn == center.GetIdx() and element.specified == Chem.StereoSpecified.Specified
        for element in Chem.FindPotentialStereo(mol)
    )


def centre_descriptor(mol, center):
    """P-93.2.2: the R or S of an azanium, phosphanium, arsanium, stibanium or bismuthanium centre with four ligands,
    cited in front of the cation name without a locant."""
    if center.GetAtomicNum() not in (7, 15, 33, 51, 83) or center.GetDegree() + center.GetTotalNumHs() != 4:
        return ""
    if not centre_is_specified_stereocentre(mol, center):
        return ""
    codes = [code for kind, where, code in specified_stereo_elements(mol) or [] if kind == "atom" and where == center.GetIdx()]
    if not codes:
        raise UnsupportedStructure("could not determine a CIP descriptor for the onium centre")
    return f"({codes[0]})-"


def onium_name(mol, center, stem):
    """`stem` ('sulfanium', ...) preceded by the prefixes of every group bonded to `center`."""
    for atom in mol.GetAtoms():
        if atom.GetIdx() == center.GetIdx():
            continue
        if atom.GetIsotope() or atom.GetNumRadicalElectrons() or (atom.GetFormalCharge() and not _nitro_part(mol, atom)):
            raise UnsupportedStructure("charged, radical or isotopically modified atoms are not supported yet")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    roots = sorted(graph[center.GetIdx()])
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, set(roots))
    with cited_branch_stereo(mol, graph, {center.GetIdx()}, roots):
        names = [
            ("phenyl", False) if root in phenyl_atoms else name_branch(graph, root, center.GetIdx(), halogens, mol=mol)
            for root in roots
        ]
    return (format_mononuclear_prefixes(names) if names else "") + stem
