"""Substituents on a mononuclear onium centre ('sulfanium', 'oxidanium', 'phosphanium'; P-73.1.1.2): a cation outranks
every neutral class (P-41), so each group of the substituents is an ordinary prefix named through `name_branch`."""

from ._common import UnsupportedStructure, adjacency, is_nitro_nitrogen, plain_phenyl_substituent_atoms
from ._substituents import format_mononuclear_prefixes, name_branch


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


def onium_name(mol, center, stem):
    """`stem` ('sulfanium', ...) preceded by the prefixes of every group bonded to `center`."""
    for atom in mol.GetAtoms():
        if atom.GetIdx() == center.GetIdx():
            continue
        if atom.GetIsotope() or atom.GetNumRadicalElectrons() or (atom.GetFormalCharge() and not _nitro_part(mol, atom)):
            raise UnsupportedStructure("charged, radical or isotopically modified atoms are not supported yet")
    graph = adjacency(mol)
    roots = sorted(graph[center.GetIdx()])
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, set(roots))
    names = [
        ("phenyl", False) if root in phenyl_atoms else name_branch(graph, root, center.GetIdx(), {}, mol=mol)
        for root in roots
    ]
    return (format_mononuclear_prefixes(names) if names else "") + stem
