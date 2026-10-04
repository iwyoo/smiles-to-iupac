"""Ring groups of carbon and halogen atoms on a mononuclear hydride parent (P-29.3.3, P-29.3.4), named by the general
ring-group namer; any other element makes a senior group, which other modules name."""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._polyfunctional import _arm_atoms
from ._substituents import name_branch


def hydride_ring_groups(mol, graph, roots, center):
    """{root: (name, is_compound, atoms of the group)} for each root that starts a ring group the namer supports."""
    groups = {}
    for root in roots:
        if not mol.GetAtomWithIdx(root).IsInRing():
            continue
        atoms = _arm_atoms(graph, root, center)
        if any(
            mol.GetAtomWithIdx(a).GetFormalCharge()
            or mol.GetAtomWithIdx(a).GetIsotope()
            or (mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and mol.GetAtomWithIdx(a).GetAtomicNum() not in HALOGEN_PREFIXES)
            for a in atoms
        ):
            continue
        try:
            name, compound = name_branch(graph, root, center, {}, mol=mol)
        except UnsupportedStructure:
            continue
        groups[root] = (name, compound, atoms)
    return groups
