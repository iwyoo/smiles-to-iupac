"""A substituent group named with the CIP descriptors of its own stereo elements cited inside the name."""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, halogen_substituents


def subtree(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def cited_group(mol, graph, root, from_atom):
    """(name, is_compound) of the group at `root` hanging off `from_atom`; every stereo element in it must be cited."""
    from ._substituents import BRANCH_STEREO, name_branch

    inside = subtree(graph, root, from_atom)
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
        name, compound = name_branch(graph, root, from_atom, halogen_substituents(mol), aromatic, mol)
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in atoms) or any(("bond", b) not in context["used"] for b in bonds):
        raise UnsupportedStructure("a stereo element of the group is not cited by any supported name")
    return name, compound
