"""Naming of simple phosphonium cations (the '-phosphanium' suffix, P-73.1.1.2): the groups on phosphorus are cited as
prefixes of 'phosphanium', the cation of the parent hydride 'phosphane' (PH3), by `onium_name` (as for `_sulfonium.py`).
P-93.2.2: a phosphorus with four different ligands is a stereogenic centre; its R or S leads the name without a locant.
Out of scope (UnsupportedStructure): a charge other than +1, more than one charged atom, isotopic modification, and a
phosphorus substituent other than carbon.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._onium_prefixes import centre_descriptor, centre_is_specified_stereocentre, only_nitro_charges_besides, onium_name


def has_phosphonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged phosphorus
    shaped like a genuine phosphonium (PH4+, a phosphorus singly bonded to
    1-3 carbons with the rest hydrogens, or a phosphorus singly bonded to
    4 carbons). Used by `core.py` to route here before `_phosphane.py`,
    which rejects any charged atom outright."""
    charged_phosphorus = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() == 1
    ]
    if len(charged_phosphorus) != 1:
        return False
    phosphorus = charged_phosphorus[0]
    if phosphorus.GetIsotope() != 0 or not only_nitro_charges_besides(mol, phosphorus):
        return False
    degree = phosphorus.GetDegree()
    if degree > 4 or phosphorus.GetTotalNumHs() + degree != 4:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in phosphorus.GetNeighbors()
    )


def name_phosphonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_phosphorus = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() != 0
    ]
    (phosphorus,) = charged_phosphorus
    if phosphorus.GetFormalCharge() != 1 or phosphorus.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "phosphonium phosphorus is supported (P-73.1.1.2)"
        )
    descriptor = centre_descriptor(mol, phosphorus)
    if centre_is_specified_stereocentre(mol, phosphorus) and not descriptor:
        raise UnsupportedStructure("the phosphonium phosphorus has no CIP descriptor")
    if any(n.GetAtomicNum() != 6 for n in phosphorus.GetNeighbors()):
        raise UnsupportedStructure(
            "a phosphonium substituent other than carbon is out of scope "
            "for this module"
        )
    if any(
        mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
        for n in phosphorus.GetNeighbors()
    ):
        raise UnsupportedStructure("the phosphonium phosphorus must be singly bonded to each substituent")

    return descriptor + onium_name(mol, phosphorus, "phosphanium")
