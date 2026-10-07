"""P-34.1.1.4: anisole is the retained PIN of methyl phenyl ether, with no substitution allowed for PINs."""

from rdkit import Chem

_ANISOLE = "COc1ccccc1"


def has_anisole_shape(mol) -> bool:
    return mol.GetNumAtoms() == 8 and Chem.MolToSmiles(mol) == _ANISOLE


def name_anisole(mol) -> str:
    return "anisole"
