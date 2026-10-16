"""Routing of diphosphoric acid (P-67.2.1) to the polynuclear oxoacid namer.

Each phosphorus of it also has the shape of a phosphate ester whose R is the other phosphorus group, so it must be
routed here before the phosphate namers.
"""

from rdkit import Chem

from ._polynuclear_oxoacid import name_polynuclear_oxoacid

_DIPHOSPHORIC_ACID = Chem.CanonSmiles("OP(=O)(O)OP(=O)(O)O")


def has_dinuclear_oxoacid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) == _DIPHOSPHORIC_ACID


def name_dinuclear_oxoacid(mol) -> str:
    return name_polynuclear_oxoacid(mol)
