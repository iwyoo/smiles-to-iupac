"""Naming 1H-indene and 2,3-dihydro-1H-indene (P-25.1.1 retained fusion parent
with indicated hydrogen, P-31.1.4.2.4), by exact whole-molecule match like
`_fluorene_parent.py`: the five-membered ring's sp3 carbon keeps the bicycle
out of every aromatic-ring module.
"""

from rdkit import Chem

_NAMES = {
    Chem.CanonSmiles("C1=Cc2ccccc2C1"): "1H-indene",
    Chem.CanonSmiles("C1Cc2ccccc2C1"): "2,3-dihydro-1H-indene",
    Chem.CanonSmiles("C1=c2ccccc2=CC1"): "2H-indene",
}


def has_indene_parent_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _NAMES


def name_indene_parent(mol) -> str:
    return _NAMES[Chem.MolToSmiles(mol)]
