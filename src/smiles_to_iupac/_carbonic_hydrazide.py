"""Hydrazides of carbonic and dicarbonic acids (P-66.3.5.1, P-66.3.5.2)."""

from rdkit import Chem

from ._common import UnsupportedStructure

_NAMES = {
    "NNC(=O)NN": "hydrazinecarbohydrazide",
    "NNC(=O)OC(=O)NN": "dicarbonic dihydrazide",
    "NNC(=O)NC(=O)NN": "2-imidodicarbonic dihydrazide",
}


def carbonic_hydrazide_name(mol):
    return _NAMES.get(Chem.MolToSmiles(mol))


def name_carbonic_hydrazide(mol) -> str:
    name = carbonic_hydrazide_name(mol)
    if name is None:
        raise UnsupportedStructure("not a hydrazide of carbonic or dicarbonic acid")
    return name
