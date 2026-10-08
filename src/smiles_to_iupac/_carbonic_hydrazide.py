"""Hydrazides, amidrazones and hydrazidines of carbonic and dicarbonic acids (P-66.3.5, P-66.4.2.2, P-66.4.3.3)."""

from rdkit import Chem

from ._common import UnsupportedStructure

_NAMES = {
    "NNC(=O)NN": "hydrazinecarbohydrazide",
    "NNC(=O)OC(=O)NN": "dicarbonic dihydrazide",
    "NNC(=O)NC(=O)NN": "2-imidodicarbonic dihydrazide",
    "NNC(=NN)NN": "hydrazinecarbohydrazonohydrazide",
    "NNC(=NN)OC(=NN)NN": "dicarbonohydrazonic dihydrazide",
    "NNC(=N)NN": "hydrazinecarboximidohydrazide",
    "NC(=N)NN": "hydrazinecarboximidamide",
    "NC(=NN)NN": "hydrazinecarbohydrazonamide",
    "NC(=NN)N": "carbonohydrazonic diamide",
    "NNC(=N)OC(=N)NN": "dicarbonimidic dihydrazide",
    "NC(=NN)OC(=NN)N": "dicarbonohydrazonic diamide",
    "NC(=NN)OC(=N)NN": "[(hydrazinecarboximidoyl)oxy]methanehydrazonamide",
}


_CANONICAL = {Chem.MolToSmiles(Chem.MolFromSmiles(smiles)): name for smiles, name in _NAMES.items()}


def carbonic_hydrazide_name(mol):
    return _CANONICAL.get(Chem.MolToSmiles(mol))


def name_carbonic_hydrazide(mol) -> str:
    name = carbonic_hydrazide_name(mol)
    if name is None:
        raise UnsupportedStructure("not a hydrazide of carbonic or dicarbonic acid")
    return name
