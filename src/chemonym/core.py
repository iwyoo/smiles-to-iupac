from rdkit import Chem

from ._acyclic import name_acyclic_alkane


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    return name_acyclic_alkane(mol)
