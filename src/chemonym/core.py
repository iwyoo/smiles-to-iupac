from rdkit import Chem

from ._acyclic import name_acyclic_alkane
from ._common import UnsupportedStructure
from ._cyclic import name_cycloalkane


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")

    num_rings = mol.GetRingInfo().NumRings()
    if num_rings == 0:
        return name_acyclic_alkane(mol)
    if num_rings == 1:
        return name_cycloalkane(mol)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
