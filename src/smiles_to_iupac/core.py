from rdkit import Chem

from ._acyclic import name_acyclic_alkane
from ._common import UnsupportedStructure, non_single_bonds
from ._cyclic import name_cycloalkane
from ._unsaturated import name_acyclic_unsaturated


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")

    num_rings = mol.GetRingInfo().NumRings()
    if num_rings == 0:
        bonds = non_single_bonds(mol)
        if not bonds:
            return name_acyclic_alkane(mol)
        if len(bonds) == 1 and bonds[0][2] in (2.0, 3.0):
            return name_acyclic_unsaturated(mol)
        raise UnsupportedStructure(
            "more than one multiple bond, or a bond order other than double "
            "or triple, is not supported yet (see P-31.1.1.1)"
        )
    if num_rings == 1:
        return name_cycloalkane(mol)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
