from rdkit import Chem

from smiles_to_iupac._fusion_numbering_general import general_peripheral_numbering
from smiles_to_iupac._quinoline_bicyclic_numbering import peripheral_numbering as bg_numbering


def test_matches_benzo_g_quinoline_numbering():
    mol = Chem.MolFromSmiles("C1=CC=C2C=C3C(=CC2=C1)C=CC=N3")
    assert general_peripheral_numbering(mol) == bg_numbering(mol)


def test_matches_benzo_g_isoquinoline_numbering():
    mol = Chem.MolFromSmiles("C1=CC=C2C=C3C=NC=CC3=CC2=C1")
    assert general_peripheral_numbering(mol) == bg_numbering(mol)


def test_single_ring_returns_none():
    assert general_peripheral_numbering(Chem.MolFromSmiles("c1ccccc1")) is None


def test_peri_fused_returns_none():
    # Pyrene: the ring-fusion graph has a cycle, not a tree.
    assert general_peripheral_numbering(Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")) is None
