from rdkit import Chem

from smiles_to_iupac._fusion_numbering_general import general_peripheral_numbering
from smiles_to_iupac._quinoline_bicyclic_numbering import peripheral_numbering as bg_numbering


def _sp3_locant(mol, locants):
    (sp3,) = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2]
    return locants[sp3]


def test_matches_1h_cyclopenta_a_naphthalene_numbering():
    # PubChem CID 11745004, "1H-cyclopenta[a]naphthalene".
    mol = Chem.MolFromSmiles("C1C=CC2=C1C3=CC=CC=C3C=C2")
    assert _sp3_locant(mol, general_peripheral_numbering(mol)) == "1"


def test_matches_3h_cyclopenta_a_naphthalene_numbering():
    # PubChem CID 11105721, "3H-cyclopenta[a]naphthalene" -- same ring
    # topology as the 1H- isomer above, different sp3 position; the
    # structurally-privileged start atom must fix this independent of
    # which numbering would otherwise give the sp3 atom a lower locant.
    mol = Chem.MolFromSmiles("C1C=CC2=C1C=CC3=CC=CC=C32")
    assert _sp3_locant(mol, general_peripheral_numbering(mol)) == "3"


def test_matches_1h_cyclopenta_b_naphthalene_numbering():
    # PubChem CID 6451436, "1H-cyclopenta[b]naphthalene".
    mol = Chem.MolFromSmiles("C1C=CC2=CC3=CC=CC=C3C=C21")
    assert _sp3_locant(mol, general_peripheral_numbering(mol)) == "1"


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
