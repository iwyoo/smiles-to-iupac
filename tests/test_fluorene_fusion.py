from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_fluorene():
    assert smiles_to_iupac("C1C2=CC=CC=C2C3=C1C4=CC=CC=C4C=C3") == "11H-benzo[a]fluorene"


def test_benzo_b_fluorene():
    assert smiles_to_iupac("C1C2=CC=CC=C2C3=CC4=CC=CC=C4C=C31") == "11H-benzo[b]fluorene"


def test_benzo_c_fluorene():
    assert smiles_to_iupac("C1C2=C(C3=CC=CC=C31)C4=CC=CC=C4C=C2") == "7H-benzo[c]fluorene"
