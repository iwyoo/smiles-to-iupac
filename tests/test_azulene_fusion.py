from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_azulene():
    assert smiles_to_iupac("C1=CC=C2C=C3C=CC=CC3=C2C=C1") == "benzo[a]azulene"


def test_benzo_e_azulene():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC=C3C2=CC=C3") == "benzo[e]azulene"
