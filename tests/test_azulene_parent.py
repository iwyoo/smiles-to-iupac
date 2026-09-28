from smiles_to_iupac import smiles_to_iupac


def test_azulene():
    # cross-checked against PubChem CID 9231 (azulene), C10H8.
    assert smiles_to_iupac("C1=CC2=CC=CC=CC2=C1") == "azulene"
