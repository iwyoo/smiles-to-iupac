from smiles_to_iupac import smiles_to_iupac


def test_fluorene():
    # cross-checked against PubChem CID 6853 (9H-fluorene), C13H10.
    assert smiles_to_iupac("C1C2=CC=CC=C2C3=CC=CC=C31") == "9H-fluorene"
