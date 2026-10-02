from smiles_to_iupac import smiles_to_iupac


def test_phosphindole():
    # PubChem CID 21811150.
    assert smiles_to_iupac("C1C=C2C=CC=CC2=P1") == "2H-phosphindole"


def test_isophosphindole():
    # PubChem CID 21667986.
    assert smiles_to_iupac("C1=CC2=CPC=C2C=C1") == "2H-isophosphindole"


def test_phosphinolizine():
    # PubChem CID 129643930.
    assert smiles_to_iupac("C1C=CC=C2P1C=CC=C2") == "4H-phosphinolizine"
