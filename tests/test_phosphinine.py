from smiles_to_iupac import smiles_to_iupac


def test_phosphinine():
    # PubChem CID 123046.
    assert smiles_to_iupac("C1=CC=PC=C1") == "phosphinine"


def test_methylphosphinine():
    assert smiles_to_iupac("CC1=CC=CC=P1") == "2-methylphosphinine"


def test_phosphinoline():
    # PubChem CID 18624333.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC=P2") == "phosphinoline"


def test_isophosphinoline():
    # PubChem CID 136065.
    assert smiles_to_iupac("C1=CC=C2C=PC=CC2=C1") == "isophosphinoline"
