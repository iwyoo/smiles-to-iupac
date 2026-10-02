from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_perylene():
    # PubChem CID 115240.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC=C4C3=C2C5=CC=CC6=C5C4=CC=C6") == "benzo[a]perylene"


def test_benzo_b_perylene():
    # PubChem CID 67455.
    assert smiles_to_iupac("C1=CC=C2C3=C4C(=CC=C3)C5=CC=CC6=C5C(=CC=C6)C4=CC2=C1") == "benzo[b]perylene"
