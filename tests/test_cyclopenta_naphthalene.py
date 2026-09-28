from smiles_to_iupac import smiles_to_iupac


def test_cyclopenta_a_naphthalene_1h():
    # PubChem CID 11745004.
    assert smiles_to_iupac("C1C=CC2=C1C3=CC=CC=C3C=C2") == "1H-cyclopenta[a]naphthalene"


def test_cyclopenta_a_naphthalene_3h():
    # PubChem CID 11105721.
    assert smiles_to_iupac("C1C=CC2=C1C=CC3=CC=CC=C32") == "3H-cyclopenta[a]naphthalene"


def test_cyclopenta_b_naphthalene_1h():
    # PubChem CID 6451436.
    assert smiles_to_iupac("C1C=CC2=CC3=CC=CC=C3C=C21") == "1H-cyclopenta[b]naphthalene"
