from smiles_to_iupac import smiles_to_iupac


def test_1h_cyclopenta_b_anthracene():
    # PubChem CID 21087838.
    assert smiles_to_iupac("C1C=CC2=CC3=CC4=CC=CC=C4C=C3C=C21") == "1H-cyclopenta[b]anthracene"


def test_2h_cyclopenta_b_anthracene():
    # PubChem CID 57983840 -- same fusion bond, the other indicated-H tautomer.
    assert smiles_to_iupac("C1C=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1") == "2H-cyclopenta[b]anthracene"


def test_1h_cyclopenta_a_anthracene():
    # PubChem CID 19845433.
    assert smiles_to_iupac("C1C=CC2=C1C3=CC4=CC=CC=C4C=C3C=C2") == "1H-cyclopenta[a]anthracene"
