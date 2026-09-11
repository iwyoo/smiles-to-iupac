from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_pyrene():
    # PubChem CID 2336 -- the well-known environmental carcinogen.
    assert smiles_to_iupac("C1=CC=C2C3=C4C(=CC2=C1)C=CC5=C4C(=CC=C5)C=C3") == "benzo[a]pyrene"


def test_benzo_e_pyrene():
    # PubChem CID 9128.
    assert smiles_to_iupac("c1ccc2c(c1)c1cccc3ccc4cccc2c4c31") == "benzo[e]pyrene"


def test_plain_pyrene_still_resolves():
    # Regression check: this module's 20-atom/5-ring gate must not
    # shadow the existing 16-atom/4-ring plain pyrene recognition.
    assert smiles_to_iupac("c1cc2ccc3cccc4ccc(c1)c2c34") == "pyrene"
