from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_pentaphene():
    # PubChem CID 21087903.
    assert (
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=CC4=C(C=C32)C5=CC6=CC=CC=C6C=C5C=C4")
        == "benzo[a]pentaphene"
    )


def test_hexaphene():
    # Letter 'b' -- hexaphene (CID 123042) -- is a retained name, cited
    # as such rather than the systematic 'benzo[b]pentaphene'.
    assert (
        smiles_to_iupac("C1=CC=C2C=C3C(=CC2=C1)C=CC4=CC5=CC6=CC=CC=C6C=C5C=C43")
        == "hexaphene"
    )


def test_benzo_c_pentaphene():
    # PubChem CID 123040.
    assert (
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=C4C=CC5=CC6=CC=CC=C6C=C5C4=C3")
        == "benzo[c]pentaphene"
    )
