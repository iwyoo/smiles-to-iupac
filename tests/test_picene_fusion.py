from smiles_to_iupac import smiles_to_iupac


def test_benzo_b_picene():
    # PubChem CID 123038.
    assert (
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC6=CC=CC=C6C=C54")
        == "benzo[b]picene"
    )


def test_benzo_c_picene():
    # PubChem CID 9168.
    assert (
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=C4C=CC6=CC=CC=C65")
        == "benzo[c]picene"
    )
