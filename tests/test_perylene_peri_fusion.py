from smiles_to_iupac import smiles_to_iupac


def test_benzo_ghi_perylene():
    # PubChem CID 9117.
    assert smiles_to_iupac("C1=CC2=C3C(=C1)C4=CC=CC5=C4C6=C(C=C5)C=CC(=C36)C=C2") == "benzo[ghi]perylene"
