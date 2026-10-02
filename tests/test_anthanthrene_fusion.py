from smiles_to_iupac import smiles_to_iupac


def test_anthanthrene():
    # PubChem CID 9118.
    assert (
        smiles_to_iupac("C1=CC2=C3C(=C1)C=C4C=CC5=C6C4=C3C(=CC6=CC=C5)C=C2")
        == "dibenzo[def,mno]chrysene"
    )
