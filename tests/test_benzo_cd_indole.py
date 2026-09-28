from smiles_to_iupac import smiles_to_iupac


def test_benzo_cd_indole():
    # PubChem CID 22146199, IUPACName 'benzo[cd]indole', C11H7N.
    assert smiles_to_iupac("C1=CC2=C3C(=C1)C=NC3=CC=C2") == "benzo[cd]indole"
