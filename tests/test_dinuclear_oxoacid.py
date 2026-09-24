from smiles_to_iupac import smiles_to_iupac


def test_diphosphoric_acid():
    # Real PubChem structure, CID 1023 (pyrophosphoric acid).
    assert smiles_to_iupac("OP(=O)(O)OP(=O)(O)O") == "diphosphoric acid"
