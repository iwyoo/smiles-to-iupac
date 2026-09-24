from smiles_to_iupac import smiles_to_iupac


def test_thiophosphoric_acid():
    # Real PubChem structure, CID 167254.
    assert smiles_to_iupac("OP(=S)(O)O") == "phosphorothioic O,O,O-acid"
