from smiles_to_iupac import smiles_to_iupac


def test_pyrimidinone_resolves():
    assert smiles_to_iupac("O=c1nccc[nH]1") == "1H-pyrimidin-2-one"  # PubChem CID 68401


def test_plain_pyrimidine_unaffected():
    assert smiles_to_iupac("c1ccncn1") == "pyrimidine"


def test_pyridinone_still_resolves():
    assert smiles_to_iupac("O=c1cccc[nH]1") == "1H-pyridin-2-one"
