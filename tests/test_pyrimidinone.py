from smiles_to_iupac import smiles_to_iupac


def test_pyrimidinone_resolves():
    assert smiles_to_iupac("O=c1nccc[nH]1") == "pyrimidin-2(1H)-one"  # PubChem CID 68401


def test_plain_pyrimidine_unaffected():
    assert smiles_to_iupac("c1ccncn1") == "pyrimidine"


def test_pyridinone_still_resolves():
    assert smiles_to_iupac("O=c1cccc[nH]1") == "pyridin-2(1H)-one"
