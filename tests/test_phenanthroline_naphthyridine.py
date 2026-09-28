from smiles_to_iupac import smiles_to_iupac


def test_1_7_phenanthroline():
    # PubChem CID 67473's canonical SMILES.
    assert smiles_to_iupac("C1=CC2=C(C3=C(C=C2)N=CC=C3)N=C1") == "1,7-phenanthroline"


def test_1_10_phenanthroline():
    # PubChem CID 1318's canonical SMILES.
    assert smiles_to_iupac("C1=CC2=C(C3=C(C=CC=N3)C=C2)N=C1") == "1,10-phenanthroline"


def test_4_7_phenanthroline():
    # PubChem CID 67472's canonical SMILES.
    assert smiles_to_iupac("C1=CC2=C(C=CC3=C2C=CC=N3)N=C1") == "4,7-phenanthroline"


def test_1_5_naphthyridine():
    # PubChem CID 136070's canonical SMILES.
    assert smiles_to_iupac("C1=CC2=C(C=CC=N2)N=C1") == "1,5-naphthyridine"


def test_1_8_naphthyridine():
    # PubChem CID 136069's canonical SMILES.
    assert smiles_to_iupac("C1=CC2=C(N=C1)N=CC=C2") == "1,8-naphthyridine"
