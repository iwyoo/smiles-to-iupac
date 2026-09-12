from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_fluoranthene():
    # PubChem CID 9146.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC=C4C3=C2C5=CC=CC=C54") == "benzo[a]fluoranthene"


def test_benzo_b_fluoranthene():
    # PubChem CID 9153.
    assert smiles_to_iupac("C1=CC=C2C3=C4C(=CC=C3)C5=CC=CC=C5C4=CC2=C1") == "benzo[b]fluoranthene"


def test_benzo_j_fluoranthene():
    # PubChem CID 9152 (this is also the compound PubChem resolves the
    # name 'benzo[l]fluoranthene' to -- 'j' is the alphabetically lower,
    # winning letter per P-25.3.1.3's own tie-break).
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC5=C4C3=CC=C5") == "benzo[j]fluoranthene"


def test_benzo_k_fluoranthene():
    # PubChem CID 9158.
    assert smiles_to_iupac("C1=CC=C2C=C3C4=CC=CC5=C4C(=CC=C5)C3=CC2=C1") == "benzo[k]fluoranthene"


def test_plain_fluoranthene_still_resolves():
    # Regression check: this module's 20-atom/5-ring gate must not
    # shadow the existing 16-atom/4-ring plain fluoranthene recognition.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C3=CC=CC4=C3C2=CC=C4") == "fluoranthene"
