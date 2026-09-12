from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_acephenanthrylene():
    # PubChem CID 11482160 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2cc3ccccc3c3c2c1cc1ccccc13") == "benzo[a]acephenanthrylene"


def test_benzo_j_acephenanthrylene():
    # PubChem CID 159608 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2cc3c4ccccc4ccc3c3cccc1c23") == "benzo[j]acephenanthrylene"


def test_benzo_k_acephenanthrylene():
    # PubChem CID 159607 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2cc3cc4ccccc4cc3c3cccc1c23") == "benzo[k]acephenanthrylene"


def test_benzo_l_acephenanthrylene():
    # PubChem CID 10848420 (indexed under the plain, unlettered synonym
    # "benzoacephenanthrylene") -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2cc3ccc4ccccc4c3c3cccc1c23") == "benzo[l]acephenanthrylene"


def test_plain_acephenanthrylene_still_resolves():
    # Regression check: this module's 20-atom/5-ring gate must not
    # shadow the existing 16-atom/4-ring plain acephenanthrylene
    # recognition.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC4=C3C2=CC=C4") == "acephenanthrylene"


def test_benzo_b_acephenanthrylene_routes_to_aceanthrylene_base():
    # Letter 'b' on acephenanthrylene is the exact same compound as
    # benzo[e]aceanthrylene (PubChem CID 114910) -- already named via the
    # more established aceanthrylene base, so this module excludes it.
    assert smiles_to_iupac("C1=Cc2c3ccccc3cc3c2c1cc1ccccc13") == "benzo[e]aceanthrylene"


def test_benzo_e_acephenanthrylene_routes_to_fluoranthene_base():
    # Letter 'e' on acephenanthrylene is the exact same compound as
    # benzo[b]fluoranthene (PubChem CID 9153) -- already named via the
    # more established fluoranthene base, so this module excludes it.
    assert smiles_to_iupac("c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13") == "benzo[b]fluoranthene"
