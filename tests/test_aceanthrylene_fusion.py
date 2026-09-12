from smiles_to_iupac import smiles_to_iupac


def test_benzo_d_aceanthrylene():
    # PubChem CID 155608 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2c3ccccc3cc3cc4ccccc4c1c23") == "benzo[d]aceanthrylene"


def test_benzo_e_aceanthrylene():
    # PubChem CID 114910 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2c3ccccc3cc3c2c1cc1ccccc13") == "benzo[e]aceanthrylene"


def test_benzo_j_aceanthrylene():
    # PubChem CID 104987, "cholanthrylene" -- InChI-verified, and matches
    # PubChem's own computed IUPAC name exactly.
    assert smiles_to_iupac("C1=Cc2c3ccc4ccccc4c3cc3cccc1c23") == "benzo[j]aceanthrylene"


def test_benzo_k_aceanthrylene():
    # PubChem CID 146307 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2c3cc4ccccc4cc3cc3cccc1c23") == "benzo[k]aceanthrylene"


def test_benzo_l_aceanthrylene():
    # PubChem CID 105092 -- InChI-verified against the real compound.
    assert smiles_to_iupac("C1=Cc2c3c1cccc3cc1ccc3ccccc3c21") == "benzo[l]aceanthrylene"


def test_plain_aceanthrylene_still_resolves():
    # Regression check: this module's 20-atom/5-ring gate must not
    # shadow the existing 16-atom/4-ring plain aceanthrylene recognition.
    assert smiles_to_iupac("C1=CC=C2C3=C4C(=CC=CC4=CC2=C1)C=C3") == "aceanthrylene"


def test_benzo_a_aceanthrylene_routes_to_fluoranthene_base():
    # Letter 'a' on aceanthrylene is the exact same compound as
    # benzo[a]fluoranthene (PubChem CID 9146) -- already named via the
    # more established fluoranthene base, so this module excludes it.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC=C4C3=C2C5=CC=CC=C54") == "benzo[a]fluoranthene"
