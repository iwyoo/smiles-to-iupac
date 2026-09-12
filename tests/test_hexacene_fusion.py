from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_hexacene():
    # PubChem CID 23505249 -- InChI-verified against the real compound.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6c(ccc7ccccc76)cc5cc4cc3cc2c1") == "benzo[a]hexacene"


def test_plain_hexacene_still_resolves():
    # Regression check: this module's 30-atom/7-ring gate must not
    # shadow the existing 26-atom/6-ring plain hexacene recognition.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1") == "hexacene"


def test_benzo_b_hexacene_routes_to_heptacene():
    # Letter 'b' on hexacene is the linear chain extension -- the exact
    # same compound as the retained name "heptacene" `_aromatic.py`
    # already recognizes directly -- so this module excludes it.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6cc7ccccc7cc6cc5cc4cc3cc2c1") == "heptacene"
