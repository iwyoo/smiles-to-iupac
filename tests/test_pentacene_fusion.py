from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_pentacene():
    # PubChem CID 67482 (CAS 239-98-5, "isohexaphene") -- InChI-verified
    # against the real compound.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5c(ccc6ccccc65)cc4cc3cc2c1") == "benzo[a]pentacene"


def test_plain_pentacene_still_resolves():
    # Regression check: this module's 26-atom/6-ring gate must not
    # shadow the existing 22-atom/5-ring plain pentacene recognition.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1") == "pentacene"


def test_benzo_b_pentacene_routes_to_hexacene():
    # Letter 'b' on pentacene is the linear chain extension -- the exact
    # same compound as the retained name "hexacene" `_aromatic.py` now
    # recognizes directly -- so this module excludes it.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1") == "hexacene"
