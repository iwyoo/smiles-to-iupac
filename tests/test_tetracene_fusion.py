from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_tetracene():
    # PubChem CID 67470 -- InChI-verified against the real compound, and
    # matches PubChem's own IUPACName property exactly.
    assert smiles_to_iupac("c1ccc2cc3cc4c(ccc5ccccc54)cc3cc2c1") == "benzo[a]tetracene"


def test_plain_tetracene_still_resolves():
    # Regression check: this module's 22-atom/5-ring gate must not
    # shadow the existing 18-atom/4-ring plain tetracene recognition.
    assert smiles_to_iupac("c1ccc2cc3cc4ccccc4cc3cc2c1") == "tetracene"


def test_benzo_b_tetracene_routes_to_pentacene():
    # Letter 'b' on tetracene is the linear chain extension -- the exact
    # same compound as the retained name "pentacene" `_aromatic.py`
    # already recognizes directly -- so this module excludes it.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1") == "pentacene"
