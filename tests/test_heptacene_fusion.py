from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_heptacene():
    # PubChem CID 21087897 -- a real registered structure (InChI-
    # verified), though no name synonym confirms "benzo[a]heptacene"
    # there (see module docstring: no methylheptacene isomer exists to
    # anchor against either, unlike the shorter acenes' own modules).
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6cc7c(ccc8ccccc87)cc6cc5cc4cc3cc2c1") == "benzo[a]heptacene"


def test_plain_heptacene_still_resolves():
    # Regression check: this module's 34-atom/8-ring gate must not
    # shadow the existing 30-atom/7-ring plain heptacene recognition.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6cc7ccccc7cc6cc5cc4cc3cc2c1") == "heptacene"


def test_benzo_b_heptacene_routes_to_octacene():
    # Letter 'b' on heptacene is the linear chain extension -- the exact
    # same compound as the retained name "octacene" `_aromatic.py`
    # already recognizes directly -- so this module excludes it.
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6cc7cc8ccccc8cc7cc6cc5cc4cc3cc2c1") == "octacene"
