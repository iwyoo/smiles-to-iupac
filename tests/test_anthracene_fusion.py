from smiles_to_iupac import smiles_to_iupac


def test_benzo_a_anthracene():
    # PubChem CID 5954 -- also the Blue Book's own worked example for
    # P-25.3.1.3's fusion-locant-letter mechanism.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32") == "benzo[a]anthracene"


def test_linear_fusion_falls_through_to_tetracene():
    # The linear (letter 'b') fusion is the same compound as tetracene,
    # already a retained name -- this module must defer to `_aromatic.py`
    # rather than emit 'benzo[b]anthracene'.
    assert smiles_to_iupac("C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1") == "tetracene"


def test_angular_four_ring_chrysene_not_claimed_here():
    # A different tetracyclic angular topology (chrysene, CID 9171) --
    # not an anthracene+benzo shape at all, so this module must not claim
    # it. It does resolve correctly, via `_phenanthrene_fusion.py`'s
    # letter 'a' instead (see that module's own tests).
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43") == "chrysene"
