import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methylpiperidine():
    # PubChem CID 12291 "1-methylpiperidine".
    assert smiles_to_iupac("CN1CCCCC1") == "1-methylpiperidine"


def test_ethylpiperidine():
    # PubChem CID 13007 "1-ethylpiperidine".
    assert smiles_to_iupac("CCN1CCCCC1") == "1-ethylpiperidine"


def test_methylpyrrolidine():
    # PubChem CID 8454 "1-methylpyrrolidine".
    assert smiles_to_iupac("CN1CCCC1") == "1-methylpyrrolidine"


def test_methylazepane():
    # 7-membered ring analogue of the piperidine/pyrrolidine cases above.
    assert smiles_to_iupac("CN1CCCCCC1") == "1-methylazepane"


def test_chloroethylpiperidine():
    # PubChem CID 74827 "1-(2-chloroethyl)piperidine".
    assert smiles_to_iupac("ClCCN1CCCCC1") == "1-(2-chloroethyl)piperidine"


def test_branched_n_substituent_is_parenthesized():
    # PubChem CID 240410's raw "1-propan-2-ylpiperidine", parenthesized
    # here per this project's usual compound-substituent convention.
    assert smiles_to_iupac("CC(C)N1CCCCC1") == "1-(propan-2-yl)piperidine"


def test_cyclopropyl_n_substituent():
    # PubChem CID 10909573 "1-cyclopropylpiperidine" -- a second,
    # unrelated ring inside the N-substituent itself must not be confused
    # with the amine ring.
    assert smiles_to_iupac("C1CC1N1CCCCC1") == "1-cyclopropylpiperidine"


def test_plain_ring_still_unaffected():
    # The unsubstituted ring itself (zero N-substituents) is a different,
    # already-working shape and must still route correctly.
    assert smiles_to_iupac("C1CCNCC1") == "piperidine"


def test_acyl_n_substituent_routes_to_hidden_amide():
    # An acyl N-substituent is `_hidden_amide_ketone.py`'s territory
    # (PR #409), not this module's -- must still route there correctly.
    assert smiles_to_iupac("CC(=O)N1CCCCC1") == "1-(piperidin-1-yl)ethan-1-one"


def test_two_heteroatom_ring_raises():
    # Morpholine (1,4-N,O) is a different axis (locant 4, not 1) -- out of
    # scope for this first slice.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN1CCOCC1")


def test_unsaturated_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CN1CCCCC1")


def test_substituted_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN1CCC(C)CC1")
