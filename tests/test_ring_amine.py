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


def test_morpholine():
    # PubChem CID 7972 "4-methylmorpholine" -- the non-nitrogen heteroatom
    # (O) always wins locant 1, so the substituted nitrogen is locant 4.
    assert smiles_to_iupac("CN1CCOCC1") == "4-methylmorpholine"


def test_ethylmorpholine():
    # PubChem CID 7525 "4-ethylmorpholine".
    assert smiles_to_iupac("CCN1CCOCC1") == "4-ethylmorpholine"


def test_propan_2_yl_morpholine_is_parenthesized():
    # PubChem CID 73998's raw "4-propan-2-ylmorpholine", parenthesized
    # here per this project's usual compound-substituent convention.
    assert smiles_to_iupac("CC(C)N1CCOCC1") == "4-(propan-2-yl)morpholine"


def test_cyclopropylmorpholine():
    # PubChem CID 52140468 "4-cyclopropylmorpholine".
    assert smiles_to_iupac("C1CC1N1CCOCC1") == "4-cyclopropylmorpholine"


def test_thiomorpholine():
    # PubChem CID 523249 "4-methylthiomorpholine" -- S also outranks N for
    # locant 1 (P-22.2.1, O > S > N).
    assert smiles_to_iupac("CN1CCSCC1") == "4-methylthiomorpholine"


def test_piperazine():
    # PubChem CID 53167 "1-methylpiperazine" -- the symmetric N,N pair
    # gives the *substituted* nitrogen the lowest locant (1), not 4.
    assert smiles_to_iupac("CN1CCNCC1") == "1-methylpiperazine"


def test_ethylpiperazine():
    # PubChem CID 79196 "1-ethylpiperazine".
    assert smiles_to_iupac("CCN1CCNCC1") == "1-ethylpiperazine"


def test_plain_morpholine_and_piperazine_unaffected():
    # The unsubstituted rings themselves are a different, already-working
    # shape and must still route correctly.
    assert smiles_to_iupac("C1COCCN1") == "morpholine"
    assert smiles_to_iupac("C1CNCCN1") == "piperazine"


def test_unsaturated_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CN1CCCCC1")


def test_substituted_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN1CCC(C)CC1")


def test_methylsulfonylpiperidine():
    # PubChem CID 273952 "1-methylsulfonylpiperidine".
    assert smiles_to_iupac("CS(=O)(=O)N1CCCCC1") == "1-methylsulfonylpiperidine"


def test_ethylsulfonylpiperidine():
    # PubChem CID 3418537 "1-ethylsulfonylpiperidine".
    assert smiles_to_iupac("CCS(=O)(=O)N1CCCCC1") == "1-ethylsulfonylpiperidine"


def test_chloromethylsulfonylpiperidine_is_parenthesized():
    # PubChem CID 1519936 "1-(chloromethylsulfonyl)piperidine" -- the
    # whole 'Rsulfonyl' group is parenthesized because R ('chloromethyl')
    # is itself a compound substituent name.
    assert smiles_to_iupac("ClCS(=O)(=O)N1CCCCC1") == "1-(chloromethylsulfonyl)piperidine"


def test_methylsulfonylmorpholine():
    # PubChem CID 519344 "4-methylsulfonylmorpholine".
    assert smiles_to_iupac("CS(=O)(=O)N1CCOCC1") == "4-methylsulfonylmorpholine"


def test_methylsulfonylpiperazine():
    # PubChem CID 709161 "1-methylsulfonylpiperazine".
    assert smiles_to_iupac("CS(=O)(=O)N1CCNCC1") == "1-methylsulfonylpiperazine"


def test_methylsulfonylthiomorpholine():
    # PubChem CID 11309885 "4-methylsulfonylthiomorpholine".
    assert smiles_to_iupac("CS(=O)(=O)N1CCSCC1") == "4-methylsulfonylthiomorpholine"


def test_unrelated_sulfonamide_unaffected():
    # A plain aromatic sulfonamide (no ring nitrogen at all) is a
    # different, already-working shape and must still route correctly.
    assert smiles_to_iupac("NS(=O)(=O)c1ccccc1") == "benzenesulfonamide"
