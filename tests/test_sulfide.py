import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_sulfide():
    # 'sulfanyl' substituent-prefix naming mirrors 'oxy' (_ether.py, see
    # test_ether.py, already cross-checked against PubChem) with -O-
    # replaced by -S-, same P-14.3.4.2(b) omitted-locant convention.
    assert smiles_to_iupac("CSC") == "methylsulfanylmethane"


def test_methylsulfanylpropane():
    assert smiles_to_iupac("CSCCC") == "1-methylsulfanylpropane"


def test_diethyl_sulfide():
    assert smiles_to_iupac("CCSCC") == "ethylsulfanylethane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'sulfanyl' outside. Structure verified against PubChem: CID 522478
    # ("1-propan-2-ylsulfanylbutane" -- PubChem elides the parentheses this
    # project's usual compound-substituent formatting keeps).
    assert smiles_to_iupac("CCCCSC(C)C") == "1-(propan-2-yl)sulfanylbutane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Same tied-parent-selection cases as `test_ether.py`, sulfanyl in
        # place of oxy -- see that module's tests for the full reasoning.
        ("CC(C)SC(C)C", "2-(propan-2-yl)sulfanylpropane"),
        ("CC(C)CSC(C)(C)C", "1-tert-butylsulfanyl-2-methylpropane"),
        ("CC(C)(C)CSCCC(C)C", "1-(2,2-dimethylpropyl)sulfanyl-3-methylbutane"),
        ("CC(C)CSC(C)CC", "2-(2-methylpropyl)sulfanylbutane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified (CID 92859631/92174457): the parent (longer)
        # chain carries the stereocenter, so the prefix mirrors
        # `_ether.py`'s P-91.3 mechanism via `winning_chain_from_carbon_graph`.
        ("CC[C@H](C)SCC", "(2S)-2-ethylsulfanylbutane"),
        ("CC[C@@H](C)SCC", "(2R)-2-ethylsulfanylbutane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CCC(C)SCC") == "2-ethylsulfanylbutane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCS[C@H](C)CC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Direct ring-sulfur bond. Structure PubChem-confirmed: CID 12144
        # 'c1ccccc1SCC' -> 'ethylsulfanylbenzene', CID 76397
        # 'c1ccccc1SC(C)C' -> 'propan-2-ylsulfanylbenzene'.
        ("c1ccccc1SCC", "ethylsulfanylbenzene"),
        ("c1ccccc1SC", "methylsulfanylbenzene"),
        # Direct ring-sulfur bond with a branched R': the same
        # '(...)sulfanyl' parenthesization as the plain two-chain path
        # above (P-63.2.2.1.1's worked example), even though PubChem's own
        # auto-generated name for the same structure omits the
        # parentheses ('c1ccccc1SC(C)C' -> 'propan-2-ylsulfanylbenzene').
        ("c1ccccc1SC(C)C", "(propan-2-yl)sulfanylbenzene"),
        # Chain spacer between the ring and the sulfide sulfur -- the
        # whole branch is parenthesized when compound, matching this
        # project's existing benzene-ring-chain convention (`_ether.py`),
        # rather than PubChem's own un-parenthesized auto-names
        # ('ethylsulfanylmethylbenzene', CID 80431).
        ("c1ccccc1CSCC", "(ethylsulfanylmethyl)benzene"),
        # A branched R' behind a chain spacer: the '(...)sulfanyl' term
        # already carries round brackets, so the outer compound-branch
        # wrap escalates to square brackets instead of nesting round ones
        # (mirrors `_ether.py`'s own bracket-escalation example).
        # Structure PubChem-confirmed ('c1ccccc1CSC(C)C' ->
        # 'propan-2-ylsulfanylmethylbenzene', CID 525513).
        ("c1ccccc1CSC(C)C", "[(propan-2-yl)sulfanylmethyl]benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CSc1ccccc1SC")


def test_benzene_ring_stereocenter_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1S[C@H](C)CC")
