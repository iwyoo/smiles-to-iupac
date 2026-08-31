import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 'peroxy' substituent-prefix naming mirrors 'oxy' (_ether.py, see
        # test_ether.py: 'methoxymethane'/'1-methoxypropane', already
        # cross-checked against PubChem) with -O- replaced by -O-O-,
        # including the same P-14.3.4.2(b) omitted-locant convention for a
        # homogeneous two-carbon chain vs. the cited locant on a longer one.
        ("COOC", "methylperoxymethane"),
        ("CCOOCC", "ethylperoxyethane"),
        ("COOCC", "methylperoxyethane"),
        ("COOCCC", "1-methylperoxypropane"),
    ],
)
def test_peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_asymmetric_prefers_longer_chain_as_parent():
    assert smiles_to_iupac("CCCCOOC") == "1-methylperoxybutane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'peroxy' outside. Structure verified against PubChem: CID 22572410
    # ("1-propan-2-ylperoxybutane" -- PubChem elides the parentheses this
    # project's usual compound-substituent formatting keeps).
    assert smiles_to_iupac("CCCCOOC(C)C") == "1-(propan-2-yl)peroxybutane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOC1CCCCC1")


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=COOC")


def test_three_oxygens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOCOC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # The achiral parent structure is PubChem-verified (CID 19875133,
        # "2-ethylperoxybutane"); PubChem has no registered peroxide
        # stereoisomer CID, so these are a structural regression check on
        # the already-verified mechanism (same pattern as `_ether.py`'s
        # `winning_chain_from_carbon_graph` reuse, PR #224).
        ("CC[C@H](C)OOCC", "(2S)-2-ethylperoxybutane"),
        ("CC[C@@H](C)OOCC", "(2R)-2-ethylperoxybutane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CCC(C)OOCC") == "2-ethylperoxybutane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCOO[C@H](C)CC")
