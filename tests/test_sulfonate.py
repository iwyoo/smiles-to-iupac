import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-72.2.2.2.1.1: 'sulfonic acid' -> 'sulfonate'. PubChem PUG REST:
        # CID 85257 (methanesulfonate), CID 3717105 (ethanesulfonate),
        # CID 4431756 (propane-1-sulfonate).
        ("CS(=O)(=O)[O-]", "methanesulfonate"),
        ("CCS(=O)(=O)[O-]", "ethanesulfonate"),
        ("CCCS(=O)(=O)[O-]", "propane-1-sulfonate"),
    ],
)
def test_sulfonate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_sulfonate_with_double_bond():
    # PubChem CID 4447605: prop-2-ene-1-sulfonate.
    assert smiles_to_iupac("C=CCS(=O)(=O)[O-]") == "prop-2-ene-1-sulfonate"


def test_sulfonate_branched_r_group():
    # PubChem CID 498034/22184257: the sulfonate carbon need not be a
    # chain terminus, unlike `_carboxylate.py`/`_thioate.py`'s fixed-C1
    # carbon -- mirrors `_sulfonic_acid.py`'s own numbering.
    assert smiles_to_iupac("CC(C)S(=O)(=O)[O-]") == "propane-2-sulfonate"
    assert smiles_to_iupac("CCC(C)S(=O)(=O)[O-]") == "butane-2-sulfonate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-14.3.4.2(b): a two-carbon chain always omits the sulfonate's
        # own locant, regardless of substituent count or position --
        # PubChem CID 19003700 (1-chloroethanesulfonate), plus two more
        # confirmed live via PubChem PUG REST.
        ("CC(Cl)S(=O)(=O)[O-]", "1-chloroethanesulfonate"),
        ("ClCS(=O)(=O)[O-]", "chloromethanesulfonate"),
        ("ClCC(Cl)S(=O)(=O)[O-]", "1,2-dichloroethanesulfonate"),
    ],
)
def test_sulfonate_two_carbon_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), PubChem CID
        # 49867167.
        ("CC[C@@H](C)S(=O)(=O)[O-]", "(2R)-butane-2-sulfonate"),
        ("CC[C@H](C)S(=O)(=O)[O-]", "(2S)-butane-2-sulfonate"),
    ],
)
def test_sulfonate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_sulfonate_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)S(=O)(=O)[O-]") == "butane-2-sulfonate"


def test_multiple_sulfonate_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]S(=O)(=O)CCCS(=O)(=O)[O-]")


def test_sulfonate_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]S(=O)(=O)C1CCCCC1")


def test_sulfonate_carbon_in_double_bond_raises():
    # PubChem itself gives 'ethenesulfonate' for this shape, but it's out
    # of scope here (mirrors `_sulfonic_acid.py`'s identical
    # `_reject_enesulfonic_carbon` rejection).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CS(=O)(=O)[O-]")


def test_sulfonate_other_heteroatom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCCS(=O)(=O)[O-]")
