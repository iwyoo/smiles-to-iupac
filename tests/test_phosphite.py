import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Primary-source worked example (P-67.1.3.2, `tmp/bluebook/P6a.txt`
        # line ~4048-4049): `P(O-CH3)3` -> "trimethyl phosphite (PIN)".
        ("COP(OC)OC", "trimethyl phosphite"),
        ("CCOP(OCC)OCC", "triethyl phosphite"),
    ],
)
def test_phosphite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mirrors `_phosphate.py`'s own mixed-alkyl coverage (#878/#890) --
        # same shared `format_ester_words` mechanism, a different anion
        # word.
        ("COP(OCC)Oc1ccccc1", "ethyl methyl phenyl phosphite"),
        ("COP(OCC)OC", "ethyl dimethyl phosphite"),
    ],
)
def test_mixed_alkyl_phosphite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_partial_hydrogen_phosphite_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COP(OC)O")


def test_phosphate_ester_unaffected():
    # A P=O bond routes to `_phosphate.py` instead, unchanged.
    assert smiles_to_iupac("COP(=O)(OC)OC") == "trimethyl phosphate"


def test_phosphane_unaffected():
    assert smiles_to_iupac("CC[P](CC)CC") == "triethylphosphane"
