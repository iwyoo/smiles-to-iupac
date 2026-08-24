import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
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


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOC1CCCCC1")


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=COOC")


def test_three_oxygens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOCOC")
