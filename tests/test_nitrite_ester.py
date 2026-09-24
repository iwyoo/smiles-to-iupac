import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-67.1.3.2) -- verified against each
        # structure's own PubChem-registered IUPAC name.
        ("CON=O", "methyl nitrite"),  # CID 12231
        ("CCON=O", "ethyl nitrite"),  # CID 8026
        ("CCCON=O", "propyl nitrite"),  # CID 10979
    ],
)
def test_nitrite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_fragments_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CON=O.C")


def test_two_nitrite_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=NOCCCON=O")


def test_nitro_unaffected():
    assert smiles_to_iupac("C[N+](=O)[O-]") == "nitromethane"


def test_nitrate_ester_unaffected():
    # A second, terminal oxygen routes to `_nitrate_ester.py` instead,
    # unchanged.
    assert smiles_to_iupac("CO[N+](=O)[O-]") == "methyl nitrate"
