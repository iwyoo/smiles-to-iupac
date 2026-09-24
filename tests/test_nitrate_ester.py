import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-67.1.3.2), charged representation
        # (RDKit's own `[N+](=O)[O-]` canonical form) -- verified against
        # each structure's own PubChem-registered IUPAC name.
        ("CO[N+](=O)[O-]", "methyl nitrate"),  # CID 11724
        ("CCO[N+](=O)[O-]", "ethyl nitrate"),  # CID 12259
        ("CCCO[N+](=O)[O-]", "propyl nitrate"),  # CID 12307
        # The neutral N(=O)=O tautomer some SMILES writers use for the
        # same compound (mirrors `_nitro.py`'s own dual-representation
        # coverage).
        ("CON(=O)=O", "methyl nitrate"),
    ],
)
def test_nitrate_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_fragments_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CO[N+](=O)[O-].C")


def test_two_nitrate_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-][N+](=O)OCCCO[N+](=O)[O-]")


def test_nitro_unaffected():
    # A direct N-C bond (no ester oxygen) routes to `_nitro.py` instead,
    # unchanged.
    assert smiles_to_iupac("C[N+](=O)[O-]") == "nitromethane"
