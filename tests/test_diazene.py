import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 123195: unsubstituted parent.
        ("N=N", "diazene"),
        # PubChem CID 123421: single substituent.
        ("CN=N", "methyldiazene"),
        # PubChem CID 10421: symmetric disubstitution.
        ("CN=NC", "dimethyldiazene"),
        # PubChem CID 526060: asymmetric disubstitution -- alphabetically
        # first ('ethyl') unparenthesized, the other ('methyl')
        # parenthesized (P-16.5.1.3.1, mirrors _phosphane.py's identical
        # rule, already verified there for 'ethyl(methyl)phosphane').
        ("CN=NCC", "ethyl(methyl)diazene"),
        # PubChem CID 13183.
        ("CCN=NCC", "diethyldiazene"),
    ],
)
def test_diazene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)N=NC")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(CC1)N=NC")


def test_aromatic_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1N=NC")
