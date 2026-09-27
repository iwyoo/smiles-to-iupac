import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CS(=O)C=C1", "thiophene 1-oxide"),  # PubChem CID 9548690
        ("O=S1CCCCC1", "thiane 1-oxide"),  # PubChem CID 534965
    ],
)
def test_hetero_ring_oxide_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_ring_unaffected():
    assert smiles_to_iupac("c1ccsc1") == "thiophene"


def test_two_heteroatom_ring_oxide_still_raises():
    with pytest.raises(Exception):
        smiles_to_iupac("O=S1C=CC=CN1")
