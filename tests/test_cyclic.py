import pytest

from chemonym import smiles_to_iupac
from chemonym._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC1", "cyclopropane"),
        ("C1CCC1", "cyclobutane"),
        ("C1CCCCC1", "cyclohexane"),
        ("CC1CCCCC1", "methylcyclohexane"),
        ("CC1CCCCC1C", "1,2-dimethylcyclohexane"),
        ("CC1CC(C)CCC1", "1,3-dimethylcyclohexane"),
    ],
)
def test_smiles_to_iupac_cyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_raises():
    # bicyclo[2.2.1]heptane (norbornane): more than one ring.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CCC1CC2")


def test_ring_compound_substituent_raises():
    # sec-butyl-like branch on the ring: forks, needs P-29.4.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(CC)C1CCCCC1")
