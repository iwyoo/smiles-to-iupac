import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


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
    # adamantane: tricyclic, more than one ring and out of scope even after
    # bicyclic support (see _bicyclic.py, test_bicyclic.py).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1C2CC3CC1CC(C2)C3")


def test_ring_compound_substituent():
    # sec-butyl-like branch on the ring (P-29.4); the only substituent on an
    # otherwise unsubstituted ring, so its locant is omitted (P-14.3.3).
    assert smiles_to_iupac("CC(CC)C1CCCCC1") == "(1-methylpropyl)cyclohexane"
