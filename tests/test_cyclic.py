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
    # pentagonal prism: a hexacyclic ring system (cyclomatic number 6, ten
    # branch atoms of degree 3, two pentagons joined by five bridging bonds):
    # out of scope for bicyclic through pentacyclic support alike (see
    # _bicyclic.py, _polycyclic.py). Built from scratch via RDKit's RWMol
    # (two independent 5-cycles plus one bond between each corresponding
    # pair of atoms), not copied from a database.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C1C2C2C3C4C12")


def test_ring_compound_substituent():
    # sec-butyl-like branch on the ring (P-29.4); the only substituent on an
    # otherwise unsubstituted ring, so its locant is omitted (P-14.3.3).
    assert smiles_to_iupac("CC(CC)C1CCCCC1") == "(1-methylpropyl)cyclohexane"
