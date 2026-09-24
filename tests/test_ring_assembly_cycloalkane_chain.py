import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Primary source's own P-28.3.1 PIN worked example
        # (tmp/bluebook/P2.txt ~7904).
        ("C1CC1C1CC1C1CC1", "11,21:22,31-tercyclopropane"),
        ("C1CCC1C1CCC1C1CCC1", "11,21:22,31-tercyclobutane"),
        ("C1CCCC1C1CCCC1C1CCCC1", "11,21:22,31-tercyclopentane"),
        ("C1CCCCC1C1CCCCC1C1CCCCC1", "11,21:22,31-tercyclohexane"),
        ("C1CC1C1CC1C1CC1C1CC1", "11,21:22,31:32,41-quatercyclopropane"),
    ],
)
def test_ring_assembly_cycloalkane_chain_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituent():
    assert smiles_to_iupac("C1C(Cl)C1C1CC1C1CC1") == "12-chloro-11,21:22,31-tercyclopropane"


def test_unsaturated_ring_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCCC1C1=CCCCC1C1=CCCCC1")


def test_mixed_ring_kinds_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC1c1ccc(-c2ccccc2)cc1")


def test_mixed_ring_sizes_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC1C1CCCC1C1CC1")


def test_two_rings_routes_to_ring_assembly_module():
    # N=2 is `_ring_assembly.py`'s own job (P-28.2.1's primed-locant
    # scheme, not this module's composite-locant one) -- this just
    # confirms the N=3-6 module here no longer claims (and mis-fails on)
    # the N=2 shape.
    assert smiles_to_iupac("C1CC1C1CC1") == "1,1'-bi(cyclopropane)"
