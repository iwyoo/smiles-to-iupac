import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single deuterium on an otherwise-plain, fully symmetric
        # monocyclic ring omits the locant entirely (P-82.6.1.1: every
        # ring position is equivalent to every other before any
        # modification, unlike a 3+-carbon chain's own chemically
        # distinct interior/terminal positions).
        ("[2H]C1CCCCC1", "(2H1)cyclohexane"),
        ("[2H]C1CCCC1", "(2H1)cyclopentane"),
        ("[2H]C1CCCCCC1", "(2H1)cycloheptane"),
        ("[2H]C1CC1", "(2H1)cyclopropane"),
        # A single skeletal carbon isotope on a plain ring: same
        # locant-omission rule.
        ("[13CH]1CCCCC1", "(13C)cyclohexane"),
        ("[14CH]1CCCC1", "(14C)cyclopentane"),
        # Two deuteriums on the same ring carbon: still one modified
        # position, still omitted.
        ("[2H]C1([2H])CCCCC1", "(2H2)cyclohexane"),
        # Two deuteriums at different ring positions: the modification
        # breaks the ring's symmetry, so locants are now necessary and
        # are minimized together (P-82.5.2's combined-locant-set rule).
        ("[2H]C1CCCC([2H])C1", "(1,3-2H2)cyclohexane"),
        # A halogen substituent alongside a single deuterium: the
        # halogen's own locant forces citation even though the
        # deuterium alone would have omitted it.
        ("ClC1CCCCC1[2H]", "1-chloro(2-2H1)cyclohexane"),
        # A skeletal carbon isotope and deuterium at two distinct ring
        # positions: both cited (two distinct modified positions).
        ("[13CH]1CC([2H])CCC1", "(1-13C,3-2H1)cyclohexane"),
    ],
)
def test_isotope_ring_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotope_ring_same_position_carbon_and_deuterium_raises():
    # A skeletal carbon isotope and deuterium both modifying the exact
    # same single ring position, with no halogen -- unconfirmed whether
    # the locant is omitted here the same way a single-kind modification
    # would be.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH]([2H])1CCCCC1")


def test_isotope_ring_tritium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[3H]C1CCCCC1")


def test_isotope_ring_substituent_branch_raises():
    # A plain alkyl branch hanging off the ring isn't one of the allowed
    # deuterium/halogen/carbon-isotope substituents.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]C1CCCCC1C")


def test_isotope_ring_mixed_carbon_isotope_nuclides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH]1C[14CH]CCC1")
