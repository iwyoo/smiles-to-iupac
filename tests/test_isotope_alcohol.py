import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-82.3.1's own confirmed worked example.
        ("[18OH]C", "methan(18O)ol"),
        # P-82.6.1.1's own confirmed 2-carbon deuterium worked example
        # ('ethan(2H)ol'), same shape with 18O instead.
        ("CC[18OH]", "ethan(18O)ol"),
        # A skeletal carbon isotope on the non-suffix carbon of ethanol:
        # front-of-name, with the otherwise-omitted -OH locant restored,
        # per the Blue Book's own confirmed '(2-13C)ethan-1-ol (PIN)'
        # (#804 correction -- an earlier pass of this module wrongly
        # spliced this before 'ol' instead).
        ("[13CH3]CO", "(2-13C)ethan-1-ol"),
        # The carbon isotope on the -OH-bearing carbon itself (locant 1),
        # same front-of-name placement.
        ("C[13CH2]O", "(1-13C)ethan-1-ol"),
        # A halogen substituent coexisting with the isotope descriptor --
        # the front-of-name descriptor sits ahead of the substituent
        # prefix too.
        ("Cl[CH2][13CH2]O", "(1-13C)2-chloroethan-1-ol"),
        # A 14C-labeled variant.
        ("[14CH3]CO", "(2-14C)ethan-1-ol"),
        # A mononuclear (methanol) carbon isotope never cites a locant,
        # mirroring P-82.2.1's own '(14C)methane'.
        ("[13CH3]O", "(13C)methanol"),
    ],
)
def test_isotope_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotope_alcohol_combined_carbon_and_oxygen_isotope_raises():
    # No confirmed worked example for this combination once an internal
    # locant is possible (see module docstring) -- an earlier pass of
    # this module guessed a merged-descriptor format here; #804 retired
    # that guess rather than keep it unverified.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[13CH2][18OH]")


def test_plain_alcohol_unaffected():
    assert smiles_to_iupac("OC") == "methanol"
    assert smiles_to_iupac("OCC") == "ethanol"


def test_isotope_alcohol_multiple_hydroxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCCC[18OH]")


def test_isotope_alcohol_longer_chain_raises():
    # A chain length of 3+ has no confirmed worked example for where the
    # isotope descriptor splices in relative to the -OH's own cited
    # locant (P-44.4.1.8), so this is deferred rather than guessed at.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[18OH]CCC")


def test_isotope_alcohol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[18OH]C=C")


def test_isotope_alcohol_mixed_carbon_isotope_nuclides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3][14CH2]O")


def test_isotope_alcohol_17o_name():
    assert smiles_to_iupac("[17OH]C") == "methan(17O)ol"


def test_isotope_alcohol_unsupported_oxygen_isotope_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[15OH]C")
