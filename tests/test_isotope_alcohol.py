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
        # the -OH's own carbon still fixes locant 1 (P-44.4.1.8), so the
        # labeled carbon is unambiguously locant 2.
        ("[13CH3]CO", "ethan(2-13C)ol"),
        # The carbon isotope on the -OH-bearing carbon itself (locant 1),
        # still cited despite the suffix's own locant being omitted --
        # mirrors this project's established 'CC(Cl)S(=O)(=O)O' ->
        # '1-chloroethanesulfonic acid' precedent for an ordinary
        # substituent locant.
        ("C[13CH2]O", "ethan(1-13C)ol"),
        # Both a carbon isotope and a hydroxyl-oxygen isotope together,
        # combined in one descriptor group, carbon first (P-82.3.1's
        # alphabetical nuclide-symbol order, 'C' before 'O').
        ("C[13CH2][18OH]", "ethan(1-13C,18O)ol"),
        # A halogen substituent coexisting with the isotope descriptor.
        ("Cl[CH2][13CH2]O", "2-chloroethan(1-13C)ol"),
        # A 14C-labeled variant.
        ("[14CH3]CO", "ethan(2-14C)ol"),
    ],
)
def test_isotope_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
