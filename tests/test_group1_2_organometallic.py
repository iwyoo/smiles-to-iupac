import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyllithium_blue_book_worked_example():
    # Blue Book P-69.3 (`tmp/bluebook/P6a.txt` line ~8843-8845): '[LiMe]'
    # -> 'methyllithium (PIN)'. PubChem is not usable to verify this
    # shape (see module docstring) - its own registered entries for these
    # compounds are fully-dissociated ionic representations, not the
    # covalent single-molecule structure the Blue Book itself names.
    assert smiles_to_iupac("[Li]C") == "methyllithium"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Na]CC", "ethylsodium"),
        ("[K]CC", "ethylpotassium"),
        ("[Li]CCC", "propyllithium"),
        ("c1ccccc1[Li]", "phenyllithium"),
    ],
)
def test_group1_organometallic_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[Mg]CC", "diethylmagnesium"),
        ("C[Mg]C", "dimethylmagnesium"),
        ("CC[Ca]CC", "diethylcalcium"),
        ("C[Ca]C", "dimethylcalcium"),
    ],
)
def test_group2_organometallic_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_substituents_on_group2_metal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Mg]CC")


def test_single_substituent_group2_metal_raises():
    # A neutral Group 2 metal with only 1 organic substituent and no
    # halide would be an open-shell radical, not a real closed-shell
    # compound -- the Grignard-shaped R-M-X case (1 organic + 1 halide)
    # is a separate milestone step, not this one.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Mg]")


def test_three_substituents_on_group2_metal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Mg](C)C")


def test_two_metal_atoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Li]C.[Na]CC")


def test_other_heteroatom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Li]CO")


def test_group13_hydride_unaffected():
    assert smiles_to_iupac("CC[Al](CC)CC") == "triethylalumane"
