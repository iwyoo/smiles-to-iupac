import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_ethane_deuterium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]CC")


def test_tritium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[3H]C")


def test_dichlorodideuteromethane_name():
    assert smiles_to_iupac("[2H]C([2H])(Cl)Cl") == "dichloro(2H2)methane"


def test_carbon_14_methane_name():
    assert smiles_to_iupac("[14CH4]") == "(14C)methane"


def test_deuterium_and_carbon_isotope_together_methane_name():
    assert smiles_to_iupac("[2H][13CH3]") == "(13C,2H1)methane"


def test_deuterium_and_carbon_isotope_together_chain_name():
    assert smiles_to_iupac("C[14CH2]C([2H])C") == "(2-14C,3-2H1)butane"


def test_isotopically_labeled_halogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]C([37Cl])")


def test_carbon_14_butane_name():
    assert smiles_to_iupac("C[14CH2]CC") == "(2-14C)butane"


def test_trifluoro_deuterio_ethane_name():
    # Blue Book P-82.6.2 worked example: "1,1,1-trifluoro(2-2H1)ethane (PIN)".
    assert smiles_to_iupac("FC(F)(F)C[2H]") == "1,1,1-trifluoro(2-2H1)ethane"


def test_single_halogen_ethane_deuterium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("FCC[2H]")


def test_mixed_carbon_isotope_nuclides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3]C[14CH2]C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[13CH2]O", "(1-13C)ethan-1-ol"),
        ("[13CH3]O", "(13C)methanol"),
    ],
)
def test_isotope_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotope_alcohol_combined_carbon_and_oxygen_isotope_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[13CH2][18OH]")


def test_isotope_alcohol_multiple_hydroxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCCC[18OH]")


def test_isotope_alcohol_longer_chain_raises():
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[18O])O", "(18O)acetic acid"),
    ],
)
def test_isotope_carboxylic_acid_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_chain_unsaturation_alongside_isotope_carboxylic_acid():
    assert smiles_to_iupac("[13CH3]C(=O)OC=C") == "ethenyl (2-13C)acetate"


def test_specified_stereocenter_alongside_isotope_carboxylic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3][C@H](Cl)C(=O)O")


def test_mixed_carbon_isotope_nuclides_raises__isotope_carboxylic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3][14CH2]C(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[18O])C", "propan-2-(18O)one"),
        ("CC(=O)[14CH3]", "(3-14C)propan-2-one"),
    ],
)
def test_isotope_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_both_isotope_kinds_together_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=[18O])[13CH3]")


def test_ring_alongside_isotope_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[13C]1CCCCC1")


def test_specified_stereocenter_alongside_isotope_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3]C(=O)[C@H](Cl)C")


def test_mixed_carbon_isotope_nuclides_raises__isotope_ketone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH3]C(=O)[14CH3]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[2H]C1CC1", "(2H1)cyclopropane"),
        ("[13CH]1CCCCC1", "(13C)cyclohexane"),
        ("[13CH]1CC([2H])CCC1", "(1-13C,3-2H1)cyclohexane"),
    ],
)
def test_isotope_ring_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotope_ring_same_position_carbon_and_deuterium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[13CH]([2H])1CCCCC1")


def test_isotope_ring_tritium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[3H]C1CCCCC1")


def test_isotope_ring_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]C1CCCCC1C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC=[18O]", "(18O)formic acid"),
    ],
)
def test_isotope_descriptor_before_retained_acid_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
