import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("[2H]CC", id="ethane_deuterium_raises"),
        pytest.param("[3H]C", id="tritium_raises"),
    ],
)
def test_ethane_deuterium_raises_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[2H]C([2H])(Cl)Cl", "dichloro(2H2)methane", id="dichlorodideuteromethane_name"),
        pytest.param("[14CH4]", "(14C)methane", id="carbon_14_methane_name"),
        pytest.param("[2H][13CH3]", "(13C,2H1)methane", id="deuterium_and_carbon_isotope_together_methane_name"),
        pytest.param("C[14CH2]C([2H])C", "(2-14C,3-2H1)butane", id="deuterium_and_carbon_isotope_together_chain_name"),
    ],
)
def test_dichlorodideuteromethane_name_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotopically_labeled_halogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]C([37Cl])")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[14CH2]CC", "(2-14C)butane", id="carbon_14_butane_name"),
        pytest.param("FC(F)(F)C[2H]", "1,1,1-trifluoro(2-2H1)ethane", id="trifluoro_deuterio_ethane_name"),
    ],
)
def test_carbon_14_butane_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("FCC[2H]", id="single_halogen_ethane_deuterium_raises"),
        pytest.param("[13CH3]C[14CH2]C", id="mixed_carbon_isotope_nuclides_raises"),
    ],
)
def test_single_halogen_ethane_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[13CH2]O", "(1-13C)ethan-1-ol"),
        ("[13CH3]O", "(13C)methanol"),
    ],
)
def test_isotope_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C[13CH2][18OH]", id="combined_carbon_and_oxygen_isotope_raises"),
        pytest.param("OCCC[18OH]", id="multiple_hydroxyls_raises"),
        pytest.param("[18OH]CCC", id="longer_chain_raises"),
        pytest.param("[18OH]C=C", id="unsaturation_raises"),
    ],
)
def test_isotope_alcohol_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("[13CH3][C@H](Cl)C(=O)O", id="specified_stereocenter_alongside_isotope_carboxylic_acid_raises"),
    ],
)
def test_specified_stereocenter_alongside_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[18O])C", "propan-2-(18O)one"),
        ("CC(=O)[14CH3]", "(3-14C)propan-2-one"),
    ],
)
def test_isotope_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(=[18O])[13CH3]", id="both_isotope_kinds_together_raises"),
        pytest.param("O=[13C]1CCCCC1", id="ring_alongside_isotope_ketone_raises"),
        pytest.param("[13CH3]C(=O)[C@H](Cl)C", id="specified_stereocenter_alongside_isotope_ketone_raises"),
    ],
)
def test_both_isotope_kinds_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("[13CH]([2H])1CCCCC1", id="same_position_carbon_and_deuterium_raises"),
        pytest.param("[3H]C1CCCCC1", id="tritium_raises"),
        pytest.param("[2H]C1CCCCC1C", id="substituent_branch_raises"),
    ],
)
def test_isotope_ring_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC=[18O]", "(18O)formic acid"),
    ],
)
def test_isotope_descriptor_before_retained_acid_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[13CH3][14CH2]O", "(2-13C,1-14C)ethanol", id="two_carbon_nuclides_in_one_series"),
        pytest.param("CC(=O)O[13CH3]", "(13C)methyl acetate", id="isotopic_methyl_group_of_an_ester"),
        pytest.param("O=C(O[13CH2]C)c1ccccc1", "(1-13C)ethyl benzoate", id="isotopic_alkyl_chain_of_an_ester"),
        pytest.param("O=C(c1ccccc1)[13CH3]", "1-phenyl(2-13C)ethan-1-one", id="isotope_in_parent_after_substituent_prefix"),
        pytest.param("Cc1cccnc1[13CH3]", "2-(13C)methyl-3-methylpyridine", id="isotopic_substituent_cited_before_unmodified_one"),
        pytest.param("CC(=O)Nc1ccc([131I])cc1", "N-[4-(131I)iodophenyl]acetamide", id="isotopic_halogen_inside_a_ring_substituent"),
    ],
)
def test_isotopic_descriptor_follows_the_unmodified_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles", ["CC(=O)Nc1ccccc1C[13CH3]", "CC(=O)c1ccc([18F])cc1"])
def test_unplaceable_isotope_label_is_never_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
