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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-69.3's own worked example (`tmp/bluebook/P6a.txt` lines
        # 8859-8863): `[MgMe]I` -> 'methylmagnesium iodide (PIN,
        # compositional name)'.
        ("C[Mg]I", "methylmagnesium iodide"),
        ("CC[Mg]Br", "ethylmagnesium bromide"),
        ("[Mg](Cl)c1ccccc1", "phenylmagnesium chloride"),
        ("CC[Ca]Cl", "ethylcalcium chloride"),
    ],
)
def test_grignard_rmx_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_inorganic_halide_of_group1_metal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Li]Cl")


def test_inorganic_dihalide_of_group2_metal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cl[Mg]Cl")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Mg]CC", "ethyl(methyl)magnesium"),
        ("C[Mg](C)C", "trimethylmagnesium"),
        ("C[Mg]", "methylmagnesium"),
    ],
)
def test_group2_additive_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_metal_atoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Li]C.[Na]CC")


def test_group1_metal_with_substituted_alkyl():
    assert smiles_to_iupac("[Li]CO") == "hydroxymethyllithium"


def test_group13_hydride_unaffected():
    assert smiles_to_iupac("CC[Al](CC)CC") == "triethylalumane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[BeH]", "ethylhydridoberyllium"),
        ("C[Be]C", "dimethylberyllium"),
    ],
)
def test_beryllium_additive_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Li][CH3][Li]", "μ-methyl-dilithium"),
        ("[Li]C[Li]", "μ-methanediyl-dilithium"),
    ],
)
def test_carbon_bridged_dilithium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
