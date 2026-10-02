import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Ti](Cl)(Cl)Cl", "trichlorido(methyl)titanium"),
        ("CC[Ti](C)(Cl)Cl", "dichlorido(ethyl)(methyl)titanium"),
        ("C[Hg]C", "dimethylmercury"),
        ("C[Hg]Cl", "chlorido(methyl)mercury"),
        ("C[Zn]C", "dimethylzinc"),
        ("C[Pt](C)(C)C", "tetramethylplatinum"),
        ("c1ccccc1[Pd]Br", "bromido(phenyl)palladium"),
        ("[Fe](C#[O+])(C#[O+])(C#[O+])(C#[O+])C#[O+]", "pentacarbonyliron"),
        ("C[Mo](C#[O+])(C#[O+])(C#[O+])", "tricarbonyl(methyl)molybdenum"),
        ("C[Pt](Cl)([NH3])([NH3])([NH3])[NH3]", "tetraamminechlorido(methyl)platinum"),
        ("C[Re]([OH2])Cl", "aquachlorido(methyl)rhenium"),
        ("[Pt](CC)(C)(P(CC)(CC)CC)P(CC)(CC)CC", "ethyl(methyl)bis(triethylphosphane)platinum"),
        ("[Ni](P(c1ccccc1)(c1ccccc1)c1ccccc1)(Cl)Cl", "dichlorido(triphenylphosphane)nickel"),
    ],
)
def test_coordination_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C[Ti](Cl)(Cl)Cl.[Na+]",
        "C[Pt](C)(C)(C)[Pd](C)(C)C",
        "C[Ti](=O)(Cl)Cl",
        "[Fe](C=C)Cl",
        "C[Hg]c1ccc(C(=O)O)cc1",
    ],
)
def test_coordination_out_of_scope_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[BiH2][GeH3]", "germylbismuthane"),
        ("[SbH2][SnH3]", "stannylstibane"),
        ("CC[Pb](CC)(CC)[Bi](C)C", "dimethyl(triethylplumbyl)bismuthane"),
        ("c1ccccc1[Bi](c1ccccc1)[Sn](C)(C)C", "diphenyl(trimethylstannyl)bismuthane"),
        ("C[Sb]([Sn](C)(C)C)[Sn](C)(C)C", "methyldi(trimethylstannyl)stibane"),
        ("CCC[Bi](C)[GeH3]", "germyl(methyl)(propyl)bismuthane"),
        ("C[Sb](C)[Ge](C)(C)[Sn](C)(C)C", "[dimethyl(trimethylstannyl)germyl]di(methyl)stibane"),
    ],
)
def test_metal_pair_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metal_pair_linked_by_carbon_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Bi](c1ccccc1)CCC[Pb](CC)(CC)CC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[Os+]([NH3])([NH3])([NH3])([NH3])[NH3].[Cl-]", "pentaammine(ethyl)osmium(1+) chloride"),
        ("Cl[Fe-](Cl)(Cl)Cl.[Na+]", "sodium tetrachloridoferrate(1-)"),
        ("Cl[Pt-2](Cl)(Cl)(Cl)(Cl)Cl.[K+].[K+]", "dipotassium hexachloridoplatinate(2-)"),
    ],
)
def test_charged_complex_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Mn]([Mn](C#[O+])(C#[O+])(C#[O+])(C#[O+])C#[O+])(C#[O+])(C#[O+])(C#[O+])(C#[O+])C#[O+]", "decacarbonyl-1κ5C,2κ5C-dimanganese(Mn—Mn)"),
        ("[Fe]([Fe](C#[O+])C#[O+])(C#[O+])(C#[O+])C#[O+]", "pentacarbonyl-1κ3C,2κ2C-diiron(Fe—Fe)"),
        ("[Ti](Cl)Cl.C1=C[CH]C=C1.C1=C[CH]C=C1", "dichloridobis(η5-cyclopenta-2,4-dien-1-yl)titanium"),
        ("[Mo+](C#[O+])(C#[O+])(C#[O+])C.C1=C[CH-]C=C1", "tricarbonyl(η5-cyclopenta-2,4-dien-1-yl)(methyl)molybdenum"),
    ],
)
def test_dinuclear_and_hapto_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1[Hg][Sb](c1ccccc1)c1ccccc1", "diphenylstibanyl(phenyl)mercury"),
        ("C[Zn][Ge](C)(C)C", "methyl(trimethylgermyl)zinc"),
    ],
)
def test_class1_metal_with_class2_metal_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
