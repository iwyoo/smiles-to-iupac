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
        "C[Pt](C)(C)(C)[Pt](C)(C)C",
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
