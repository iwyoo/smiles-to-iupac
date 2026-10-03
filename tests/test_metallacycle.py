import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-69.4's own worked skeletal-replacement pattern: the metal
        # always sits at locant 1, cited explicitly (unlike the analogous
        # single-heteroatom Hantzsch-Widman case).
        ("[Ni]1CCCC1", "1-nickelacyclopentane"),
        ("[Ti]1CCC1", "1-titanacyclobutane"),
        ("[Fe]1CCCC1", "1-ferracyclopentane"),
        ("[Pt]1CCC1", "1-platinacyclobutane"),
    ],
)
def test_metallacycle_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC1=C(C)[Pt](Cl)(Cl)C(C)=C1C", "1,1-dichloro-2,3,4,5-tetramethyl-1-platinacyclopenta-2,4-diene"),
        ("CC1C[Pt](P(CC)(CC)CC)(P(CC)(CC)CC)CC1", "3-methyl-1,1-bis(triethylphosphane)-1-platinacyclopentane"),
        ("[Ni]1(Cl)CCCC1", "1-chloro-1-nickelacyclopentane"),
        ("[Fe]1(C#[O+])(C#[O+])CCC1", "1,1-dicarbonyl-1-ferracyclobutane"),
    ],
)
def test_metallacycle_with_ligands_and_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metallacycle_with_unsaturated_ligand_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ni]1(C=C)CCCC1")


def test_metallacycle_with_two_metals_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ni]1CC[Ni]C1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "[Si]1(Cl)(Cl)[Fe](C#[O+])(C#[O+])(C#[O+])(C#[O+])CCC1",
            "2,2,2,2-tetracarbonyl-1,1-dichloro-1-sila-2-ferracyclopentane",
        ),
        (
            "CC1=CC(C)=C[Ir](C#[O+])(P(CC)(CC)CC)(P(CC)(CC)CC)=C1",
            "1-carbonyl-3,5-dimethyl-1,1-bis(triethylphosphane)-1-iridabenzene",
        ),
        ("C[Ir]1(C#[O+])=CC(C)=CC(C)=C1", "1-carbonyl-1,3,5-trimethyl-1-iridabenzene"),
    ],
)
def test_hetero_and_metallabenzene_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
