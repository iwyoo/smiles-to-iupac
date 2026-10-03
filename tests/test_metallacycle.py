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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C1CC2C(C1)C[Ti]2.C1=C[CH]C=C1.C1=C[CH]C=C1",
            "6,6-di(η5-cyclopenta-2,4-dien-1-yl)-6-titanabicyclo[3.2.0]heptane",
        ),
        (
            "[Pt]1(P(CC)(CC)CC)(P(CC)(CC)CC)C2CC1C(OC)CCC2OC",
            "2,5-dimethoxy-7,7-bis(triethylphosphane)-7-platinabicyclo[4.1.1]octane",
        ),
        ("C1CC2C(C1)C[Ti]2(C)C", "6,6-dimethyl-6-titanabicyclo[3.2.0]heptane"),
    ],
)
def test_bicyclic_metallacycle_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CC2C(C1)C[Ti]2(C)C", "6,6-dimethyl-6-titanabicyclo[3.2.0]hept-3-ene"),
        ("C1=CC2C=CC1[Pt]2", "7-platinabicyclo[2.2.1]hepta-2,5-diene"),
    ],
)
def test_unsaturated_bicyclic_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1c2ccccc2[Pt](P(C)(C)C)(P(C)(C)C)c2ccccc12", "9,9-bis(trimethylphosphane)-10H-9-platinaanthracene"),
        ("C1c2cccc(C)c2[Pt](Cl)(Cl)c2ccccc12", "9,9-dichloro-1-methyl-10H-9-platinaanthracene"),
    ],
)
def test_metallaanthracene_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C1c2ccccc2[Pt]2(P(C)(C)CP2(C)C)c2ccccc12",
            "9,9-[methylenebis(dimethylphosphane)]-10H-9-platinaanthracene",
        ),
        (
            "[Pt]1(P(c2ccccc2)(c2ccccc2)c2ccccc2)(P(c2ccccc2)(c2ccccc2)c2ccccc2)C2CC1C(OC)CCC2OC",
            "2,5-dimethoxy-7,7-bis(triphenylphosphane)-7-platinabicyclo[4.1.1]octane",
        ),
        (
            "CC1=C(C)[Pt](P(c2ccccc2)(c2ccccc2)c2ccccc2)(P(c2ccccc2)(c2ccccc2)c2ccccc2)C(C)=C1C",
            "2,3,4,5-tetramethyl-1,1-bis(triphenylphosphane)-1-platinacyclopenta-2,4-diene",
        ),
    ],
)
def test_metallacycles_with_ring_and_chelating_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1c2ccccc2[Pt](Cl)(Cl)C=C1", "1,1-dichloro-4H-1-platinanaphthalene"),
        ("C1=Cc2ccccc2C[Pt]1(Cl)Cl", "2,2-dichloro-1H-2-platinanaphthalene"),
        ("C1c2cccc(C)c2[Pt](Cl)(Cl)C=C1", "1,1-dichloro-8-methyl-4H-1-platinanaphthalene"),
        ("C1c2ccccc2-c2ccccc2[Pt]1(Cl)Cl", "9,9-dichloro-10H-9-platinaphenanthrene"),
        ("C1=Cc2ccccc2[Pt]1(Cl)Cl", "1,1-dichloro-1-platinaindene"),
        ("C1c2ccccc2-c2ccccc2[Pt]1(Cl)Cl".replace("C1c2ccccc2-c2ccccc2", "C1c2ccccc2-c2ccccc2"), "9,9-dichloro-10H-9-platinaphenanthrene"),
    ],
)
def test_fused_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Ti]1(CCC1)2CCC2", "4-titanaspiro[3.3]heptane"),
        ("C1CC2(C1)C[Ti]C2", "2-titanaspiro[3.3]heptane"),
        ("C1CCC2(CC1)C[Pt]C2", "2-platinaspiro[3.5]nonane"),
        ("C1C2CC3C1C[Ti]3C2", "3-titanatricyclo[3.2.1.0^3,6]octane"),
        ("C1C2C(C)CC1C[Ti]2(Cl)Cl", "2,2-dichloro-6-methyl-2-titanabicyclo[2.2.1]heptane"),
    ],
)
def test_spiro_and_polycyclic_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
