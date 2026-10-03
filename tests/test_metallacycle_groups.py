import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC2CC[Pt]1(Cl)(Cl)C2c1ccccc1", "1,1-dichloro-7-phenyl-1-platinabicyclo[2.2.1]heptane"),
        ("C1CC2CC[Pt]1(Cl)(Cl)C2N(C)C", "1,1-dichloro-N,N-dimethyl-1-platinabicyclo[2.2.1]heptan-7-amine"),
        ("C1CC2CC[Pt]1(Cl)(Cl)C2O", "1,1-dichloro-1-platinabicyclo[2.2.1]heptan-7-ol"),
        ("C1CC(O)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentan-3-ol"),
        ("C1CC(N)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentan-3-amine"),
        ("C1CC(C(=O)O)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentane-3-carboxylic acid"),
        ("C1CC(=O)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentan-3-one"),
    ],
)
def test_metallacycle_ring_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1Cc2ccccc2C[Pt]1(Cl)Cl", "2,2-dichloro-3,4-dihydro-1H-2-platinanaphthalene"),
        ("C1CCc2ccccc2[Pt]1(Cl)Cl", "1,1-dichloro-3,4-dihydro-2H-1-platinanaphthalene"),
        ("C1Cc2ccccc2[Pt]1(Cl)Cl", "1,1-dichloro-2,3-dihydro-1-platinaindene"),
        ("C1CCC2CCC[Pt](Cl)(Cl)C2C1", "1,1-dichloro-3,4,4a,5,6,7,8,8a-octahydro-2H-1-platinanaphthalene"),
    ],
)
def test_hydro_fused_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=C1CC[Pt]1(Cl)Cl", "1,1-dichloro-2-methylidene-1-platinacyclobutane"),
        ("CC=C1CCC[Pt]1(Cl)Cl", "1,1-dichloro-2-ethylidene-1-platinacyclopentane"),
        ("C1CC[Pt]1([H])Cl", "1-chloro-1-hydrido-1-platinacyclobutane"),
        ("C1CC[Pt]1([H])[H]", "1,1-dihydrido-1-platinacyclobutane"),
        ("C1CC[Pt+]1(Cl)Cl", "1,1-dichloro-1-platinacyclobutan-1-ium"),
        ("C1CC[Pt-]1(Cl)Cl", "1,1-dichloro-1-platinacyclobutan-1-ide"),
        ("C1CC[Pd]1(Cl)Cl", "1,1-dichloro-1-palladacyclobutane"),
        ("C1CC[Zr]1(Cl)Cl", "1,1-dichloro-1-zirconacyclobutane"),
        ("C1CCC[Rh]1(Cl)", "1-chloro-1-rhodacyclopentane"),
    ],
)
def test_metallacycle_ylidene_hydrido_ionic_and_metals(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@H]1CC[Pt]1(Cl)Cl", "(2S)-1,1-dichloro-2-methyl-1-platinacyclobutane"),
        ("C[C@@H]1CC[Pt]1(Cl)Cl", "(2R)-1,1-dichloro-2-methyl-1-platinacyclobutane"),
        ("C[C@H]1C[C@@H](C)C[Pt]1(Cl)Cl", "(2S,4R)-1,1-dichloro-2,4-dimethyl-1-platinacyclopentane"),
    ],
)
def test_metallacycle_ring_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metallacycle_stereo_outside_ring_is_rejected():
    with pytest.raises(Exception):
        smiles_to_iupac("C1CC[Pt]1(Cl)(Cl)[C@H](C)CC")


def test_metallacycle_two_chelates():
    smiles = "C[P]1(C)CC[P](C)(C)[Pt]123([CH2]C[CH2]2)[P](C)(C)CC[P]3(C)C"
    assert smiles_to_iupac(smiles) == "1,1,1,1-bis[ethane-1,2-diylbis(dimethylphosphane)]-1-platinacyclobutane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CN[Pt]1(Cl)Cl", "2,2-dichloro-1-aza-2-platinacyclobutane"),
        ("C1CO[Pt]1(Cl)Cl", "2,2-dichloro-1-oxa-2-platinacyclobutane"),
        ("C1CC[Pt]1=C", "1-methylidene-1-platinacyclobutane"),
    ],
)
def test_heteroatom_rings_and_ylidene_on_ring_metal(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
