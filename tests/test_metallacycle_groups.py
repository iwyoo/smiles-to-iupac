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
