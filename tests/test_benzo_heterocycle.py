import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2[nH]cnc2c1", "1H-benzimidazole"),
        ("c1ccc2ocnc2c1", "1,3-benzoxazole"),
        ("c1ccc2scnc2c1", "1,3-benzothiazole"),
        ("c1ccc2c(c1)OCO2", "2H-1,3-benzodioxole"),
        ("C1=COc2ccccc2O1", "1,4-benzodioxine"),
        ("c1ccc2c(c1)OCCO2", "2,3-dihydro-1,4-benzodioxine"),
        ("c1ccc2[nH]ncc2c1", "1H-indazole"),
        ("c1ccc2n[nH]cc2c1", "2H-indazole"),
        ("c1ccc2[nH]nnc2c1", "1H-benzotriazole"),
        ("c1ccc2c[nH]cc2c1", "2H-isoindole"),
        ("C1=Cc2ccccc2OC1", "2H-1-benzopyran"),
        ("O1C=CC=Cc2ccccc12", "1-benzoxepine"),
        ("Cc1nc2ccccc2[nH]1", "2-methyl-1H-benzimidazole"),
        ("Nc1ccc2[nH]cnc2c1", "1H-benzimidazol-5-amine"),
        ("Cc1ccc2c(c1)OCO2", "5-methyl-2H-1,3-benzodioxole"),
        ("O=C(O)c1ccc2OCOc2c1", "2H-1,3-benzodioxole-5-carboxylic acid"),
        ("O=c1[nH]c2ccccc2o1", "1,3-benzoxazol-2(3H)-one"),
        ("O=C1NC(=O)c2ccccc12", "1H-isoindole-1,3(2H)-dione"),
    ],
)
def test_benzo_heterocycle_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc(-c2nc3ccccc3o2)cc1", "2-phenyl-1,3-benzoxazole"),
        ("c1ccc(-c2ccc3ccccc3n2)cc1", "2-phenylquinoline"),
        ("c1ccc(Cc2ccc3ccccc3c2)cc1", "2-benzylnaphthalene"),
        ("Clc1ccc2[nH]c(-c3ccccc3)cc2c1", "5-chloro-2-phenyl-1H-indole"),
        ("CN1CCc2ccccc2C1c1ccccc1", "2-methyl-1-phenyl-1,2,3,4-tetrahydroisoquinoline"),
    ],
)
def test_senior_fused_system_is_the_parent_over_other_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
