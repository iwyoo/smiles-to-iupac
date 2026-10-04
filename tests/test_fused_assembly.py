import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1cc2ccccc2cc1-c1ccc2ccccc2c1-c1ccc2ccccc2c1", "12,21:22,32-ternaphthalene"),
        ("c1ccc2cc(ccc2c1)-c1ccc2cc(ccc2c1)-c1ccc2ccccc2c1", "12,22:26,32-ternaphthalene"),
        ("Oc1ccc2cc(ccc2c1)-c1ccc2cc(ccc2c1)-c1ccc2ccccc2c1", "[12,22:26,32-ternaphthalen]-16-ol"),
        ("Brc1ccc2cc(ccc2c1)-c1ccc2cc(ccc2c1)-c1ccc2ccccc2c1", "16-bromo-12,22:26,32-ternaphthalene"),
        ("OC(=O)c1cc2ccccc2cc1-c1ccc2ccccc2c1-c1ccc2ccccc2c1", "[12,21:22,32-ternaphthalene]-33-carboxylic acid"),
        ("c1ccc2[nH]c(cc2c1)-c1cc2ccccc2n1-c1cc2ccccc2[nH]1", "11H,31H-12,21:22,32-terindole"),
        ("c1ccc2ncc(cc2c1)-c1cnc2ccccc2c1-c1cnc2ccccc2c1", "13,23:24,33-terquinoline"),
    ],
)
def test_assemblies_of_three_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccc2c(c1)ccn2-c1cc2ccccc2[nH]1", "1'H-1,2'-biindole"),
        ("C1=Cc2ccccc2C1-C1C=Cc2ccccc21", "1H,1'H-1,1'-biindene"),
        ("C1Cc2ccccc2C1C1Cc2ccccc2C1", "2,2',3,3'-tetrahydro-1H,1'H-1,2'-biindene"),
        ("C1CCc2ccccc2C1C1CCc2ccccc2C1", "1,1',2,2',3,3',4,4'-octahydro-1,2'-binaphthalene"),
        ("OC1Cc2ccccc2C1C1Cc2ccccc2C1", "2,2',3,3'-tetrahydro-1H,1'H-[1,2'-biinden]-2-ol"),
        ("C1CCC2CCCCC2C1-C1CCC2CCCCC2C1", "icosahydro-1,2'-binaphthalene"),
        ("C1Cc2ccccc2C1C1Cc2ccccc2C1C1Cc2ccccc2C1", "12,13,22,23,32,33-hexahydro-11H,21H,31H-11,22:21,32-terindene"),
    ],
)
def test_hydro_and_indicated_hydrogen_of_fused_assembly_components(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1(=C2c3ccccc3-c3ccccc23)c2ccccc2-c2ccccc12", "9,9'-bifluorenylidene"),
        ("C1(=C2c3ccccc3-c3ccccc23)c2ccccc2-c2cc(O)ccc12", "[9,9'-bifluorenylidene]-3-ol"),
        ("O=C1c2ccccc2C(=C2c3ccccc3C(=O)c3ccccc23)c2ccccc12", "[9,9'-bianthracenylidene]-10,10'-dione"),
    ],
)
def test_double_bond_junction_between_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
