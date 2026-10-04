import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CC2=CC=C(CCC3=CC=C1C=C3)C=C2", "2-hydroxy-1,4(1,4)-dibenzenacyclohexaphane"),
        ("OCC1=CC2=CC=C1CCC3=CC=C(CC2)C=C3", "12-(hydroxymethyl)-1,4(1,4)-dibenzenacyclohexaphane"),
    ],
)
def test_cyclophane_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclophane_as_substituent():
    name = smiles_to_iupac("OC(=O)CC1=CC2=CC=C1CCC3=CC=C(CC2)C=C3")
    assert "dibenzenacyclohexaphan" in name


def test_biphenyl_assembly_stays_assembly():
    assert smiles_to_iupac("c1ccc(-c2ccc(-c3ccccn3)cc2)cc1") == "2-([1,1'-biphenyl]-4-yl)pyridine"
