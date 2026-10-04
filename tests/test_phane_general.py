import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("COC1CC2=CC=C(CCC3=CC=C1C=C3)C=C2", "2-methoxy-1,4(1,4)-dibenzenacyclohexaphane"),
        ("COCC1=CC2=CC=C1CCC3=CC=C(CC2)C=C3", "12-(methoxymethyl)-1,4(1,4)-dibenzenacyclohexaphane"),
    ],
)
def test_cyclophane_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclophane_as_substituent():
    name = smiles_to_iupac("OC(=O)CC1=CC2=CC=C1CCC3=CC=C(CC2)C=C3")
    assert "dibenzenacyclohexaphan" in name


def test_biphenyl_assembly_stays_assembly():
    assert smiles_to_iupac("c1ccc(-c2ccc(-c3ccccn3)cc2)cc1") == "2-([1,1'-biphenyl]-4-yl)pyridine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1c2cccc(n2)Cc2cccc(n2)Cc2cccc(n2)Cc2cccc1n2", "1,3,5,7(2,6)-tetrapyridinacyclooctaphane"),
        (
            "OC(=O)C=C1c2cccc(n2)Cc2cccc(n2)Cc2cccc(n2)Cc2cccc1n2",
            "2-[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-2-ylidene]ethanoic acid",
        ),
        (
            "OC(=O)C=C1C=C2Cc3cccc(n3)Cc3cccc(n3)Cc3cccc(n3)CC(=C1)N2",
            "2-[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-14(11H)-ylidene]ethanoic acid",
        ),
    ],
)
def test_pyridine_amplificants_and_ylidene_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
