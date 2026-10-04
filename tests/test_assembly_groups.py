import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)Cc1ccc(cc1)-c1cc(-c3ccccc3)c(CC(=O)O)cc1",
            "2,2'-([11,21:23,31-terphenyl]-14,24-diyl)diethanoic acid",
        ),
        (
            "OC(=O)Cc1ccc(cc1)-c1cc(CC(=O)O)c(-c3ccccc3)cc1",
            "2,2'-([11,21:24,31-terphenyl]-14,23-diyl)diethanoic acid",
        ),
        (
            "c1(-c3ccc4ccccc4c3CC(=O)O)cccc2c(CC(=O)O)cccc12",
            "2,2'-([1,2'-binaphthalene]-1',5-diyl)diethanoic acid",
        ),
        ("OC(=O)Cc1ccc2ccccc2c1-c1cccc2ccccc12", "2-([1,1'-binaphthalen]-2-yl)ethanoic acid"),
        ("OC(=O)C=C2C=CC(c3ccc4ccccc4c3)c3ccccc23", "2-[[1,2'-binaphthalen]-4(1H)-ylidene]ethanoic acid"),
    ],
)
def test_assembly_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
