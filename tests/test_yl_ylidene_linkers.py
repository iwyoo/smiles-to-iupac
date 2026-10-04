import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCCCC1=CCC2CCCCC2", "1,1'-(ethan-1-yl-2-ylidene)dicyclohexane"),
        ("C1CCCCC1=CCC(C)C2CCCCC2", "1,1'-(butan-3-yl-1-ylidene)dicyclohexane"),
        ("C1CCCCC1=C(CC)C2CCCCC2", "1,1'-(propan-1-yl-1-ylidene)dicyclohexane"),
        ("C1CCCCC1=C(C)C(C)C2CCCCC2", "1,1'-(butan-2-yl-3-ylidene)dicyclohexane"),
        ("C1CCCCC1=CC=C2CCCCC2", "1,1'-(ethane-1,2-diylidene)dicyclohexane"),
        (
            "OC(=O)C1CCC(CC1)=CCC2CCC(C(=O)O)CC2",
            "4,4'-(ethan-1-yl-2-ylidene)di(cyclohexane-1-carboxylic acid)",
        ),
    ],
)
def test_yl_ylidene_linker_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
