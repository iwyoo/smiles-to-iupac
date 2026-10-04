import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("COCCOCCOCCOCCCc1ccc(C(=O)O)cc1", "4-(2,5,8,11-tetraoxatetradecan-14-yl)benzoic acid"),
        ("OC(=O)CC1CCCCCOCCCCC1", "2-(1-oxacyclododecan-7-yl)ethanoic acid"),
        ("OC(=O)CC1CCCCCOCCCCC1C", "2-(6-methyl-1-oxacyclododecan-7-yl)ethanoic acid"),
        ("OC(=O)CC1COCCOCCOCCOCCOCCO1", "2-(1,4,7,10,13,16-hexaoxacyclooctadecan-2-yl)ethanoic acid"),
    ],
)
def test_skeletal_replacement_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
