import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCC[Se]C", "2-(methylselanyl)ethanol"),
        ("NCC[Te]CC", "2-(ethyltellanyl)ethanamine"),
        ("OC(=O)C[Se]C", "2-(methylselanyl)ethanoic acid"),
        ("OC(=O)c1ccc(cc1)[SeH]", "4-selanylbenzoic acid"),
        ("OC(=O)c1ccc(cc1)[Se]C", "4-(methylselanyl)benzoic acid"),
        ("OC(=O)C1CCC(=[Se])CC1", "4-selanylidenecyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(=[Te])CC1", "4-tellanylidenecyclohexane-1-carboxylic acid"),
        ("CC(=O)C[Se]C", "1-(methylselanyl)propan-2-one"),
        ("OCC[Se]C(C)C", "2-[(propan-2-yl)selanyl]ethanol"),
    ],
)
def test_selanyl_and_tellanyl_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
