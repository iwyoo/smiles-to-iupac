import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCCc4ccc(c7ccc8ccccc8c7)cc4", "2-[4-(naphthalen-2-yl)phenyl]ethanol"),
        ("OC(=O)Cc4ccc(c7ccc8ccccc8c7)cc4", "2-[4-(naphthalen-2-yl)phenyl]ethanoic acid"),
        ("OCCc4ccc(c7ccc8ccccc8c7)nc4", "2-[6-(naphthalen-2-yl)pyridin-3-yl]ethanol"),
        ("OCCc4ccc(C7CCCCO7)nc4", "2-[6-(oxan-2-yl)pyridin-3-yl]ethanol"),
        ("OCCOC(C)(C)C(Cl)Cl", "2-[(1,1-dichloro-2-methylpropan-2-yl)oxy]ethanol"),
        ("OCCOC1CCCCC1Cl", "2-[(2-chlorocyclohexyl)oxy]ethanol"),
        ("OCCOc1ccc(Cl)nc1", "2-[(6-chloropyridin-3-yl)oxy]ethanol"),
    ],
)
def test_complex_substituted_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
