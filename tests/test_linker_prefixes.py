import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)c1ccc(cc1)[SiH2][SiH2]C[SiH2][SiH2]c1ccc(cc1)C(=O)O",
            "4,4'-[methylenebis(disilane-2,1-diyl)]dibenzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)CCOP(C)OCCc1ccc(cc1)C(=O)O",
            "4,4'-[(methylphosphanediyl)bis(oxyethane-2,1-diyl)]dibenzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)CN(Cl)Cc1ccc(cc1)C(=O)O",
            "4,4'-[(chloroazanediyl)bis(methylene)]dibenzoic acid",
        ),
        ("OC(=O)C=C1C=C2C=CC=CC2=CC1=CC(=O)O", "2,2'-(naphthalene-2,3-diylidene)diethanoic acid"),
        ("OC(=O)C=C1CC(=CC(=O)O)Cc2ccccc12", "2,2'-[naphthalene-1,3(2H,4H)-diylidene]diethanoic acid"),
        ("OC(=O)C=C1CCC(=CC(=O)O)CC1", "2,2'-(cyclohexane-1,4-diylidene)diethanoic acid"),
        ("OC(=O)C=C1C(Cl)=C2C=CC=CC2=CC1=CC(=O)O", "2,2'-(1-chloronaphthalene-2,3-diylidene)diethanoic acid"),
    ],
)
def test_linker_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
