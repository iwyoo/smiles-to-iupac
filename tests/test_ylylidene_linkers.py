import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)c1ccc(cc1)C=NOCON=Cc1ccc(cc1)C(=O)O",
            "4,4'-[methylenebis(oxyazanylylidenemethanylylidene)]dibenzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)C=NON=Cc1ccc(cc1)C(=O)O",
            "4,4'-[oxybis(azanylylidenemethanylylidene)]dibenzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)C(Cl)=NCCN=C(Cl)c1ccc(cc1)C(=O)O",
            "4,4'-{ethane-1,2-diylbis[azanylylidene(chloromethanylylidene)]}dibenzoic acid",
        ),
        (
            "OC(=O)CN(CC(=O)O)COCN(CC(=O)O)CC(=O)O",
            "2,2',2'',2'''-[oxybis(methylenenitrilo)]tetraethanoic acid",
        ),
    ],
)
def test_ylylidene_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_imine_linker_without_senior_unit_group_is_not_multiplicative():
    with pytest.raises(Exception):
        smiles_to_iupac("c1ccccc1C=NCCN=Cc1ccccc1")
