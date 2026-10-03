import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)C(c1ccccc1)c1ccccc1", "2,2-diphenylethanoic acid"),
        ("OC(=O)CC(c1ccccc1)c1ccccc1", "3,3-diphenylpropanoic acid"),
        ("O=CC(c1ccccc1)c1ccccc1", "2,2-diphenylethanal"),
        ("N#CC(c1ccccc1)c1ccccc1", "2,2-diphenylethanenitrile"),
        ("SC(c1ccccc1)c1ccccc1", "diphenylmethanethiol"),
        ("NC(=O)C(c1ccccc1)c1ccccc1", "2,2-diphenylethanamide"),
        ("OS(=O)(=O)C(c1ccccc1)c1ccccc1", "diphenylmethanesulfonic acid"),
    ],
)
def test_chain_parent_cites_each_aromatic_ring_as_a_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "OC(=O)c1ccc(Cc2ccccc2)cc1",
        "O=Cc1ccc(Cc2ccccc2)cc1",
        "N#Cc1ccc(Cc2ccccc2)cc1",
        "Sc1ccc(Cc2ccccc2)cc1",
    ],
)
def test_functional_group_on_a_ring_is_not_misnamed_as_a_chain_parent(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
