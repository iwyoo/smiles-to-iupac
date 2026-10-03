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
    "smiles,expected",
    [
        ("OC(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzoic acid"),
        ("O=Cc1ccc(Cc2ccccc2)cc1", "4-benzylbenzaldehyde"),
        ("N#Cc1ccc(Cc2ccccc2)cc1", "4-benzylbenzonitrile"),
        ("Sc1ccc(Cc2ccccc2)cc1", "4-benzylbenzenethiol"),
        ("NC(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzamide"),
        ("OS(=O)(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzenesulfonic acid"),
        ("NS(=O)(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzenesulfonamide"),
        ("Nc1ccc(Cc2ccccc2)cc1", "4-benzylaniline"),
        ("Oc1ccc(Cc2ccccc2)cc1", "4-benzylphenol"),
        ("Oc1cc(Cl)ccc1Cc1ccccc1", "2-benzyl-5-chlorophenol"),
        ("OC(=O)c1ccc(Cc2ccccn2)cc1", "4-[(pyridin-2-yl)methyl]benzoic acid"),
    ],
)
def test_ring_parent_cites_other_aromatic_rings_as_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_functional_group_on_a_heteroaromatic_ring_beside_another_ring_is_rejected_is_named():
    assert smiles_to_iupac("OC(=O)c1ccncc1Cc1ccccc1") == "3-benzylpyridine-4-carboxylic acid"
