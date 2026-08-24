import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=COC=O", "methanoic anhydride"),
        ("CC(=O)OC(=O)C", "ethanoic anhydride"),
        ("CCC(=O)OC(=O)CC", "propanoic anhydride"),
        ("CCCC(=O)OC(=O)CCC", "butanoic anhydride"),
    ],
)
def test_symmetric_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_symmetric_anhydride():
    assert smiles_to_iupac("CC(C)C(=O)OC(=O)C(C)C") == "2-methylpropanoic anhydride"


def test_halogen_substituent():
    assert smiles_to_iupac("ClCC(=O)OC(=O)CCl") == "2-chloroethanoic anhydride"


def test_unsymmetric_anhydride_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OC(=O)CC")


def test_unsaturated_anhydride_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)OC(=O)C=C")


def test_cyclic_anhydride_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCC(=O)O1")


def test_two_anhydride_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(OC=O)CC(=O)OC=O")
