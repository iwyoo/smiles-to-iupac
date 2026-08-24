import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC(N)=O", "methyl carbamate"),
        ("CCOC(N)=O", "ethyl carbamate"),
        ("CCCOC(N)=O", "propyl carbamate"),
        ("CCCCOC(N)=O", "butyl carbamate"),
    ],
)
def test_carbamate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)OC(N)=O")


def test_n_substituted_carbamate_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)NC")


def test_free_carbamic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(N)=O")


def test_cyclic_r_group_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)OC1CCCC1")


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCOC(N)=O")
