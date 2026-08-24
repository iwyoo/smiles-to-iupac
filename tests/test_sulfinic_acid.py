import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CS(=O)O", "methanesulfinic acid"),
        ("CCS(=O)O", "ethanesulfinic acid"),
        ("CCCS(=O)O", "propane-1-sulfinic acid"),
        ("CC(S(=O)O)C", "propane-2-sulfinic acid"),
        ("CCCCS(=O)O", "butane-1-sulfinic acid"),
    ],
)
def test_saturated_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_sulfinic_acid():
    assert smiles_to_iupac("C=CCS(=O)O") == "prop-2-ene-1-sulfinic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)S(=O)O") == "1-chloroethane-1-sulfinic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C(S(=O)O)C")


def test_two_sulfinic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)CS(=O)O")


def test_sulfonic_acid_not_confused_with_sulfinic():
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_cyclic_sulfinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CCCCC1")
