import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-seleninic acid' mirrors '-selenonic acid' (see
        # test_selenonic_acid.py) with one fewer oxygen; same locant rules
        # as '-sulfinic acid'.
        ("C[Se](=O)O", "methaneseleninic acid"),
        ("CC[Se](=O)O", "ethaneseleninic acid"),
        ("CCC[Se](=O)O", "propane-1-seleninic acid"),
        ("CC([Se](=O)O)C", "propane-2-seleninic acid"),
        ("CCCC[Se](=O)O", "butane-1-seleninic acid"),
    ],
)
def test_saturated_seleninic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_seleninic_acid():
    assert smiles_to_iupac("C=CC[Se](=O)O") == "prop-2-ene-1-seleninic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)[Se](=O)O") == "1-chloroethane-1-seleninic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C([Se](=O)O)C")


def test_two_seleninic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)C[Se](=O)O")


def test_sulfinic_acid_not_confused_with_seleninic_acid():
    assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"


def test_selenonic_acid_not_confused_with_seleninic_acid():
    assert smiles_to_iupac("C[Se](=O)(=O)O") == "methaneselenonic acid"


def test_ring_seleninic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)C1CCCCC1")


def test_seleninic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)CCO")
