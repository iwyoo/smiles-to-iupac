import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-tellurinic acid' mirrors '-tellurionic acid'/'-seleninic acid'
        # (see test_telluronic_acid.py/test_seleninic_acid.py) with one
        # fewer oxygen; same locant rules as '-sulfinic acid'. Only the
        # methane case is registered on PubChem (ethane/propane come back
        # as CID 0) -- the same sparse-Te-data gap test_telluronic_acid.py
        # already documented.
        ("C[Te](=O)O", "methanetellurinic acid"),
        ("CC[Te](=O)O", "ethanetellurinic acid"),
        ("CCC[Te](=O)O", "propane-1-tellurinic acid"),
        ("CC([Te](=O)O)C", "propane-2-tellurinic acid"),
        ("CCCC[Te](=O)O", "butane-1-tellurinic acid"),
    ],
)
def test_saturated_tellurinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_tellurinic_acid():
    assert smiles_to_iupac("C=CC[Te](=O)O") == "prop-2-ene-1-tellurinic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)[Te](=O)O") == "1-chloroethane-1-tellurinic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C([Te](=O)O)C")


def test_two_tellurinic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)C[Te](=O)O")


def test_telluronic_acid_not_confused_with_tellurinic_acid():
    assert smiles_to_iupac("C[Te](=O)(=O)O") == "methanetelluronic acid"


def test_ring_tellurinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)C1CCCCC1")


def test_tellurinic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)CCO")
