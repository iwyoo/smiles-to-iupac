import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-sulfinic acid' mirrors '-sulfonic acid' (see
        # test_sulfonic_acid.py) with one fewer oxygen; same locant rules.
        # 'methanesulfinic acid'/'ethanesulfinic acid' are also
        # independently verifiable real compound names.
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


def test_cyclohexanesulfinic_acid():
    # PubChem CID 3302189.
    assert smiles_to_iupac("OS(=O)C1CCCCC1") == "cyclohexanesulfinic acid"


def test_2_methylcyclohexane_1_sulfinic_acid():
    # PubChem CID 67183560.
    assert smiles_to_iupac("OS(=O)C1CCCCC1C") == "2-methylcyclohexane-1-sulfinic acid"


def test_cyclopentanesulfinic_acid():
    # PubChem CID 14138538.
    assert smiles_to_iupac("OS(=O)C1CCCC1") == "cyclopentanesulfinic acid"


def test_2_chlorocyclohexane_1_sulfinic_acid():
    # PubChem CID 67182750.
    assert smiles_to_iupac("OS(=O)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfinic acid"


def test_polycyclic_sulfinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CCCC=C1")


def test_sulfinic_acid_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)CC1CCCCC1")
