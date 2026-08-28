import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanesulfonic_acid():
    # '-sulfonic acid' suffix construction mirrors '-thiol'/'-ol' locant
    # rules (see test_thiol.py/test_alcohol.py). 'methanesulfonic acid' is
    # also a well-known real compound (a common industrial acid catalyst),
    # independently verifiable.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_ethanesulfonic_acid():
    assert smiles_to_iupac("CCS(=O)(=O)O") == "ethanesulfonic acid"


def test_propane_1_sulfonic_acid():
    assert smiles_to_iupac("CCCS(=O)(=O)O") == "propane-1-sulfonic acid"


def test_propane_2_sulfonic_acid():
    assert smiles_to_iupac("CC(S(=O)(=O)O)C") == "propane-2-sulfonic acid"


def test_chlorobutanesulfonic_acid():
    assert smiles_to_iupac("ClCCCCS(=O)(=O)O") == "4-chlorobutane-1-sulfonic acid"


def test_pent_4_ene_1_sulfonic_acid():
    assert smiles_to_iupac("C=CCCCS(=O)(=O)O") == "pent-4-ene-1-sulfonic acid"


def test_disulfonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CCS(=O)(=O)O")


def test_cyclohexanesulfonic_acid():
    # PubChem CID 428836.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCCC1") == "cyclohexanesulfonic acid"


def test_2_methylcyclohexane_1_sulfonic_acid():
    # PubChem CID 121004858.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCCC1C") == "2-methylcyclohexane-1-sulfonic acid"


def test_cyclopentanesulfonic_acid():
    # PubChem CID 15707015.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCC1") == "cyclopentanesulfonic acid"


def test_2_chlorocyclohexane_1_sulfonic_acid():
    # PubChem CID 129994011.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfonic acid"


def test_polycyclic_sulfonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC=C1")


def test_sulfonic_acid_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CC1CCCCC1")


def test_sulfonic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CCO")
