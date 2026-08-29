import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanetelluronic_acid():
    # PubChem structure match: "methanetelluronic acid" -- the only
    # telluronic acid PubChem has registered at all (ethane/propane/
    # branched candidates all come back as CID 0). The chain-locant
    # mechanism itself is inherited, unverified for tellurium specifically,
    # from the identical mechanism `_selenonic_acid.py`/`_sulfonic_acid.py`
    # already confirm independently.
    assert smiles_to_iupac("C[Te](=O)(=O)O") == "methanetelluronic acid"


def test_ethanetelluronic_acid():
    assert smiles_to_iupac("CC[Te](=O)(=O)O") == "ethanetelluronic acid"


def test_propane_1_telluronic_acid():
    assert smiles_to_iupac("CCC[Te](=O)(=O)O") == "propane-1-telluronic acid"


def test_propane_2_telluronic_acid():
    assert smiles_to_iupac("CC([Te](=O)(=O)O)C") == "propane-2-telluronic acid"


def test_chlorobutanetelluronic_acid():
    assert smiles_to_iupac("ClCCCC[Te](=O)(=O)O") == "4-chlorobutane-1-telluronic acid"


def test_pent_4_ene_1_telluronic_acid():
    assert smiles_to_iupac("C=CCCC[Te](=O)(=O)O") == "pent-4-ene-1-telluronic acid"


def test_ditelluronic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)(=O)C[Te](=O)(=O)O")


def test_selenonic_acid_not_confused_with_telluronic_acid():
    assert smiles_to_iupac("C[Se](=O)(=O)O") == "methaneselenonic acid"


def test_ring_telluronic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)(=O)C1CCCCC1")


def test_telluronic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)(=O)CCO")
