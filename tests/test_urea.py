import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_urea():
    # PubChem structure match: "urea".
    assert smiles_to_iupac("NC(=O)N") == "urea"


def test_n_substituted_urea_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=O)N")


def test_n_n_disubstituted_urea_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=O)NC")


def test_thiourea_not_confused_with_urea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=S)N")
