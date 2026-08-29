import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_guanidine():
    # PubChem structure match: "guanidine".
    assert smiles_to_iupac("NC(=N)N") == "guanidine"


def test_n_substituted_guanidine_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=N)N")


def test_amidine_not_confused_with_guanidine():
    assert smiles_to_iupac("CC(=N)N") == "ethanimidamide"


def test_urea_not_confused_with_guanidine():
    assert smiles_to_iupac("NC(=O)N") == "urea"
