import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_thiocyanate():
    # PubChem structure match: "methyl thiocyanate".
    assert smiles_to_iupac("CSC#N") == "methyl thiocyanate"


def test_ethyl_thiocyanate():
    # PubChem structure match: "ethyl thiocyanate".
    assert smiles_to_iupac("CCSC#N") == "ethyl thiocyanate"


def test_propyl_thiocyanate():
    # PubChem structure match: "propyl thiocyanate".
    assert smiles_to_iupac("CCCSC#N") == "propyl thiocyanate"


def test_branched_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)SC#N")


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CSC#N")


def test_cyclic_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1SC#N")


def test_two_thiocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#CSCCSC#N")


def test_nitrile_not_confused_with_thiocyanate():
    assert smiles_to_iupac("CC#N") == "ethanenitrile"


def test_isothiocyanate_not_confused_with_thiocyanate():
    assert smiles_to_iupac("CN=C=S") == "isothiocyanatomethane"
