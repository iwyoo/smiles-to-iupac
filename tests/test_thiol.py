import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanethiol():
    assert smiles_to_iupac("SC") == "methanethiol"


def test_ethanethiol():
    assert smiles_to_iupac("SCC") == "ethanethiol"


def test_propane_1_thiol():
    assert smiles_to_iupac("SCCC") == "propane-1-thiol"


def test_propane_2_thiol():
    assert smiles_to_iupac("CC(S)C") == "propane-2-thiol"


def test_butane_1_thiol_with_chloro_substituent():
    assert smiles_to_iupac("SCCCCCl") == "4-chlorobutane-1-thiol"


def test_pent_4_ene_1_thiol():
    assert smiles_to_iupac("SCCCC=C") == "pent-4-ene-1-thiol"


def test_dithiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCCS")


def test_cyclic_thiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCCC1")


def test_thiol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCCO")
