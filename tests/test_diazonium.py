import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanediazonium():
    # PubChem structure match: "methanediazonium".
    assert smiles_to_iupac("C[N+]#N") == "methanediazonium"


def test_ethanediazonium():
    # PubChem structure match: "ethanediazonium".
    assert smiles_to_iupac("CC[N+]#N") == "ethanediazonium"


def test_propane_1_diazonium():
    # PubChem structure match: "propane-1-diazonium".
    assert smiles_to_iupac("CCC[N+]#N") == "propane-1-diazonium"


def test_propane_2_diazonium():
    # PubChem structure match: "propane-2-diazonium" (a branched chain).
    assert smiles_to_iupac("CC(C)[N+]#N") == "propane-2-diazonium"


def test_chloroethane_diazonium():
    # Locant cited once a substituent is present on a 2-carbon chain, the
    # same project-wide convention as '2-chloroethane-1-selenol'
    # (PubChem's own generated name, "2-chloroethanediazonium", omits it).
    assert smiles_to_iupac("ClCC[N+]#N") == "2-chloroethane-1-diazonium"


def test_pent_4_ene_1_diazonium():
    assert smiles_to_iupac("C=CCCC[N+]#N") == "pent-4-ene-1-diazonium"


def test_two_diazonium_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#[N+]CC[N+]#N")


def test_ring_diazonium_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1[N+]#N")


def test_diazonium_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[N+]#N")


def test_ammonium_not_confused_with_diazonium():
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"
