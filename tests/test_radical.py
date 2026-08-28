import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_radical_name():
    # Blue Book P-71.2.1.1 worked example: "*CH3 -> methyl (PIN)". PubChem
    # cannot structurally cross-check radical SMILES (it silently strips
    # the radical electron and normalizes to the closed-shell parent
    # alkane, e.g. [CH3] -> CID 297 "methane") -- confirmed by this
    # project's own scoping pass, so this and the other radical tests
    # below rely solely on the Blue Book's own cited worked examples.
    assert smiles_to_iupac("[CH3]") == "methyl"


def test_ethyl_radical_name():
    assert smiles_to_iupac("[CH2]C") == "ethyl"


def test_propyl_radical_name():
    # Blue Book P-71.2.1.1 worked example: a terminal radical on propane
    # -> "propyl (PIN)".
    assert smiles_to_iupac("[CH2]CC") == "propyl"


def test_pentyl_radical_name():
    assert smiles_to_iupac("[CH2]CCCC") == "pentyl"


def test_cyclobutyl_radical_name():
    # Blue Book P-71.2.1.1 worked example: a cyclobutane ring radical ->
    # "cyclobutyl (PIN)".
    assert smiles_to_iupac("C1C[CH]C1") == "cyclobutyl"


def test_cyclopentyl_radical_name():
    assert smiles_to_iupac("C1CC[CH]C1") == "cyclopentyl"


def test_branched_chain_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C(C)C")


def test_nonterminal_chain_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[CH](C)")


def test_substituted_ring_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC[CH]C1")


def test_multiple_radical_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2][CH2]")


def test_halogen_substituted_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C(Cl)")
