import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methylium_name():
    # Blue Book P-73.2.2.1.1 worked example: "[CH3]+ -> methylium (PIN)".
    assert smiles_to_iupac("[CH3+]") == "methylium"


def test_ethylium_name():
    assert smiles_to_iupac("[CH2+]C") == "ethylium"


def test_propylium_name():
    # Blue Book P-73.2.2.1.1 worked example: a terminal cation on propane
    # -> "propylium (PIN)".
    assert smiles_to_iupac("CC[CH2+]") == "propylium"


def test_pentylium_name():
    assert smiles_to_iupac("CCCC[CH2+]") == "pentylium"


def test_cyclobutylium_name():
    # Blue Book P-73.2.2.1.1 worked example: a cyclobutane ring cation ->
    # "cyclobutylium (PIN)".
    assert smiles_to_iupac("C1C[CH+]C1") == "cyclobutylium"


def test_cyclopentylium_name():
    assert smiles_to_iupac("C1CC[CH+]C1") == "cyclopentylium"


def test_branch_point_carbenium_raises():
    # A cation carbon that is itself a branch point (isopropylium) needs
    # P-73.2.2.1.2's "general method" (out of scope for this module).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[CH+]C")


def test_branched_chain_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C(C)C")


def test_substituted_ring_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC[CH+]C1")


def test_polycyclic_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[C+]12CCC1CC2")


def test_multiple_carbenium_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+][CH2+]")


def test_halogen_substituted_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C(Cl)")


def test_unsaturated_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C=C")


def test_radical_not_confused_with_carbenium():
    assert smiles_to_iupac("[CH2]CC") == "propyl"


def test_ammonium_not_confused_with_carbenium():
    assert smiles_to_iupac("[NH4+]") == "azanium"
