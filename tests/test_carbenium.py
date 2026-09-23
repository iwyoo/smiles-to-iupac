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


def test_branch_point_carbenium_propan_2_ylium_name():
    # isopropylium (P-73.2.2.1.2's "general method"), mirroring
    # `_radical.py`'s confirmed P-29.3.2.2 "propan-2-yl" mechanism with
    # 'ylium' instead of 'yl'. The cation carbon is itself the branch
    # point (not a chain terminus), unlike the test_*ylium_name cases
    # above.
    assert smiles_to_iupac("C[CH+]C") == "propan-2-ylium"


def test_branch_point_carbenium_butan_2_ylium_name():
    assert smiles_to_iupac("C[CH+]CC") == "butan-2-ylium"


def test_branch_point_carbenium_three_branches_raises():
    # Scoped out for now: no confirmed carbon-only PIN worked example
    # settles this (extra substituent prefix vs. a possible 'tert-butyl'-
    # style retained-name exception), and PubChem doesn't reliably
    # verify these cationic structures either (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C+](C)C")


def test_branch_point_carbenium_with_further_branching_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[CH+]C(C)C")


def test_branched_chain_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C(C)C")


def test_substituted_ring_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC[CH+]C1")


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
