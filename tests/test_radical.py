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


def test_branch_point_radical_propan_2_yl_name():
    # isopropyl radical, Blue Book P-29.3.2.2 worked example:
    # "propan-2-yl (preferred prefix) (not prop-2-yl)". The radical carbon
    # is itself the branch point (not a chain terminus, unlike the
    # test_*_radical_name cases above) -- this used to raise before
    # tasks/radical-branch-point-naming.md's extension.
    assert smiles_to_iupac("C[CH](C)") == "propan-2-yl"


def test_branch_point_radical_butan_2_yl_name():
    # sec-butyl radical, Blue Book worked example: "butan-2-yl (preferred
    # prefix) (not but-2-yl)" -- 'sec-butyl' itself is general
    # nomenclature only, not the PIN (P-29.6.2.2).
    assert smiles_to_iupac("C[CH]CC") == "butan-2-yl"


def test_branch_point_radical_tert_butyl_name():
    # tert-butyl radical: P-29.6.1's sole retained-name exception -- the
    # unsubstituted (CH3)3C- radical keeps "tert-butyl" as its PIN rather
    # than the general rule's own "2-methylpropan-2-yl".
    assert smiles_to_iupac("[C](C)(C)C") == "tert-butyl"


def test_branch_point_radical_2_methylbutan_2_yl_name():
    # tert-pentyl radical, Blue Book worked example:
    # "2-methylbutan-2-yl (preferred prefix) (not tert-pentyl)" -- unlike
    # tert-butyl, this retained name is NOT a PIN exception.
    assert smiles_to_iupac("CC[C](C)C") == "2-methylbutan-2-yl"


def test_branch_point_radical_with_further_branching_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH](C(C)C)C")


def test_substituted_ring_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC[CH]C1")


def test_multiple_radical_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2][CH2]")


def test_halogen_substituted_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C(Cl)")
