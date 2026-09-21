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
    # test_*_radical_name cases above).
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


def test_two_radical_centers_diyl_name():
    assert smiles_to_iupac("[CH2][CH2]") == "ethane-1,2-diyl"


def test_three_radical_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2][CH][CH2]")


def test_halogen_substituted_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C(Cl)")


def test_methylidene_radical_name():
    # Blue Book P-71.2.2.1 worked example: methylidene (mononuclear
    # divalent radical, degree 0, like methyl's own mononuclear case).
    assert smiles_to_iupac("[CH2]") == "methylidene"


def test_methylidyne_radical_name():
    # Blue Book P-71.2.2.1 worked example: methylidyne (mononuclear
    # trivalent radical).
    assert smiles_to_iupac("[CH]") == "methylidyne"


def test_ethylidene_radical_name():
    # Blue Book P-71.2.2.1 worked example: a divalent radical on a chain
    # terminus, one 'idene' appended after the '-yl' name.
    assert smiles_to_iupac("[CH]C") == "ethylidene"


def test_ethylidyne_radical_name():
    # Blue Book P-71.2.2.1 worked example: a trivalent radical on a chain
    # terminus, 'idyne' appended after the '-yl' name.
    assert smiles_to_iupac("[C]C") == "ethylidyne"


def test_propylidene_radical_name():
    assert smiles_to_iupac("[CH]CC") == "propylidene"


def test_propylidyne_radical_name():
    assert smiles_to_iupac("[C]CC") == "propylidyne"


def test_cyclohexylidene_radical_name():
    # Blue Book P-71.2.2.1 worked example: a divalent radical on a
    # monocyclic ring.
    assert smiles_to_iupac("[C]1CCCCC1") == "cyclohexylidene"


def test_cyclobutylidene_radical_name():
    assert smiles_to_iupac("[C]1CCC1") == "cyclobutylidene"


def test_branch_point_divalent_radical_raises():
    # A trivalent radical carbon can never itself be a branch point
    # (3 unpaired electrons leave room for at most one more bond), so only
    # the divalent case needs this check.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[C](C)C")
