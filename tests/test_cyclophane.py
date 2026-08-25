import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_paracyclophane():
    # [2.2]paracyclophane's PIN, confirmed directly against the Blue Book
    # text (P-26.3.2.1 and P-26.4.1.4, both marked "(PIN)"). Two equally
    # valid SMILES (this project's own, and Wikipedia's) both resolve to
    # the same canonical structure, C16H16.
    assert smiles_to_iupac("C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2") == "1,4(1,4)-dibenzenacyclohexaphane"
    assert smiles_to_iupac("C=1C=C2C=CC1CCC3=CC=C(C=C3)CC2") == "1,4(1,4)-dibenzenacyclohexaphane"


def test_substituted_paracyclophane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC2=CC=C(CCC3=CC=C1C=C3)C=C2")


def test_different_bridge_length_raises():
    # [3.3]paracyclophane: same overall shape, three-carbon bridges instead
    # of two-carbon -- a different, unsupported phane, not this exact
    # compound.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2=CC=C(CCCC3=CC=C1C=C3)C=C2")
