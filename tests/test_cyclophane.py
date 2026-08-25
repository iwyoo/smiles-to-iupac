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


def test_metacyclophane():
    # [2.2]metacyclophane's PIN, confirmed directly against the Blue Book
    # text (P-26.4.1.4's own paired example alongside the para isomer
    # above, both marked "(PIN)"). PubChem CID 137543's own canonical
    # SMILES for "(2.2)Metacyclophane" is used here and independently
    # matches this project's own SMILES for the same structure, C16H16.
    assert smiles_to_iupac("C1CC2=CC(=CC=C2)CCC3=CC=CC1=C3") == "1,4(1,3)-dibenzenacyclohexaphane"
    assert smiles_to_iupac("c1cc2cc(c1)CCc1cccc(c1)CC2") == "1,4(1,3)-dibenzenacyclohexaphane"


def test_tetrabenzenacyclooctaphane():
    # [1.1.1.1]metacyclophane's PIN, confirmed directly against the Blue
    # Book text (P-26.4.1.4's third worked example, printed immediately
    # alongside the para/meta pair above, all marked "(PIN)"). This
    # project's own independently-built SMILES resolves (via InChIKey) to
    # PubChem CID 11740710, same formula and connectivity, C28H24.
    assert (
        smiles_to_iupac("C1c2cccc(c2)Cc2cccc(c2)Cc2cccc(c2)Cc2cccc1c2")
        == "1,3,5,7(1,3)-tetrabenzenacyclooctaphane"
    )


def test_substituted_tetrabenzenacyclooctaphane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cccc2c1CC1=CC=CC(=C1)CC1=CC=CC(=C1)CC1=CC=CC(=C1)C2")


def test_substituted_paracyclophane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC2=CC=C(CCC3=CC=C1C=C3)C=C2")


def test_substituted_metacyclophane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC2=CC(=CC=C2)CCC3=CC=CC1=C3")


def test_different_bridge_length_raises():
    # [3.3]paracyclophane: same overall shape, three-carbon bridges instead
    # of two-carbon -- a different, unsupported phane, not this exact
    # compound.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2=CC=C(CCCC3=CC=C1C=C3)C=C2")
