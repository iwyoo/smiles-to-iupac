import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_epoxycholestane():
    # PubChem CID 281912; CAS-style synonym "5,6-epoxycholestane" (also
    # listed as the fully stereo-specified "5.alpha.,6.alpha.-Epoxycholestane").
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC4C5(C3(CCCC5)C)O4)C")
        == "5,6-epoxycholestane"
    )


def test_plain_cholestane_unaffected():
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"


def test_transannular_bridge_not_claimed():
    # A one-atom -O- bridge whose two neighbors are NOT already directly
    # bonded (a genuine 1,4-type transannular span, not this module's
    # ortho-fused-epoxide shape) must not be claimed here -- confirmed by
    # the whole molecule falling through to `UnsupportedStructure` rather
    # than any module (this one included) producing a wrong "epoxy" name.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC12CCCC1C1CCC3CC4CCCCC4(O3)C1CC2")
