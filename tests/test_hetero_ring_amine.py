import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_aziridine_amine():
    # PubChem CID 20071107 "aziridin-2-amine".
    assert smiles_to_iupac("NC1CN1") == "aziridin-2-amine"


def test_azetidine_amine():
    # PubChem CID 21661087 "azetidin-2-amine".
    assert smiles_to_iupac("NC1CCN1") == "azetidin-2-amine"


def test_pyrrolidine_amine_2():
    # PubChem CID 14298876 "pyrrolidin-2-amine".
    assert smiles_to_iupac("NC1CCCN1") == "pyrrolidin-2-amine"


def test_pyrrolidine_amine_3():
    # PubChem CID 164401 "pyrrolidin-3-amine".
    assert smiles_to_iupac("NC1CCNC1") == "pyrrolidin-3-amine"


def test_piperidine_amine_2():
    # PubChem CID 421842 "piperidin-2-amine".
    assert smiles_to_iupac("NC1CCCCN1") == "piperidin-2-amine"


def test_piperidine_amine_3():
    # PubChem CID 148119 "piperidin-3-amine".
    assert smiles_to_iupac("NC1CCCNC1") == "piperidin-3-amine"


def test_piperidine_amine_4():
    # PubChem CID 424361 "piperidin-4-amine" -- the amine sits directly
    # across the ring from the nitrogen, so both traversal directions
    # give the same locant (4); regression check for a bug where
    # reversing the whole rotated ring list (instead of keeping the ring
    # nitrogen fixed at position 1) moved the nitrogen off locant 1 and
    # produced the wrong locant (3) for this symmetric case.
    assert smiles_to_iupac("NC1CCNCC1") == "piperidin-4-amine"


def test_azepane_amine_2():
    # PubChem CID 12040473 "azepan-2-amine".
    assert smiles_to_iupac("NC1CCCCCN1") == "azepan-2-amine"


def test_azepane_amine_3():
    # PubChem CID 2756445 "azepan-3-amine".
    assert smiles_to_iupac("NC1CCCCNC1") == "azepan-3-amine"


def test_ring_nitrogen_substituent_raises():
    # The ring nitrogen itself carrying a substituent is `_ring_amine.py`'s
    # shape, not this module's -- out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN1CCC(N)CC1")


def test_secondary_exocyclic_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC1CCNCC1")


def test_two_exocyclic_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCNC(N)C1")
