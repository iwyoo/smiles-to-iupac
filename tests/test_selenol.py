import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methaneselenol():
    # PubChem PUG REST CID 440764, auto-generated name matches exactly.
    assert smiles_to_iupac("C[SeH]") == "methaneselenol"


def test_ethaneselenol():
    # PubChem PUG REST CID 5252527, auto-generated name matches exactly --
    # the same P-14.3.4.2(b) two-carbon locant omission as 'ethanethiol'.
    assert smiles_to_iupac("CC[SeH]") == "ethaneselenol"


def test_propane_1_selenol():
    # PubChem PUG REST CID 71373846, auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[SeH]") == "propane-1-selenol"


def test_2_methylpropane_1_selenol():
    # A branched chain. PubChem PUG REST CID 157368643, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CC(C)C[SeH]") == "2-methylpropane-1-selenol"


def test_prop_2_ene_1_selenol():
    # An unsaturated chain. PubChem PUG REST CID 15821407, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("C=CC[SeH]") == "prop-2-ene-1-selenol"


def test_2_chloroethane_1_selenol():
    # Halogen coexistence (P-35.2.1) -- no PubChem-listed compound found
    # for this specific structure (every halogenated selenol SMILES tried
    # came back as CID 0), so this is a reviewed result, not an
    # independently verified one: the mechanism itself already has
    # independent confirmation via `_thiol.py`'s own identical
    # '2-chloroethane-1-thiol' case (same locant-citation rule once a
    # substituent is present on the two-carbon chain).
    assert smiles_to_iupac("ClCC[SeH]") == "2-chloroethane-1-selenol"


def test_diselenol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]CC[SeH]")


def test_selenide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se]C")


def test_ring_selenol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1[SeH]")


def test_selenol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[SeH]")
