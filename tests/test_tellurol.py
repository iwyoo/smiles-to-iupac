import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanetellurol():
    # PubChem PUG REST CID 356643, auto-generated name matches exactly.
    assert smiles_to_iupac("C[TeH]") == "methanetellurol"


def test_ethanetellurol():
    # PubChem PUG REST CID 71407255, auto-generated name matches exactly --
    # the same P-14.3.4.2(b) two-carbon locant omission as 'ethaneselenol'/
    # 'ethanethiol'.
    assert smiles_to_iupac("CC[TeH]") == "ethanetellurol"


def test_propane_1_tellurol():
    # PubChem PUG REST CID 71405781, auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[TeH]") == "propane-1-tellurol"


def test_2_methylpropane_1_tellurol():
    # A branched chain. PubChem PUG REST CID 101718199, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CC(C)C[TeH]") == "2-methylpropane-1-tellurol"


def test_prop_2_ene_1_tellurol():
    # An unsaturated chain. PubChem PUG REST CID 101841305, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("C=CC[TeH]") == "prop-2-ene-1-tellurol"


def test_2_chloroethane_1_tellurol():
    # Halogen coexistence (P-35.2.1) -- no PubChem-listed compound found
    # for this specific structure (CID 0), so this is a reviewed result,
    # not an independently verified one: the mechanism itself already has
    # independent confirmation via `_selenol.py`'s/`_thiol.py`'s own
    # identical '2-chloroethane-1-selenol'/'2-chloroethane-1-thiol' cases.
    assert smiles_to_iupac("ClCC[TeH]") == "2-chloroethane-1-tellurol"


def test_ethane_1_2_ditellurol():
    # `_selenol.py`'s/`_thiol.py`'s own precedent generalizes group-count
    # support beyond one without a specific 2+-group PubChem worked
    # example of its own (the shared locant/suffix machinery already
    # handles an arbitrary-length locant list) -- same reasoning applied
    # here.
    assert smiles_to_iupac("[TeH]CC[TeH]") == "ethane-1,2-ditellurol"


def test_propane_1_3_ditellurol():
    assert smiles_to_iupac("[TeH]CCC[TeH]") == "propane-1,3-ditellurol"


def test_cyclohexanetellurol():
    # PubChem structure match: "cyclohexanetellurol".
    assert smiles_to_iupac("C1CCCCC1[TeH]") == "cyclohexanetellurol"


def test_substituted_ring_tellurol_locant_cited():
    assert smiles_to_iupac("CC1CCCCC1[TeH]") == "2-methylcyclohexane-1-tellurol"


def test_polycyclic_tellurol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12(CCC(CC1)CC2)[TeH]")


def test_tellurol_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C[TeH]")


def test_tellurol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[TeH]")
