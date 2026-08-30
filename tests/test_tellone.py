import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_propane_2_tellone():
    # PubChem structure match: "propane-2-tellone" (same structure as the
    # Blue Book's own "propane-2-thione (PIN)" for the sulfur case).
    assert smiles_to_iupac("CC(=[Te])C") == "propane-2-tellone"


def test_pentane_2_4_ditellone():
    # `_thione.py`'s/`_selone.py`'s own precedent generalizes group-count
    # support without a specific PubChem worked example for every
    # chalcogen (PubChem has no computed IUPACName for this exact
    # structure) -- the shared locant/suffix machinery already handles it.
    assert smiles_to_iupac("CC(=[Te])CC(=[Te])C") == "pentane-2,4-ditellone"


def test_cyclohexanetellone():
    assert smiles_to_iupac("C1CCC(=[Te])CC1") == "cyclohexanetellone"


def test_telluroaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Te]")


def test_polycyclic_tellone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCC2(CCCCC2)CC1")


def test_tellurol_not_confused_with_tellone():
    assert smiles_to_iupac("C[TeH]") == "methanetellurol"


def test_telluride_not_confused_with_tellone():
    assert smiles_to_iupac("C[Te]C") == "methyltellanylmethane"


def test_selone_not_confused_with_tellone():
    assert smiles_to_iupac("CC(=[Se])C") == "propane-2-selone"


def test_thione_not_confused_with_tellone():
    assert smiles_to_iupac("CC(=S)C") == "propane-2-thione"
