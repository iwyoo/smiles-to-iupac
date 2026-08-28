import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_deuteromethane_name():
    # methane-d1, PubChem CID 12669 (structure-verified via
    # ConnectivitySMILES "[2H]C"; PubChem's own auto-generated name
    # "deuteriomethane" isn't in PIN format, so the name itself is
    # justified by the Blue Book's own worked example, P-82.2.1:
    # "CH3-2H -> (2H1)methane (PIN)").
    assert smiles_to_iupac("[2H]C") == "(2H1)methane"


def test_plain_methane_unaffected():
    assert smiles_to_iupac("C") == "methane"


def test_ethane_deuterium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]CC")


def test_multiple_deuteriums_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]C[2H]")


def test_tritium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[3H]C")
