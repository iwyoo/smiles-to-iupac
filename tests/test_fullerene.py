import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._fullerene import _FULLERENE_C60_SMILES


def test_buckminsterfullerene():
    # C60, cross-checked against PubChem CID 123591's own connectivity
    # SMILES: 60 all-carbon atoms, every atom degree 3, ring perception of
    # exactly 12 five-membered and 20 six-membered rings.
    assert smiles_to_iupac(_FULLERENE_C60_SMILES) == "[60]fullerene"


def test_benzene_still_resolves():
    # a sanity check that the fullerene shape check (60 atoms, all
    # degree-3 carbon) doesn't misfire on an ordinary small ring.
    assert smiles_to_iupac("c1ccccc1") == "benzene"


def test_smaller_cage_raises():
    # a smaller all-carbon cage-like polycyclic (not 60 atoms) must not
    # match; falls through to the ordinary "not supported" path.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C1C2C2C3C4C12")
