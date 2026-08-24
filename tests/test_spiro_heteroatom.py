import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # spiro[4.5]decane with the heteroatom adjacent to the spiro atom
        # in the smaller (5-membered) ring -- cross-checked against
        # PubChem CID 79743 (8-oxaspiro[4.5]decane-7,9-dione), which
        # numbers its own ring-1 oxa position the same way.
        ("C1CCC2(CCCCC2)O1", "1-oxaspiro[4.5]decane"),
        # same skeleton, heteroatom adjacent to the spiro atom in the
        # larger (6-membered) ring instead -- ring-size numbering priority
        # is unaffected by the heteroatom (P-24.2.1: "never modified by
        # the introduction of heteroatoms"), so this is locant 6, not 10;
        # cross-checked against PubChem CID 12630235
        # (6-Oxaspiro[4,5]decane).
        ("C1CCC2(OCCCC2)C1", "6-oxaspiro[4.5]decane"),
        ("C1CCC2(NCCCC2)C1", "6-azaspiro[4.5]decane"),
        ("C1CCC2(SCCCC2)C1", "6-thiaspiro[4.5]decane"),
        # equal-size rings: the heteroatom's locant breaks the
        # otherwise-free choice of which ring is numbered first, ahead of
        # any substituent (P-14.4/P-45.2's usual tie-break, superseded
        # here the same way it is in _von_baeyer_heteroatom.py) -- the
        # O-bearing ring is numbered first, giving O locant 1 rather than
        # numbering the all-carbon ring first (which would put O at 6).
        ("C1CCC2(OCCC2)C1", "1-oxaspiro[4.4]nonane"),
    ],
)
def test_smiles_to_iupac_spiro_heteroatom(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_heteroatoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2(OCCC2)OC1")


def test_heteroatom_outside_ring_raises():
    # a plain hydrocarbon spiro skeleton with an exocyclic -OH
    # substituent: the heteroatom isn't a *ring* atom, so this is a
    # different module's territory (deferred, see _alcohol.py).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2(CCCCC2)C1O")


def test_unsupported_heteroatom_element_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2(PCCC2)C1")


def test_unsaturated_heteroatom_spiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2(OC=CCC2)C1")
