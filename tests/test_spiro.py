import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # spiro[n.m]alkane names (P-24.2.1) are real, independently-named
        # compounds; the substituted variants below are covered by
        # in-line rule citations (P-24.2.1 numbering start point,
        # P-14.4/P-45.2 tiebreak) rather than a per-case PubChem lookup.
        ("C1CC12CC2", "spiro[2.2]pentane"),
        ("C1CCC12CCC2", "spiro[3.3]heptane"),
        ("C1CCCC12CCCCC2", "spiro[4.5]decane"),
        # methyl on the ring atom next to the spiro atom, in the smaller ring
        # (P-24.2.1: numbering starts there) -> locant 1.
        ("CC1CCCC12CCCCC2", "1-methylspiro[4.5]decane"),
        # methyl on the ring atom next to the spiro atom, in the larger ring
        # (numbered right after the spiro atom, locant a+1+1 = 6).
        ("C1CCCC12C(C)CCCC2", "6-methylspiro[4.5]decane"),
        # equal-size rings: whichever physical ring carries the substituent
        # is numbered first, so the locant is still 1 (P-14.4/P-45.2).
        ("CC1CCC12CCC2", "1-methylspiro[3.3]heptane"),
        ("C1CC(C)C12CCC2", "1-methylspiro[3.3]heptane"),
    ],
)
def test_smiles_to_iupac_spiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bridged_bicyclic_is_not_spiro():
    # bicyclo[2.2.2]octane: bridged, not spiro; handled by _bicyclic.py, not
    # this module (see test_bicyclic.py for the expected name).
    assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"


def test_fused_bicyclic_is_not_spiro():
    # decahydronaphthalene (decalin): two fused six-membered rings sharing
    # one bond (two atoms), not a single spiro atom; handled by _bicyclic.py,
    # not this module (see test_bicyclic.py for the expected name).
    assert smiles_to_iupac("C1CCC2CCCCC2C1") == "bicyclo[4.4.0]decane"


def test_two_separate_rings_raises():
    # two cyclohexane rings joined by a single bond: two rings, but they
    # share no atom at all, so this is not a spiro system either.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C1CCCCC1")


def test_unsaturated_spiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC12CCCC=C2")


def test_two_ring_heteroatoms_still_raises():
    # two ring heteroatoms (O in each ring) -- out of scope for the
    # single-heteroatom spiro skeletal-replacement module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2(OCCC2)OC1")
