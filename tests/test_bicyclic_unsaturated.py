import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-31.1.4.1's own worked example: 'bicyclo[3.2.1]oct-2-ene
        # (PIN)'.
        ("C1=CC2CCC(C1)C2", "bicyclo[3.2.1]oct-2-ene"),
        # A double bond spanning locants 2-3 of a bicyclo[2.2.2]octane skeleton.
        ("C1CC2CCC1C=C2", "bicyclo[2.2.2]oct-2-ene"),
        # A bridgehead alkene (double bond touches locant 1 itself).
        ("C1CC2CC=C1C2", "bicyclo[2.2.1]hept-1-ene"),
        # A halogen substituent coexisting with the ring double bond -- the
        # bond locant still wins the numbering tie-break the same way
        # `_von_baeyer_heteroatom.py`'s heteroatom locant does.
        ("ClC1CC2CCC1C=C2", "5-chlorobicyclo[2.2.2]oct-2-ene"),
        # A triple bond on the same bicyclo[3.2.1]octane skeleton as the
        # pilot example above.
        ("C1CC2CC#CC1C2", "bicyclo[3.2.1]oct-2-yne"),
        # Two double bonds (a bicyclic diene), multiplying prefix 'di' plus
        # the euphonic stem 'a' insertion (P-31.1.1.2).
        ("C1CC2C=CCC=C1C2", "bicyclo[4.2.1]nona-1,4-diene"),
    ],
)
def test_bicyclic_unsaturated_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # The double bond sits between a bridgehead and the lone atom of
        # the uniquely shortest (1-atom) bridge, which P-23.2.3's own
        # numbering rule always places last (never adjacent to locant 1
        # under any candidate numbering) -- a genuine compound-locant case
        # (P-31.1.4.2(1), the higher locant cited in parentheses).
        # PubChem-confirmed, CID 53949421.
        ("C1=C2CCC1CC2", "bicyclo[2.2.1]hept-1(7)-ene"),
        # Same shape on larger bridge-length combinations, each PubChem-
        # confirmed (CID 85836581, 85912210, 17876085, 154531328,
        # 149179167 respectively).
        ("C1=C2CCCC1CC2", "bicyclo[3.2.1]oct-1(8)-ene"),
        ("C1=C2CCCCC1CC2", "bicyclo[4.2.1]non-1(9)-ene"),
        ("C1=C2CCCC1CCC2", "bicyclo[3.3.1]non-1(9)-ene"),
        ("C1=C2CCCCC1CCC2", "bicyclo[4.3.1]dec-1(10)-ene"),
        ("C1=C2CCCCCC1CC2", "bicyclo[5.2.1]dec-1(10)-ene"),
        # The Blue Book's own P-31.1.4.2 worked example: 'bicyclo[4.2.0]
        # oct-6-ene (PIN) [not bicyclo[4.2.0]oct-1(8)-ene]' -- confirms
        # criterion (1) (fewest compound locants) correctly prefers this
        # molecule's plain '-6-' numbering over the compound-locant
        # alternative. Notably, PubChem's own auto-generated name for
        # this exact structure (CID 12729262) is the non-PIN
        # 'bicyclo[4.2.0]oct-1(8)-ene' the Blue Book explicitly rejects --
        # this project follows the Blue Book's own worked example
        # instead, the same kind of documented PubChem-generator
        # divergence already noted elsewhere in this project's tests.
        ("C1=C2CCCCC2C1", "bicyclo[4.2.0]oct-6-ene"),
        # A diene mixing one compound locant and one plain locant in the
        # same name (P-31.1.4.2's citation joins both kinds the same way
        # `_unsaturation_suffix_from_citations` formats any diene).
        # PubChem-confirmed, CID 154530450.
        ("C1=CC2C=C(C1)CC2", "bicyclo[3.2.1]octa-1(8),3-diene"),
    ],
)
def test_bicyclic_compound_locant_unsaturation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bicyclic_unsaturated_saturated_case_unaffected():
    assert smiles_to_iupac("C1CCC2CC1CC2") == "bicyclo[3.2.1]octane"
