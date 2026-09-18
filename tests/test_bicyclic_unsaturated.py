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


def test_bicyclic_compound_locant_unsaturation_raises():
    # The double bond sits between a bridgehead and the lone atom of the
    # uniquely shortest (1-atom) bridge, which P-23.2.3's own numbering
    # rule always places last (never adjacent to locant 1 under any
    # candidate numbering) -- a genuine compound-locant case, deferred to
    # a separate follow-up milestone.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=C2CCC1CC2")


def test_bicyclic_unsaturated_saturated_case_unaffected():
    assert smiles_to_iupac("C1CCC2CC1CC2") == "bicyclo[3.2.1]octane"
