import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 21483273, 15719642, 53485543
        # respectively).
        ("NC1CCC2(CC1)CCCC2", "spiro[4.5]decan-8-amine"),
        ("NC1CCC2(CC1)CCCCC2", "spiro[5.5]undecan-3-amine"),
        ("NC1CCC2(CCC2)CC1", "spiro[3.5]nonan-7-amine"),
        # Equal-size rings: whichever physical ring carries the amine is
        # numbered first, giving the lowest available locant (2, not the
        # higher locant the other ring-first choice would force) -- same
        # known-quirk divergence from PubChem's own generated name
        # ('spiro[5.5]undecan-4-amine' for this same real structure,
        # PubChem CID 21483267) already documented for the identical
        # skeleton position in `test_spiro_alcohol.py`'s
        # 'OC1CCCC2(C1)CCCCC2' -> 'spiro[5.5]undecan-2-ol' case.
        ("NC1CCCC2(C1)CCCCC2", "spiro[5.5]undecan-2-amine"),
    ],
)
def test_spiro_amine_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_spiro_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC2(C1)CCCCC2N")


def test_spiro_amine_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC1CCCC12CCCCC2")


def test_unsaturated_spiro_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC12C=CCCC2")
