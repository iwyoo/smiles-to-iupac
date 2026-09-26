import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 14626805, 574172, 15091318 respectively).
        ("O=C1CCC2(CC1)CCCC2", "spiro[4.5]decan-8-one"),
        ("O=C1CCC2(CC1)CCCCC2", "spiro[5.5]undecan-3-one"),
        ("O=C1CCC2(CCC2)CC1", "spiro[3.5]nonan-7-one"),
        # Equal-size rings: whichever physical ring carries the ketone is
        # numbered first, giving the lowest available locant (2, not the
        # higher locant the other ring-first choice would force) -- same
        # known-quirk divergence from PubChem's own generated name
        # ('spiro[5.5]undecan-4-one' for this same real structure,
        # PubChem CID 579588) already documented for the identical
        # skeleton position in `test_spiro_alcohol.py`'s
        # 'OC1CCCC2(C1)CCCCC2' -> 'spiro[5.5]undecan-2-ol' case and
        # `test_spiro_amine.py`'s equivalent amine case.
        ("O=C1CCCC2(C1)CCCCC2", "spiro[5.5]undecan-2-one"),
    ],
)
def test_spiro_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_spiro_ketones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC2(C1)CCCC(=O)C2")


def test_unsaturated_spiro_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC12C=CCCC2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No real registered PubChem structure matching "monospiro ketone,
        # single specified stereocenter, no coexisting hydroxyl" turned up
        # in a reasonable search (#1078, M5 step 1) -- both CIP labels
        # below were independently cross-checked with `rdCIPLabeler`
        # directly against atom 4 (the methyl-bearing carbon), not just
        # trusted from this module's own output.
        ("O=C1CC[C@H](C)C12CCCCC2", "(4S)-4-methylspiro[4.5]decan-1-one"),
        ("O=C1CC[C@@H](C)C12CCCCC2", "(4R)-4-methylspiro[4.5]decan-1-one"),
    ],
)
def test_spiro_ketone_specified_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
