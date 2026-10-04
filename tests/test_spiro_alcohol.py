import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # -OH on the ring atom next to the spiro atom, in the smaller
        # ring (P-24.2.1 numbering starts there) -> locant 1, same
        # tie-break `test_spiro.py`'s own methyl case already confirms.
        ("OC1CCCC12CCCCC2", "spiro[4.5]decan-1-ol"),
        # -OH on the ring atom next to the spiro atom, in the larger
        # ring -> locant 6 (matches `test_spiro.py`'s
        # '6-methylspiro[4.5]decane' precedent for the same position).
        # A real registered structure (PubChem CID 12447317,
        # 'C12(CCCC1)C(O)CCCC2') -- note PubChem's own generated name
        # ('spiro[4.5]decan-10-ol') doesn't take the lower of the two
        # available locants for this position the way this project's
        # own pre-existing, already-tested spiro numbering does (see
        # `test_spiro.py`'s identical '6-chlorospiro[4.5]decane' vs.
        # PubChem's own '10-chlorospiro[4.5]decane' for the unmodified
        # module, confirming this is a known PubChem generator quirk,
        # not something new introduced by this module).
        ("C12(CCCC1)C(O)CCCC2", "spiro[4.5]decan-6-ol"),
        # A hydroxyl coexisting with a plain alkyl substituent, both on
        # the ring -- the -OH locant wins the tie-break ahead of the
        # substituent's (mirrors `_von_baeyer_alcohol.py`'s identical
        # suffix-before-substituent rank).
        ("OC1(C)CCCC12CCCCC2", "1-methylspiro[4.5]decan-1-ol"),
        # Equal-size rings: whichever physical ring carries the -OH is
        # numbered first (P-14.4/P-45.2), giving the lowest available
        # locant for that position (2, not the higher locant the other
        # ring-first choice would force) -- same real structure as
        # PubChem CID 90762054 ('OC1CCCC2(C1)CCCCC2'), whose own
        # generated name ('spiro[5.5]undecan-4-ol') has the same known
        # tie-break quirk noted above.
        ("OC1CCCC2(C1)CCCCC2", "spiro[5.5]undecan-2-ol"),
    ],
)
def test_spiro_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_spiro_hydroxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCC12CCCCC2O")


def test_spiro_hydroxyl_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("OCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanol"


def test_unsaturated_spiro_alcohol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCC12C=CCCC2")


def test_stereocenter_alongside_spiro_alcohol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[C@H]1CCCC12CCCCC2")
