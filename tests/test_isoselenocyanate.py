import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 12618729.
        ("CC[N]=C=[Se]", "isoselenocyanatoethane"),
        # PubChem CID 134989633.
        ("CCC[N]=C=[Se]", "1-isoselenocyanatopropane"),
    ],
)
def test_isoselenocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mononuclear_case_follows_blue_book_prefix_rule():
    # Same reasoning as `_isocyanate.py`/`_isothiocyanate.py`'s identical
    # mononuclear case (see module docstring): PubChem's own
    # auto-generated name for this exact structure (CID 138232) is
    # "methylimino(selanylidene)methane", a different, non-prefix parent
    # selection -- treated as the same PubChem-autoname-vs-PIN mismatch.
    assert smiles_to_iupac("C[N]=C=[Se]") == "isoselenocyanatomethane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N]=C=[Se])CC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N]=C=[Se]")


def test_two_isoselenocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C=NCN=C=[Se]")


def test_isothiocyanate_still_works():
    # Sanity check: the sulfur analogue must not be misrouted here.
    assert smiles_to_iupac("CCN=C=S") == "isothiocyanatoethane"
