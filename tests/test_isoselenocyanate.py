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


def test_phenyl_isoselenocyanate_direct_bond():
    # P-44.1.2.2 rule (1): 'isoselenocyanato' has no suffix form, so the
    # ring is always senior to a chain of the same class -- confirmed by
    # PubChem CID 555335.
    assert smiles_to_iupac("c1ccccc1[N]=C=[Se]") == "isoselenocyanatobenzene"


def test_phenyl_isoselenocyanate_chain():
    # PubChem CID 12506035 gives "isoselenocyanatomethylbenzene" (no
    # parentheses), but this codebase follows the sibling modules'
    # identical, Blue-Book-verified rule instead.
    assert smiles_to_iupac("c1ccccc1C[N]=C=[Se]") == "(isoselenocyanatomethyl)benzene"


def test_phenyl_isoselenocyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[N]=C=[Se]")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N]=C=[Se]")


def test_two_isoselenocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C=NCN=C=[Se]")


def test_isothiocyanate_still_works():
    # Sanity check: the sulfur analogue must not be misrouted here.
    assert smiles_to_iupac("CCN=C=S") == "isothiocyanatoethane"
