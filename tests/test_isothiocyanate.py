import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 10966.
        ("CCN=C=S", "isothiocyanatoethane"),
        # PubChem CID 69403.
        ("CCCN=C=S", "1-isothiocyanatopropane"),
    ],
)
def test_isothiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mononuclear_case_follows_blue_book_prefix_rule():
    # Same reasoning as `_isocyanate.py`'s identical mononuclear case (see
    # module docstring): PubChem's own auto-generated name for this exact
    # structure (CID 11167) is "methylimino(sulfanylidene)methane", a
    # different, non-prefix parent selection -- treated as the same
    # PubChem-autoname-vs-PIN mismatch already settled for the oxygen
    # analogue via P-61.8's 'isocyanatoborane (PIN)' worked example.
    assert smiles_to_iupac("CN=C=S") == "isothiocyanatomethane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(N=C=S)CC1")


def test_phenyl_isothiocyanate_direct_bond():
    # P-44.1.2.2 rule (1): 'isothiocyanato' has no suffix form, so the
    # ring is always senior to a chain of the same class, mirroring
    # `_isocyanate.py`'s identical rule -- confirmed by PubChem CID 7673.
    assert smiles_to_iupac("c1ccccc1N=C=S") == "isothiocyanatobenzene"


def test_phenyl_isothiocyanate_chain():
    # PubChem CID 2346 gives "isothiocyanatomethylbenzene" (no
    # parentheses), but this codebase follows `_isocyanate.py`'s/
    # `_nitro.py`'s identical, Blue-Book-verified rule instead.
    assert smiles_to_iupac("c1ccccc1CN=C=S") == "(isothiocyanatomethyl)benzene"


def test_phenyl_isothiocyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1N=C=S")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCN=C=S")


def test_two_isothiocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C=NCN=C=S")


def test_isocyanate_still_works():
    # Sanity check: the oxygen analogue must not be misrouted here.
    assert smiles_to_iupac("CCN=C=O") == "isocyanatoethane"
