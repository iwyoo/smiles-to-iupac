import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Retained names, Blue Book P-72.2.2.2.2's own text.
        ("C[O-]", "methoxide"),
        ("CC[O-]", "ethoxide"),
        ("CCC[O-]", "propoxide"),
        ("CCCC[O-]", "butoxide"),
        # Non-terminal oxygen falls through to the systematic '-olate'
        # suffix, even for a chain length that has a retained name for its
        # terminal position: 'propan-2-olate (PIN)', not 'isopropoxide'
        # (Blue Book's own worked example; PubChem CID 3260420 structure
        # match, and PubChem's own generated name agrees exactly here).
        ("CC(C)[O-]", "propan-2-olate"),
        # Longer unbranched chain, no retained name at that length --
        # PubChem CID 13076479, name matches exactly.
        ("CCCCC[O-]", "pentan-1-olate"),
    ],
)
def test_alkoxide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_alkoxide_with_halogen_substituent():
    # PubChem CID 101079710, "fluoromethanolate" -- a halogen substituent
    # on the retained methoxide's plain shape falls through to the
    # systematic mononuclear-parent path (P-14.3.4.2(a): locant omitted).
    assert smiles_to_iupac("FC[O-]") == "fluoromethanolate"


def test_alkoxide_with_halogen_on_two_carbon_chain():
    # PubChem CID 17767529 (structure match; PubChem's own generated name
    # is "2-fluoroethanolate" -- this project's `_alcohol.py` already
    # established citing the '-1-' locant here for the neutral alcohol
    # ('2-fluoroethan-1-ol', PubChem CID divergence accepted project-wide),
    # so the anion mirrors that same convention with 'ate' appended.
    assert smiles_to_iupac("FCC[O-]") == "2-fluoroethan-1-olate"


def test_alkoxide_with_unsaturation():
    # PubChem CID 21252372, "prop-2-en-1-olate" -- exact match.
    assert smiles_to_iupac("C=CC[O-]") == "prop-2-en-1-olate"


def test_tert_butoxide():
    # Blue Book P-72.2.2.2.2's own text names 'tert-butoxide' as the
    # retained PIN for (CH3)3C-O(-) (structure confirmed via PubChem,
    # which itself returns the systematic '2-methylpropan-2-olate').
    assert smiles_to_iupac("CC(C)(C)[O-]") == "tert-butoxide"


def test_other_branched_alkoxide_raises():
    # A branched skeleton other than tert-butoxide's own fixed shape is
    # still out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C[O-]")


def test_phenoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]c1ccccc1")


def test_two_alkoxide_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]CC[O-]")


def test_ether_oxygen_alongside_alkoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]CCOC")
