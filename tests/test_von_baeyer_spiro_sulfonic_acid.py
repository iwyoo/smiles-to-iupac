import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_von_baeyer_sulfonic_acid_name():
    # PubChem-confirmed: CID 14579878, 'bicyclo[2.2.1]heptane-2-sulfonic
    # acid' -- the parent stem's final 'e' is kept ('sulfonic acid' begins
    # with a consonant, P-16.3.3, same as `_thiol.py`'s 'thiol').
    assert smiles_to_iupac("OS(=O)(=O)C1CC2CCC1C2") == "bicyclo[2.2.1]heptane-2-sulfonic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No further real PubChem-registered polycyclic/spiro sulfonic
        # acid structure was found beyond the one PubChem-confirmed case
        # above -- a mechanical single-axis extension of the identical,
        # already-independently-verified mechanism (`_alcohol.py`/
        # `_amine.py`/`_ketone.py`/`_thiol.py`/`_sulfinamide.py`/
        # `_selone.py`/`_sulfonamide.py`/`_sulfinic_acid.py` all share
        # it), reviewed rather than independently PubChem-confirmed (see
        # `implementation-notes.md`'s test-writing policy).
        ("OS(=O)(=O)C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-sulfonic acid"),
        ("OS(=O)(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfonic acid"),
    ],
)
def test_von_baeyer_spiro_sulfonic_acid_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_sulfonic_acid_ring_unsaturation():
    # Ring unsaturation composes with the sulfonic acid suffix locant on
    # the bicyclic/polycyclic branch, mirroring `_ketone.py`'s identical
    # extension -- reviewed rather than independently PubChem-confirmed.
    assert smiles_to_iupac("OS(=O)(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfonic acid"


def test_unsaturated_monospiro_sulfonic_acid_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC2(C1)C=CCCC2")
