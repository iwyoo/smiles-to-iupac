import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 101086150). This project's adamantane
        # already uses the systematic 'tricyclo[3.3.1.1^3,7]decane' name
        # rather than the retained 'adamantane' one (see
        # `test_von_baeyer_thiol.py`'s identical 1-adamantanethiol
        # precedent), so this module's ring_count>=3 path follows the
        # same systematic convention rather than PubChem's own retained-
        # name 'adamantane-2-selenol'.
        ("[SeH]C1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.1^3,7]decane-2-selenol"),
    ],
)
def test_von_baeyer_selenol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No further real PubChem-registered polycyclic/spiro selenol
        # structure was found beyond the one above -- a mechanical
        # single-axis extension of the identical, already-independently-
        # verified mechanism (`_alcohol.py`/`_amine.py`/`_ketone.py`/
        # `_thiol.py`/`_sulfinamide.py`/`_selone.py`/`_sulfonamide.py`/
        # `_sulfinic_acid.py`/`_sulfonic_acid.py`/`_thione.py` all share
        # it), reviewed rather than independently PubChem-confirmed (see
        # `implementation-notes.md`'s test-writing policy).
        ("[SeH]C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2-selenol"),
        ("[SeH]C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-selenol"),
    ],
)
def test_von_baeyer_spiro_selenol_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_selenols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CC2CCC1C([SeH])C2")


def test_von_baeyer_selenol_ring_unsaturation():
    # Ring unsaturation composes with the selenol suffix locant on the
    # bicyclic/polycyclic branch, mirroring `_ketone.py`'s identical
    # extension -- no real PubChem-registered example was found, a
    # mechanical single-axis extension of the already-verified mechanism
    # (reviewed rather than independently PubChem-confirmed).
    assert smiles_to_iupac("[SeH]C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-selenol"


def test_unsaturated_monospiro_selenol_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CCCC2(C1)C=CCCC2")
