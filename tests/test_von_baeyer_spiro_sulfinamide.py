import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_von_baeyer_sulfinamide_name():
    # PubChem-confirmed: CID 165459501, 'bicyclo[2.2.1]heptane-2-
    # sulfinamide' -- the parent stem's final 'e' is kept ('sulfinamide'
    # begins with a consonant, P-16.3.3, same as `_thiol.py`'s 'thiol').
    assert smiles_to_iupac("NS(=O)C1CC2CCC1C2") == "bicyclo[2.2.1]heptane-2-sulfinamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No further real PubChem-registered polycyclic/spiro sulfinamide
        # structure was found beyond the one PubChem-confirmed case above
        # (sulfinamides are rare in the registry at this ring-system
        # scale) -- these are a mechanical single-axis extension of the
        # identical, already-independently-verified mechanism
        # `_alcohol.py`/`_amine.py`/`_ketone.py`/`_thiol.py` all share
        # (only the suffix word itself varies), reviewed rather than
        # independently PubChem-confirmed (see `implementation-notes.md`'s
        # test-writing policy).
        ("NS(=O)C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-sulfinamide"),
        ("NS(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfinamide"),
    ],
)
def test_von_baeyer_spiro_sulfinamide_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_sulfinamide_on_polycyclic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNS(=O)C1CC2CCC1C2")


def test_multiple_ring_sulfinamides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)C1CC2CCC1C2S(N)=O")


def test_von_baeyer_sulfinamide_ring_unsaturation():
    # Ring unsaturation composes with the sulfinamide suffix locant on
    # the bicyclic/polycyclic branch, mirroring `_ketone.py`'s identical
    # extension -- reviewed rather than independently PubChem-confirmed.
    assert smiles_to_iupac("NS(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfinamide"


def test_unsaturated_monospiro_sulfinamide_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)C1CCCC2(C1)C=CCCC2")
