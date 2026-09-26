import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No real PubChem-registered polycyclic/spiro tellone structure
        # was found (sparse tellone coverage in general, same caveat
        # `_tellone.py`'s own module docstring and existing monocyclic
        # tests already document) -- verified instead by structural
        # analogy against the identical, already-independently-verified
        # mechanism on the same ring shapes for the sulfur/selenium
        # analogues (`test_von_baeyer_spiro_thione.py`'s
        # `S=C1CC2CCC1C2`/CID 12859496, `test_von_baeyer_spiro_selone.py`'s
        # matching bicyclic case), reviewed rather than independently
        # PubChem-confirmed (see `implementation-notes.md`'s
        # test-writing policy).
        ("[Te]=C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2-tellone"),
        ("[Te]=C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-tellone"),
        ("[Te]=C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-tellone"),
        ("[Te]=C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-tellone"),
    ],
)
def test_von_baeyer_spiro_tellone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_tellones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CC2CCC1C(=[Te])C2")


def test_von_baeyer_tellone_ring_unsaturation():
    # Ring unsaturation composes with the tellone suffix locant on the
    # bicyclic/polycyclic branch, mirroring `_ketone.py`'s identical
    # extension -- reviewed rather than independently PubChem-confirmed.
    assert smiles_to_iupac("[Te]=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-tellone"


def test_unsaturated_monospiro_tellone_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC2(C1)C=CCCC2")
