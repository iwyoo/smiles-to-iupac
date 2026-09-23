import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_von_baeyer_selone_name():
    # PubChem-confirmed: CID 175998577, 'bicyclo[2.2.1]heptane-2-selone'
    # -- the parent stem's final 'e' is kept ('selone' begins with a
    # consonant, P-16.3.3, same as `_thiol.py`'s 'thiol').
    assert smiles_to_iupac("[Se]=C1CC2CCC1C2") == "bicyclo[2.2.1]heptane-2-selone"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No further real PubChem-registered polycyclic/spiro selone
        # structure was found beyond the one PubChem-confirmed case above
        # (selones are rare in the registry at this ring-system scale) --
        # these are a mechanical single-axis extension of the identical,
        # already-independently-verified mechanism `_alcohol.py`/
        # `_amine.py`/`_ketone.py`/`_thiol.py`/`_sulfinamide.py` all share
        # (only the suffix word itself varies), reviewed rather than
        # independently PubChem-confirmed (see `implementation-notes.md`'s
        # test-writing policy).
        ("[Se]=C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-selone"),
        ("[Se]=C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-selone"),
    ],
)
def test_von_baeyer_spiro_selone_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_selones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CC2CCC1C(=[Se])C2")


def test_unsaturated_von_baeyer_selone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CC2C=CC1C2")
