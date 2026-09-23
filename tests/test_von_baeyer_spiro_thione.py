import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 12859496, 126746952 respectively). The
        # parent stem's final 'e' is kept ('thione' begins with a
        # consonant, P-16.3.3, same as `_thiol.py`'s 'thiol').
        ("S=C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2-thione"),
        ("S=C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-thione"),
    ],
)
def test_von_baeyer_spiro_thione_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No further real PubChem-registered polycyclic/spiro thione
        # structure was found beyond the two above -- a mechanical
        # single-axis extension of the identical, already-independently-
        # verified mechanism (`_alcohol.py`/`_amine.py`/`_ketone.py`/
        # `_thiol.py`/`_sulfinamide.py`/`_selone.py`/`_sulfonamide.py`/
        # `_sulfinic_acid.py`/`_sulfonic_acid.py` all share it), reviewed
        # rather than independently PubChem-confirmed (see
        # `implementation-notes.md`'s test-writing policy).
        ("S=C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-thione"),
        ("S=C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-thione"),
    ],
)
def test_von_baeyer_spiro_thione_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_thiones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CC2CCC1C(=S)C2")


def test_unsaturated_von_baeyer_thione_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CC2C=CC1C2")
