import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Same skeletons/locants as `test_von_baeyer_carbenium.py`'s
        # cases, with the ring carbon itself carrying a monovalent radical
        # (one fewer hydrogen) instead of a +1 charge -- these radical
        # structures aren't reliably registered/named by PubChem (its own
        # auto-namer drops the radical entirely and returns the plain
        # parent-hydride name, see `_radical.py`'s wiring commit), so
        # verification here cross-checks the neutral parent-hydride
        # skeleton name plus the shared suffix-locant mechanism reused
        # unchanged from `_carbenium.py`.
        ("[CH]1CC2CCC1CC2", "bicyclo[2.2.2]octan-2-yl"),
        ("[CH]1CCC2CCC1C2", "bicyclo[3.2.1]octan-2-yl"),
        ("[CH]1CC2CCC1C2", "bicyclo[2.2.1]heptan-2-yl"),
        ("[CH]1CC2CCC(C1)C2", "bicyclo[3.2.1]octan-3-yl"),
        ("[CH]1CC2CCCC(C1)C2", "bicyclo[3.3.1]nonan-3-yl"),
        ("[CH]1CC2CCC2C1", "bicyclo[3.2.0]heptan-3-yl"),
        # A bridgehead radical (degree 3, no hydrogen) rather than a
        # secondary ring position -- same tricyclic (adamantane) skeleton
        # as `test_von_baeyer_carbenium.py`'s bridgehead case.
        ("[C]12CC3CC(CC(C3)C1)C2", "tricyclo[3.3.1.1^3,7]decan-1-yl"),
    ],
)
def test_von_baeyer_radical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_divalent_von_baeyer_radical_idene():
    # P-71.2.2.1: 'idene' is appended after the already-computed '-yl'
    # name for a divalent radical center on a von Baeyer ring system too.
    assert smiles_to_iupac("[C]1CC2CCC1C2") == "bicyclo[2.2.1]heptan-2-ylidene"


def test_multiple_ring_radical_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C1CC2CCC1C2[CH2]")


def test_radical_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C1CC2CCC1C2")


def test_unsaturated_von_baeyer_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH]1CC2C=CC1C2")


def test_stereo_von_baeyer_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH]1C[C@@H]2CCC1C2")
