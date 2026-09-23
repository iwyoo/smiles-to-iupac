import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 10345, 107363, 139709, 276875, 137688,
        # 12550422 respectively).
        ("O=C1CC2CCC1C2", "bicyclo[2.2.1]heptan-2-one"),
        ("O=C1CCC2CCC1C2", "bicyclo[3.2.1]octan-2-one"),
        ("O=C1CC2CCC(C1)C2", "bicyclo[3.2.1]octan-3-one"),
        ("O=C1CC2CCCC(C1)C2", "bicyclo[3.3.1]nonan-3-one"),
        ("O=C1CC2CCC1CC2", "bicyclo[2.2.2]octan-2-one"),
        ("O=C1CC2CCC2C1", "bicyclo[3.2.0]heptan-3-one"),
        # 1-adamantanone -- this project's adamantane already uses the
        # systematic 'tricyclo[3.3.1.1^3,7]decane' name rather than the
        # retained 'adamantane' one (see `test_von_baeyer_alcohol.py`'s
        # identical 1-/2-adamantanol precedent), so this module's
        # ring_count>=3 path follows the same systematic convention
        # rather than PubChem's own retained-name 'adamantan-2-one'.
        ("O=C1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.1^3,7]decan-2-one"),
    ],
)
def test_von_baeyer_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_ketones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2CCC1C(=O)C2")


def test_ketone_hydroxyl_combination_on_ring_system_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2CCC1C2O")


def test_unsaturated_von_baeyer_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2C=CC1C2")
