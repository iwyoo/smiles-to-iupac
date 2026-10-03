import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 13239288, 13239278, 98480, 13311924,
        # 21522752, 53402094 respectively).
        ("NC1CC2CCC1CC2", "bicyclo[2.2.2]octan-2-amine"),
        ("NC1CCC2CCC1C2", "bicyclo[3.2.1]octan-2-amine"),
        ("NC1CC2CCC1C2", "bicyclo[2.2.1]heptan-2-amine"),
        ("NC1CC2CCC(C1)C2", "bicyclo[3.2.1]octan-3-amine"),
        ("NC1CC2CCCC(C1)C2", "bicyclo[3.3.1]nonan-3-amine"),
        ("NC1CC2CCC2C1", "bicyclo[3.2.0]heptan-3-amine"),
        # 1-adamantanamine -- this project's adamantane already uses the
        # systematic 'adamantane' name rather than the
        # retained 'adamantane' one (see `test_von_baeyer_alcohol.py`'s
        # identical 1-/2-adamantanol precedent), so this module's
        # ring_count>=3 path follows the same systematic convention
        # rather than PubChem's own retained-name 'adamantan-1-amine'.
        ("NC12CC3CC(CC(C3)C1)C2", "adamantan-1-amine"),
    ],
)
def test_von_baeyer_amine_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CC2CCC1C2N")


def test_amine_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC1CC2CCC1C2")


def test_von_baeyer_amine_ring_unsaturation():
    # Ring unsaturation composes with the amine suffix locant, mirroring
    # `_ketone.py`'s identical extension -- PubChem CID 12841951, matches
    # its own IUPACName exactly.
    assert smiles_to_iupac("C1C2CC(C1C=C2)N") == "bicyclo[2.2.1]hept-5-en-2-amine"


def test_von_baeyer_amine_specified_stereocenter():
    # PubChem CID 10240785; matches its own IUPACName exactly.
    assert smiles_to_iupac("C1C[C@@H]2C[C@H]1C[C@H]2N") == "(1R,2R,4S)-bicyclo[2.2.1]heptan-2-amine"
