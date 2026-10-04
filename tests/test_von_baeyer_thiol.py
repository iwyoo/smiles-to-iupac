import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 164174492, 160890585, 156181711,
        # 146941873 respectively). The parent stem's final 'e' is kept
        # throughout ('thiol' begins with a consonant, P-16.3.3), unlike
        # `_alcohol.py`'s/`_amine.py`'s/`_ketone.py`'s vowel-initial
        # suffixes.
        ("SC1CC2CCC1C2", "bicyclo[2.2.1]heptane-2-thiol"),
        ("SC1CCC2CCC1C2", "bicyclo[3.2.1]octane-2-thiol"),
        ("SC1CC2CCC(C1)C2", "bicyclo[3.2.1]octane-3-thiol"),
        ("SC1CC2CCCC(C1)C2", "bicyclo[3.3.1]nonane-3-thiol"),
        ("SC1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-thiol"),
        ("SC1CC2CCC2C1", "bicyclo[3.2.0]heptane-3-thiol"),
        # 1-adamantanethiol -- this project's adamantane already uses the
        # systematic 'adamantane' name rather than the
        # retained 'adamantane' one (see `test_von_baeyer_alcohol.py`'s
        # identical 1-/2-adamantanol precedent), so this module's
        # ring_count>=3 path follows the same systematic convention
        # rather than PubChem's own retained-name 'adamantane-1-thiol'.
        ("SC12CC3CC(CC(C3)C1)C2", "adamantane-1-thiol"),
    ],
)
def test_von_baeyer_thiol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_thiols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CC2CCC1C2S")


def test_thiol_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("SCC1CC2CCC1C2") == "(bicyclo[2.2.1]heptan-2-yl)methanethiol"


def test_von_baeyer_thiol_ring_unsaturation():
    # Ring unsaturation composes with the thiol suffix locant, mirroring
    # `_ketone.py`'s identical extension -- PubChem CID 73114971, matches
    # its own IUPACName exactly.
    assert smiles_to_iupac("C1C2CC(C1C=C2)S") == "bicyclo[2.2.1]hept-5-ene-2-thiol"
