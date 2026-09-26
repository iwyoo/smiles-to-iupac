import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 45080580, 45080610, 150349005, 114820018
        # respectively). The parent stem's final 'e' is kept throughout
        # ('sulfonamide' begins with a consonant, P-16.3.3, same as
        # `_thiol.py`'s 'thiol').
        ("NS(=O)(=O)C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2-sulfonamide"),
        ("NS(=O)(=O)C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2-sulfonamide"),
        ("NS(=O)(=O)C1CC2CCC(C1)C2", "bicyclo[3.2.1]octane-3-sulfonamide"),
        ("NS(=O)(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfonamide"),
    ],
)
def test_von_baeyer_spiro_sulfonamide_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_sulfonamide_on_polycyclic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNS(=O)(=O)C1CC2CCC1C2")


def test_multiple_ring_sulfonamides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)C1CC2CCC1C2S(N)(=O)=O")


def test_von_baeyer_sulfonamide_ring_unsaturation():
    # Ring unsaturation composes with the sulfonamide suffix locant on
    # the bicyclic/polycyclic branch, mirroring `_ketone.py`'s identical
    # extension -- reviewed rather than independently PubChem-confirmed.
    assert smiles_to_iupac("NS(=O)(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfonamide"


def test_unsaturated_monospiro_sulfonamide_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)C1CCCC2(C1)C=CCCC2")
