import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-65.1.5.2's own worked example for CH3-CH2-CH2-CH2-CH2-CS-OH:
        # 'hexanethioic O-acid (PIN)' -- the oxygen carries the acidic H.
        ("CCCCCC(=S)O", "hexanethioic O-acid"),
        # Blue Book P-65.1.5.2's own worked example for CH3-CS-OH:
        # 'ethanethioic O-acid (PIN)'.
        ("CC(=S)O", "ethanethioic O-acid"),
        # Blue Book P-65.1.5.2's own worked example for HCO-SH (formic acid
        # analogue, chain length 1, 'methane' stem not 'carbo-'):
        # 'methanethioic S-acid (PIN)' -- the sulfur carries the acidic H.
        ("O=CS", "methanethioic S-acid"),
        # Same tautomer direction as the formic case above, but with a
        # longer chain: the sulfur carries the acidic H, so 'S-acid'.
        ("CCCCCC(=O)S", "hexanethioic S-acid"),
    ],
)
def test_thioic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)S")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCCC1S")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)S")


def test_two_thioic_acid_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC(=O)CC(=O)S")
