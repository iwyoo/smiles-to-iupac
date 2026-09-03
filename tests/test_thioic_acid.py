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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269-#279's
        # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/
        # acyl-halide/sulfonic-acid/thiol/sulfinic-acid chains): the ring
        # is cited as a "phenyl" substituent prefix. This module never
        # supports any other substituent, so the ring is always at the
        # chain's far terminus with no locant tie-break needed.
        ("c1ccccc1CC(=O)S", "2-phenylethanethioic S-acid"),
        ("c1ccccc1CCC(=O)S", "3-phenylpropanethioic S-acid"),
        ("c1ccccc1CC(=S)O", "2-phenylethanethioic O-acid"),
    ],
)
def test_phenyl_chain_thioic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_thioic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)S")


def test_phenyl_substituted_benzene_ring_thioic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)S")
