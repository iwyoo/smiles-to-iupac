import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Formic acid analogue (chain length 1, 'methane' stem), Te-acid
        # direction (tellurium carries the acidic H), structure confirmed
        # on PubChem (CID 173268937, auto-generated name 'methanetelluroic
        # Te-acid').
        ("O=C[TeH]", "methanetelluroic Te-acid"),
        # Longer chain, same tautomer direction -- no PubChem entry exists
        # for this exact structure (CID 0), so this is a reviewed, not
        # directly verified, mechanical extension of the pattern confirmed
        # by the formic case above (see module docstring).
        ("CC(=O)[TeH]", "ethanetelluroic Te-acid"),
        # O-acid direction (the oxygen carries the acidic H, C=Te retained)
        # -- PubChem has a bare CID for this structure (85820967) but no
        # computed name, so this is also a reviewed mechanical extension.
        ("CC(=[Te])O", "ethanetelluroic O-acid"),
    ],
)
def test_telluroic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)[TeH]")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCCC1[TeH]")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)[TeH]")


def test_two_telluroic_acid_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C(=O)CC(=O)[TeH]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269-#281's
        # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/
        # acyl-halide/sulfonic-acid/thiol/sulfinic-acid/thioic-acid/
        # selenoic-acid chains): the ring is cited as a "phenyl"
        # substituent prefix. This module never supports any other
        # substituent, so the ring is always at the chain's far terminus
        # with no locant tie-break needed.
        ("c1ccccc1CC(=O)[TeH]", "2-phenylethanetelluroic Te-acid"),
        ("c1ccccc1CCC(=O)[TeH]", "3-phenylpropanetelluroic Te-acid"),
        ("c1ccccc1CC(=[Te])O", "2-phenylethanetelluroic O-acid"),
    ],
)
def test_phenyl_chain_telluroic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_telluroic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[TeH]")


def test_phenyl_substituted_benzene_ring_telluroic_acid_methyl():
    # A plain alkyl ring substituent is now supported alongside the chain
    # (mirroring `_thioic_acid.py`'s `plain_alkyl_ring_substituents`
    # rollout).
    assert smiles_to_iupac("Cc1ccccc1CC(=O)[TeH]") == "2-(2-methylphenyl)ethanetelluroic Te-acid"


def test_phenyl_substituent_telluroic_acid_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CC(=O)[TeH])cc1") == "2-(4-chlorophenyl)ethanetelluroic Te-acid"


def test_phenyl_substituent_telluroic_acid_branched_chain():
    assert (
        smiles_to_iupac("CC(C)Cc1ccc(cc1)C(C)C(=O)[TeH]")
        == "2-[4-(2-methylpropyl)phenyl]propanetelluroic Te-acid"
    )
