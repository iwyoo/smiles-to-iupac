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
