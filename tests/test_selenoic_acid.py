import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # S-acid-direction tautomer (the selenium carries the acidic H),
        # structure confirmed on PubChem (CID 57348413, auto-generated name
        # 'ethaneselenoic Se-acid').
        ("CC(=O)[SeH]", "ethaneselenoic Se-acid"),
        # O-acid-direction tautomer (the oxygen carries the acidic H, C=Se
        # retained) -- no PubChem-computed name exists for this tautomer
        # (its CID has only a bare connectivity record), so this is a
        # reviewed, not directly verified, mechanical extension of the
        # identical labeling rule confirmed in `_thioic_acid.py` (see
        # module docstring).
        ("CC(=[Se])O", "ethaneselenoic O-acid"),
        # Formic acid analogue (chain length 1, 'methane' stem), same
        # tautomer direction as the first case above -- a reviewed
        # mechanical extension (PubChem has no entry at all for this exact
        # structure).
        ("O=C[SeH]", "methaneselenoic Se-acid"),
        # Longer chain, S-acid direction -- reviewed mechanical extension of
        # the first case (PubChem has no entry for this exact structure).
        ("CCCCCC(=O)[SeH]", "hexaneselenoic Se-acid"),
    ],
)
def test_selenoic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)[SeH]")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCCC1[SeH]")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)[SeH]")


def test_two_selenoic_acid_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C(=O)CC(=O)[SeH]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269-#280's
        # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/
        # acyl-halide/sulfonic-acid/thiol/sulfinic-acid/thioic-acid
        # chains): the ring is cited as a "phenyl" substituent prefix.
        # This module never supports any other substituent, so the ring
        # is always at the chain's far terminus with no locant tie-break
        # needed.
        ("c1ccccc1CC(=O)[SeH]", "2-phenylethaneselenoic Se-acid"),
        ("c1ccccc1CCC(=O)[SeH]", "3-phenylpropaneselenoic Se-acid"),
        ("c1ccccc1CC(=[Se])O", "2-phenylethaneselenoic O-acid"),
    ],
)
def test_phenyl_chain_selenoic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_selenoic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[SeH]")


def test_phenyl_substituted_benzene_ring_selenoic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)[SeH]")
