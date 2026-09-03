import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 136810730: both drawn tautomers collapse to the same
        # structure and name (no O/Se letter, unlike the neutral acid) --
        # "CC(=O)[Se-]"/"CC(=[Se])[O-]" both -> "ethaneselenoate".
        ("CC(=O)[Se-]", "ethaneselenoate"),
        ("CC(=[Se])[O-]", "ethaneselenoate"),
        # Formic acid analogue (chain length 1) -- PubChem CID 148549672,
        # "methaneselenoate".
        ("[O-]C=[Se]", "methaneselenoate"),
    ],
)
def test_selenoate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_selenoate_with_halogen_substituent():
    # Not registered in PubChem (queried, CID 0) -- mechanical extension of
    # the PubChem-verified `_thioate.py` halogen case
    # ('ClCC(=O)[S-]' -> '2-chloroethanethioate') to the selenium analogue,
    # reviewed but not independently verified.
    assert smiles_to_iupac("ClCC(=O)[Se-]") == "2-chloroethaneselenoate"


def test_selenoate_with_double_bond():
    # Not registered in PubChem (queried, CID 0) -- mechanical extension of
    # the PubChem-verified `_thioate.py` alkene case
    # ('CC=CC(=O)[S-]' -> 'but-2-enethioate'), reviewed but not
    # independently verified.
    assert smiles_to_iupac("CC=CC(=O)[Se-]") == "but-2-eneselenoate"


def test_multiple_selenoate_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se-]C(=O)CCC(=O)[Se-]")


def test_selenoate_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se-]C(=O)C1CCCCC1")


def test_selenoate_branched_r_group():
    # Not registered in PubChem (queried, CID 0) -- mechanical extension of
    # the PubChem-verified `_thioate.py` branched case
    # ('CC(C)C(=O)[S-]' -> '2-methylpropanethioate'), reviewed but not
    # independently verified.
    assert smiles_to_iupac("CC(C)C(=O)[Se-]") == "2-methylpropaneselenoate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): the selenoate
        # carbon is always chain-terminal (fixed C1), same pattern as
        # `_thioate.py`/`_carboxylate.py`. Not registered in PubChem
        # (queried, CID 0) -- mechanical extension of the PubChem-verified
        # `_thioate.py` stereocenter case, reviewed but not independently
        # verified.
        ("CC[C@@H](C)C(=O)[Se-]", "(2R)-2-methylbutaneselenoate"),
        ("CC[C@H](C)C(=O)[Se-]", "(2S)-2-methylbutaneselenoate"),
        ("CC[C@@H](C)C(=[Se])[O-]", "(2R)-2-methylbutaneselenoate"),
    ],
)
def test_selenoate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_selenoate_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)C(=O)[Se-]") == "2-methylbutaneselenoate"


def test_phenyl_chain_selenoate():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_carboxylate.py`'s
    # phenyl-chain path): the ring is cited as a "phenyl" substituent
    # prefix. PubChem PUG REST: "2-phenylethaneselenoate".
    assert smiles_to_iupac("c1ccccc1CC(=O)[Se-]") == "2-phenylethaneselenoate"
    assert smiles_to_iupac("c1ccccc1CCC(=O)[Se-]") == "3-phenylpropaneselenoate"


def test_phenyl_directly_attached_selenoate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[Se-]")


def test_phenyl_substituted_benzene_ring_selenoate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)[Se-]")


def test_phenyl_chain_selenoate_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)[Se-]")
