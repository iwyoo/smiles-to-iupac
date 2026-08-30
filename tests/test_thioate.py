import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-72.2.2.2.1.1 worked examples, both drawn tautomers -- the
        # primary source gives both the *same* name (no O/S letter,
        # unlike the neutral acid): "CH3-CH2-CO-S- <-> CH3-CH2-CS-O-" ->
        # "propanethioate (PIN)", "CH3-CO-S- <-> CH3-CS-O-" ->
        # "ethanethioate (PIN)".
        ("CCC(=O)[S-]", "propanethioate"),
        ("CCC(=S)[O-]", "propanethioate"),
        ("CC(=O)[S-]", "ethanethioate"),
        ("CC(=S)[O-]", "ethanethioate"),
        # Formic acid analogue (chain length 1) uses the 'methane' stem --
        # PubChem CID 10905395, "methanethioate".
        ("[O-]C=S", "methanethioate"),
    ],
)
def test_thioate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_thioate_with_halogen_substituent():
    # PubChem CID 23063092: 2-chloroethanethioate ('ClCC(=O)[S-]').
    assert smiles_to_iupac("ClCC(=O)[S-]") == "2-chloroethanethioate"


def test_thioate_with_double_bond():
    # PubChem CID 149437485: but-2-enethioate ('CC=CC(=O)[S-]') -- no
    # stereo marker, since E/Z on the module's own double bond is out of
    # scope here (mirrors `_carboxylate.py`'s identical lack of E/Z
    # handling).
    assert smiles_to_iupac("CC=CC(=O)[S-]") == "but-2-enethioate"


def test_multiple_thioate_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[S-]C(=O)CCC(=O)[S-]")


def test_thioate_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[S-]C(=O)C1CCCCC1")


def test_thioate_branched_r_group():
    # PubChem CID 20063546: 2-methylpropanethioate.
    assert smiles_to_iupac("CC(C)C(=O)[S-]") == "2-methylpropanethioate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): the thioate
        # carbon is always chain-terminal (fixed C1), same pattern as
        # `_carboxylate.py`. PubChem CID 154057341.
        ("CC[C@@H](C)C(=O)[S-]", "(2R)-2-methylbutanethioate"),
        ("CC[C@H](C)C(=O)[S-]", "(2S)-2-methylbutanethioate"),
        # The other drawn tautomer gives the identical name (module
        # docstring: the anion charge is delocalized across both
        # chalcogens).
        ("CC[C@@H](C)C(=S)[O-]", "(2R)-2-methylbutanethioate"),
    ],
)
def test_thioate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_thioate_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)C(=O)[S-]") == "2-methylbutanethioate"
