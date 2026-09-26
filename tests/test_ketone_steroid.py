import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 102028, 5alpha-androstan-3-one -- real, fully
        # stereo-specified structure; the recognized steroid skeleton
        # (androstane) plus its own fixed-numbering locant map identifies
        # the ketone at C3.
        (
            "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CC(=O)CC[C@@]4([C@H]3CC2)C",
            "androstan-3-one",
        ),
        # PubChem CID 101927, 5alpha-androstan-17-one -- same skeleton,
        # ketone at C17 instead (the non-fusion D-ring position).
        (
            "C[C@]12CCCC[C@@H]1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CCC4=O)C",
            "androstan-17-one",
        ),
        # PubChem CID 22213548, 5alpha-pregnan-3-one -- pregnane skeleton
        # (bigger side chain, exercises constitution matching that ignores
        # the side chain's own stereocenter).
        (
            "CC[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCC(=O)C4)C)C",
            "pregnan-3-one",
        ),
        # PubChem CID 92128, 5alpha-cholestan-3-one -- cholestane skeleton
        # (longest recognized side chain).
        (
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCC(=O)C4)C)C",
            "cholestan-3-one",
        ),
        # Unspecified stereochemistry is also accepted -- the steroid check
        # runs before the specified-stereocenter gate, but doesn't require
        # stereo to be present at all.
        ("CC12CCC3C(C1CCC2)CCC4CCC(=O)CC34C", "androstan-2-one"),
    ],
)
def test_steroid_ketone_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hydroxyl_alongside_steroid_ketone_still_raises():
    # Testosterone: androst-4-en-17-ol-3-one -- a coexisting hydroxyl is
    # still out of scope (WS2/M2, #1029), and the ring's own C4-C5
    # unsaturation is separately out of scope too; either gate must still
    # fire exactly as before this change.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34")


def test_non_steroid_von_baeyer_ketone_with_stereo_now_resolves():
    # A specified stereocenter on a non-steroid von Baeyer polycyclic
    # ketone is no longer rejected (#1078, M5 step 1) -- the steroid check
    # only ever bypassed the stereo gate for a real steroid skeleton
    # match; the non-steroid path now cites the stereocenter instead of
    # rejecting it, same as this module's steroid branch already did.
    assert smiles_to_iupac("O=C1C[C@@H]2CC[C@H]1[C@@H]2C") == "(1S,4S,7R)-7-methylbicyclo[2.2.1]heptan-2-one"


def test_non_steroid_von_baeyer_ketone_without_stereo_still_resolves():
    assert smiles_to_iupac("O=C1CC2CCC1CC2") == "bicyclo[2.2.2]octan-2-one"
