import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 92877, 5alpha-androstan-3beta-ol -- real, fully
        # stereo-specified structure; the recognized steroid skeleton
        # (androstane) plus its own fixed-numbering locant map identifies
        # the hydroxyl at C3.
        (
            "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4C[C@H](CC[C@@]4([C@H]3CC2)C)O",
            "5α-androstan-3β-ol",
        ),
        # PubChem CID 449196, 5alpha-androstan-3alpha-ol -- same skeleton
        # and locant, opposite configuration at the hydroxyl carbon itself
        # (which the steroid check strips away and re-derives without
        # relying on, since a steroid name doesn't cite this project's
        # stereodescriptors yet).
        (
            "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4C[C@@H](CC[C@@]4([C@H]3CC2)C)O",
            "5α-androstan-3α-ol",
        ),
        # PubChem CID 235071, 5alpha-androstan-17beta-ol -- hydroxyl at
        # the non-fusion D-ring position instead.
        (
            "C[C@]12CCCC[C@@H]1CC[C@@H]3[C@@H]2CC[C@]4([C@H]3CC[C@@H]4O)C",
            "5α-androstan-17β-ol",
        ),
        # PubChem CID 6665, 5alpha-cholestan-3beta-ol -- cholestane
        # skeleton (longest recognized side chain).
        (
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C",
            "5α-cholestan-3β-ol",
        ),
        # Unspecified stereochemistry is also accepted -- the steroid check
        # runs before the specified-stereocenter gate, but doesn't require
        # stereo to be present at all. PubChem CID 439222/623919/616360
        # (androstan-3-ol/androstan-17-ol/pregnan-3-ol) are each registered
        # under their generic name with no stereo specified at all.
        ("CC12CCCC1C3CCC4CC(CCC4(C3CC2)C)O", "androstan-3-ol"),
        ("CC12CCCCC1CCC3C2CCC4(C3CCC4O)C", "androstan-17-ol"),
        ("CCC1CCC2C1(CCC3C2CCC4C3(CCC(C4)O)C)C", "pregnan-3-ol"),
    ],
)
def test_steroid_alcohol_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_non_steroid_von_baeyer_alcohol_with_stereo_now_resolves():
    # A specified stereocenter on a non-steroid von Baeyer polycyclic
    # alcohol is no longer rejected (#1079, M5 step 2) -- the steroid
    # check only ever bypassed the stereo gate for a real steroid
    # skeleton match; the non-steroid path now cites the stereocenter
    # instead of rejecting it, same as this module's steroid branch
    # already did.
    assert smiles_to_iupac("O[C@H]1C[C@@H]2CC[C@]1(C)C2") == "(1R,2S,4R)-1-methylbicyclo[2.2.1]heptan-2-ol"
