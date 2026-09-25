import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem open-chain (aldehydo-) aldose structures, D and L
        # series, trioses through hexoses (#1039 M1 step 1).
        ("C([C@H](C=O)O)O", "D-glyceraldehyde"),  # CID 79014
        ("C([C@@H](C=O)O)O", "L-glyceraldehyde"),  # CID 439723
        ("C([C@H]([C@H](C=O)O)O)O", "D-erythrose"),  # CID 94176
        ("C([C@H]([C@@H](C=O)O)O)O", "D-threose"),  # CID 439665
        ("C([C@@H]([C@@H](C=O)O)O)O", "L-erythrose"),  # CID 5460674
        ("C([C@@H]([C@H](C=O)O)O)O", "L-threose"),  # CID 5460672
        ("C([C@H]([C@H]([C@H](C=O)O)O)O)O", "D-ribose"),  # CID 5311110
        ("C([C@H]([C@@H]([C@H](C=O)O)O)O)O", "D-xylose"),  # CID 644160
        ("C([C@H]([C@H]([C@@H](C=O)O)O)O)O", "D-arabinose"),  # CID 66308
        ("C([C@H]([C@@H]([C@@H](C=O)O)O)O)O", "D-lyxose"),  # CID 65550
        ("C([C@@H]([C@@H]([C@H](C=O)O)O)O)O", "L-arabinose"),  # CID 5460291
        ("C([C@@H]([C@H]([C@@H](C=O)O)O)O)O", "L-xylose"),  # CID 95259
        (
            "C([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O",
            "D-glucose",
        ),  # CID 107526
        (
            "C([C@H]([C@H]([C@@H]([C@@H](C=O)O)O)O)O)O",
            "D-mannose",
        ),  # CID 161658
        (
            "C([C@H]([C@@H]([C@@H]([C@H](C=O)O)O)O)O)O",
            "D-galactose",
        ),  # CID 3037556
        ("C([C@H]([C@H]([C@H]([C@H](C=O)O)O)O)O)O", "D-allose"),  # CID 102288
        (
            "C([C@H]([C@H]([C@H]([C@@H](C=O)O)O)O)O)O",
            "D-altrose",
        ),  # CID 94780
        (
            "C([C@H]([C@@H]([C@H]([C@H](C=O)O)O)O)O)O",
            "D-gulose",
        ),  # CID 167792
        ("C([C@H]([C@@H]([C@H]([C@@H](C=O)O)O)O)O)O", "D-idose"),  # CID 111123
        (
            "C([C@H]([C@@H]([C@@H]([C@@H](C=O)O)O)O)O)O",
            "D-talose",
        ),  # CID 99459
        (
            "C([C@@H]([C@@H]([C@H]([C@H](C=O)O)O)O)O)O",
            "L-mannose",
        ),  # CID 82308
    ],
)
def test_open_chain_aldose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ketose_still_falls_through_unchanged():
    # keto-D-fructose (PubChem CID 5984) -- a ketose, not an aldose, so
    # it's still out of scope here (later M1 step, #1041) and must keep
    # getting the existing generic substitutive ketone name unchanged.
    assert smiles_to_iupac("C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O") == (
        "(3S,4R,5R)-1,3,4,5,6-pentahydroxyhexan-2-one"
    )


def test_heptose_still_falls_through_unchanged():
    # D-glycero-D-gluco-heptose (PubChem CID 87131842) -- 7 carbons, out
    # of scope here (P-102.5.1.1.2's multi-prefix systematic composition,
    # later M1 step, #1042).
    assert smiles_to_iupac("C([C@H]([C@H]([C@H]([C@@H]([C@H](C=O)O)O)O)O)O)O") == (
        "(2R,3S,4R,5R,6R)-2,3,4,5,6,7-hexahydroxyheptanal"
    )


def test_unspecified_stereo_aldose_still_falls_through_unchanged():
    # A hexose with no stereo specified at all has no determinable D/L
    # series (P-102.3.2), so it keeps the existing generic substitutive
    # name unchanged, same as any other stereo-unspecified molecule.
    assert smiles_to_iupac("C(C(C(C(C(C=O)O)O)O)O)O") == "2,3,4,5,6-pentahydroxyhexanal"
