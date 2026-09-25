import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem open-chain (keto-) 2-ketose structures, D and L
        # series (#1039 M1 step 2).
        ("C([C@H](C(=O)CO)O)O", "D-erythrulose"),  # CID 5460177
        ("C([C@H]([C@H](C(=O)CO)O)O)O", "D-ribulose"),  # CID 151261
        ("C([C@H]([C@@H](C(=O)CO)O)O)O", "D-xylulose"),  # CID 5289590
        ("C([C@@H]([C@H](C(=O)CO)O)O)O", "L-xylulose"),  # CID 22253
        ("C([C@@H]([C@@H](C(=O)CO)O)O)O", "L-ribulose"),  # CID 644111
        ("C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O", "D-fructose"),  # CID 5984
        ("C([C@H]([C@@H]([C@H](C(=O)CO)O)O)O)O", "D-sorbose"),  # CID 107428
        ("C([C@H]([C@H]([C@H](C(=O)CO)O)O)O)O", "D-psicose"),  # CID 90008
        ("C([C@H]([C@@H]([C@@H](C(=O)CO)O)O)O)O", "D-tagatose"),  # CID 92092
        ("C([C@@H]([C@H]([C@@H](C(=O)CO)O)O)O)O", "L-sorbose"),  # CID 6904
        ("C([C@@H]([C@@H]([C@H](C(=O)CO)O)O)O)O", "L-fructose"),  # CID 5460024
        ("C([C@@H]([C@@H]([C@@H](C(=O)CO)O)O)O)O", "L-psicose"),  # CID 11961810
        ("C([C@@H]([C@H]([C@H](C(=O)CO)O)O)O)O", "L-tagatose"),  # CID 10965117
    ],
)
def test_open_chain_2_ketose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_dihydroxyacetone_still_falls_through_unchanged():
    # 1,3-dihydroxyacetone (PubChem CID 670) has no chirality center at
    # all (both sides of the carbonyl are plain -CH2OH), so there's no
    # D/L series to assign and this stays out of scope -- keeps the
    # existing generic substitutive name unchanged.
    assert smiles_to_iupac("OCC(=O)CO") == "1,3-dihydroxypropan-2-one"


def test_unspecified_stereo_2_ketose_still_falls_through_unchanged():
    # A 2-ketohexose with no stereo specified at all has no determinable
    # D/L series (P-102.3.2), so it keeps the existing generic
    # substitutive name unchanged, same as any other stereo-unspecified
    # molecule.
    assert smiles_to_iupac("C(C(C(C(C(=O)CO)O)O)O)O") == "1,3,4,5,6-pentahydroxyhexan-2-one"
