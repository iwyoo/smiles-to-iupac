import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem D-aldohexopyranose structures, both anomers,
        # spanning all 8 retained hexose stem names (#1039 M2 step 1).
        ("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", "alpha-D-glucopyranose"),  # CID 79025
        ("C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O", "beta-D-glucopyranose"),  # CID 64689
        ("C([C@@H]1[C@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O)O", "alpha-D-mannopyranose"),  # CID 185698
        (
            "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O",
            "beta-D-galactopyranose",
        ),  # CID 439353
        ("C([C@@H]1[C@H]([C@H]([C@H]([C@H](O1)O)O)O)O)O", "alpha-D-allopyranose"),  # CID 3034742
        ("C([C@@H]1[C@H]([C@H]([C@H]([C@@H](O1)O)O)O)O)O", "beta-D-allopyranose"),  # CID 448388
        (
            "C([C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)O)O)O)O)O",
            "alpha-D-altropyranose",
        ),  # CID 7098663
        ("C([C@@H]1[C@H]([C@H]([C@@H]([C@@H](O1)O)O)O)O)O", "beta-D-altropyranose"),  # CID 448702
        ("C([C@@H]1[C@@H]([C@H]([C@H]([C@H](O1)O)O)O)O)O", "alpha-D-gulopyranose"),  # CID 7044038
        ("C([C@@H]1[C@@H]([C@H]([C@H]([C@@H](O1)O)O)O)O)O", "beta-D-gulopyranose"),  # CID 6102790
        ("C([C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O)O", "alpha-D-idopyranose"),  # CID 7098664
        ("C([C@@H]1[C@@H]([C@H]([C@@H]([C@@H](O1)O)O)O)O)O", "beta-D-idopyranose"),  # CID 7018164
        ("C([C@@H]1[C@@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O)O", "alpha-D-talopyranose"),  # CID 81696
        (
            "C([C@@H]1[C@@H]([C@@H]([C@@H]([C@@H](O1)O)O)O)O)O",
            "beta-D-talopyranose",
        ),  # CID 5319264
    ],
)
def test_cyclic_aldopyranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hetero_ring_ketone_still_resolves():
    # A plain saturated single-heteroatom ring bearing an actual ring
    # ketone (not a pyranose hemiacetal) must still route to `_ketone.py`
    # unchanged.
    assert smiles_to_iupac("O=C1CCCCO1") == "oxan-2-one"


def test_furanose_still_falls_through_unchanged():
    # A 5-membered cyclic hemiacetal (furanose) is a different ring size,
    # out of scope here (later M2 step) -- still collides with
    # `_ketone.py`'s hetero-ring-ketone path unchanged, same as a pyranose
    # did before this feature.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C([C@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O")
