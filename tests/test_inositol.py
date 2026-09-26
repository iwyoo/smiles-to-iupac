import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Each isomeric SMILES was built directly from P-104.2.1's own
        # face-notation fraction (numerator locants' hydroxyl 'up', the
        # rest 'down', ring numbered clockwise as viewed from the 'up'
        # face) and cross-checked against P-104.2.3's own worked
        # '1L-1,2/3,5-cyclohexanetetrol' == '(1R,2R,3R,5R)-cyclohexane-
        # 1,2,3,5-tetrol' example before being trusted here (module
        # docstring).
        ("O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O", "cis-inositol"),  # (1,2,3,4,5,6/0)
        ("O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O", "epi-inositol"),  # (1,2,3,4,5/6-)
        ("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O", "allo-inositol"),  # (1,2,3,4/5,6-)
        ("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "myo-inositol"),  # (1,2,3,5/4,6-)
        ("O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O", "muco-inositol"),  # (1,2,4,5/3,6-)
        ("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1O", "neo-inositol"),  # (1,2,3/4,5,6-)
        ("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "scyllo-inositol"),  # (1,3,5/2,4,6-)
        ("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O", "1L-chiro-inositol"),  # (1,2,4/3,5,6-)
        ("O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O", "1D-chiro-inositol"),  # (1,2,4/3,5,6-)
    ],
)
def test_inositol_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_non_stereo_hexol_uses_systematic_name():
    # No stereo specified at all -- must keep resolving to the generic
    # systematic name, not any retained inositol name.
    assert smiles_to_iupac("OC1C(O)C(O)C(O)C(O)C1O") == "cyclohexane-1,2,3,4,5,6-hexol"


def test_partially_specified_hexol_does_not_match_retained_name():
    # Only some stereocenters specified -- doesn't match any of the 9
    # fully-specified retained-name entries, so it must not be misnamed
    # as one of them.
    assert smiles_to_iupac("O[C@H]1C(O)C(O)C(O)C(O)C1O") != "myo-inositol"
