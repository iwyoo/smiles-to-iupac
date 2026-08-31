import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._substituents import name_branch


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # isopropyl branch: root forks into two equal-length methyls, so
        # (P-29.3.2.2) the free valence sits internal to its own principal
        # chain -- propane, at locant 2 -- rather than at a chain terminus.
        ("CCCC(C(C)C)CCC", "4-(propan-2-yl)heptane"),
        # tert-butyl branch: root forks into three methyls -- P-29.6.1's
        # retained name, not the general rule's own name.
        ("CCCC(C(C)(C)C)CCC", "4-tert-butylheptane"),
        # sec-butyl branch: root forks into a methyl and an ethyl
        # continuation; the longer (ethyl) side wins the chain, giving
        # butane with the free valence at locant 2 (P-29.2: lowest locant
        # consistent with the chain).
        ("CCCCC(C(C)CC)CCCCC", "5-(butan-2-yl)decane"),
        # a compound substituent alongside a simple one on the same chain;
        # alphanumerical order (P-14.5.2) puts 'methyl' before the compound
        # substituent's own key ('propanyl'), so it gets the lower locant.
        ("CCCC(C)C(C(C)C)CCC", "4-methyl-5-(propan-2-yl)octane"),
        # a compound substituent on a ring, alongside a simple one;
        # alphanumerical order now puts the compound substituent's own key
        # ('butanyl') before 'methyl', so it gets the lower locant.
        ("CC1CCCCC1C(C)CC", "1-(butan-2-yl)-2-methylcyclohexane"),
        # two identical compound substituents: 'bis', not 'di' (P-14.2.2).
        (
            "CCCCC(C(C)CC)CCCCC(C(C)CC)CCCCC",
            "5,10-bis(butan-2-yl)pentadecane",
        ),
    ],
)
def test_compound_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_substituent_raises():
    # Two separate cyclopropane rings joined by a chain: more than one ring
    # overall, so this is rejected before a cyclic substituent could even be
    # considered (see smiles_to_iupac.core's polycyclic check).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC1CCCCC1CC1")


def test_simple_ring_substituent():
    # A plain, unsubstituted 3-membered ring (1-2-3-1) hanging off atom 0 is
    # named directly as "cyclopropyl" (see tasks/ring-substituent-chain-
    # suffix.md, 2026-08-25) rather than rejected as a cyclic substituent.
    graph = {0: [1], 1: [0, 2, 3], 2: [1, 3], 3: [1, 2]}
    assert name_branch(graph, 1, 0) == ("cyclopropyl", False)


def test_substituted_ring_substituent_raises():
    # Same 3-membered ring, but one ring atom (2) also carries its own
    # exocyclic branch (4) -- no longer the plain unsubstituted shape
    # `_simple_ring_substituent` recognizes, so this still falls through to
    # the ordinary chain-walk's cycle-detection rejection (P-29.3.3).
    graph = {0: [1], 1: [0, 2, 3], 2: [1, 3, 4], 3: [1, 2], 4: [2]}
    with pytest.raises(UnsupportedStructure):
        name_branch(graph, 1, 0)


def test_ring_substituent_with_own_hydroxyl():
    # A 6-membered ring (1-2-3-4-5-6-1) hanging off atom 0, attached at
    # ring atom 1, with a hydroxyl oxygen (7) on the ring atom directly
    # opposite the attachment point (position 4, unambiguous either way
    # around the ring) -- see tasks/ring-substituent-own-hydroxyl-naming.md,
    # 2026-08-26; cross-checked end-to-end against PubChem CID 21395558 via
    # the full molecule in tests/test_alcohol.py.
    graph = {0: [1], 1: [0, 2, 6], 2: [1, 3], 3: [2, 4], 4: [3, 5, 7], 5: [4, 6], 6: [5, 1], 7: [4]}
    assert name_branch(graph, 1, 0, {7: "hydroxy"}) == ("4-hydroxycyclohexyl", True)


def test_ring_substituent_with_own_hydroxyl_picks_lower_locant():
    # Same ring, hydroxyl on the atom immediately adjacent to the
    # attachment point -- going one way around gives it locant 2, the
    # other way locant 6; the lower one (2) must win.
    graph = {0: [1], 1: [0, 2, 6], 2: [1, 3, 7], 3: [2, 4], 4: [3, 5], 5: [4, 6], 6: [5, 1], 7: [2]}
    assert name_branch(graph, 1, 0, {7: "hydroxy"}) == ("2-hydroxycyclohexyl", True)


def test_ring_substituent_with_two_hydroxyls():
    # Two hydroxyls on the same ring substituent (on ring atoms 2 and 4,
    # relative to the attachment at atom 1) -- generalized from the
    # single-hydroxyl case above by tasks/ring-vs-chain-alcohol-multi-oh.md,
    # 2026-08-28; cited together with a "di" multiplying prefix. Going one
    # way around the ring gives locants (2, 4), the other way (4, 6); the
    # lower set (2, 4) must win -- hand-verified, no PubChem cross-check at
    # this atom-index level (see tests/test_alcohol.py for an end-to-end
    # molecule-level case).
    graph = {
        0: [1],
        1: [0, 2, 6],
        2: [1, 3, 7],
        3: [2, 4],
        4: [3, 5, 8],
        5: [4, 6],
        6: [5, 1],
        7: [2],
        8: [4],
    }
    assert name_branch(graph, 1, 0, {7: "hydroxy", 8: "hydroxy"}) == ("2,4-dihydroxycyclohexyl", True)
