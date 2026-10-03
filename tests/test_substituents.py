import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure, adjacency
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
        # NOT tert-butyl: root forks three ways with each arm still a
        # length-1 chain, but one arm is -CH2Br, not a bare methyl --
        # taking the tert-butyl shortcut here would silently drop the
        # bromine (found via real-data testing:
        # 'CC(C)c1cccc(C(C)(C)CBr)c1' was misnamed
        # '1-tert-butyl-3-(propan-2-yl)benzene').
        ("CCCC(C(C)(C)CBr)CCC", "4-(1-bromo-2-methylpropan-2-yl)heptane"),
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
        # a branch-point root forking three ways, not all methyls (unlike
        # tert-butyl above) -- found via real-data testing: the "extra"
        # third branch (not part of the two-branch spine) was counted
        # twice, once by the spine walk and once by a redundant separate
        # pass, over-counting it as a fourth methyl
        # ("2,2,4,4-tetramethylpentan-2-yl" instead of the correct
        # "2,4,4-trimethylpentan-2-yl", PubChem CID 252905's ester analog).
        (
            "CCCCCC(C(C)(C)CC(C)(C)C)CCCCC",
            "6-(2,4,4-trimethylpentan-2-yl)undecane",
        ),
    ],
)
def test_compound_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_disjoint_ring_pair_substituent_resolves():
    # Two separate cyclopropane rings joined by a chain (PubChem CID
    # 524617): a disjoint-ring-pair shape, not a cyclic substituent.
    assert smiles_to_iupac("C1CC1CCCCC1CC1") == "1,1'-(butane-1,4-diyl)dicyclopropane"


def test_polycyclic_substituent_raises():
    # A plain ring joined by a chain to a genuinely polycyclic (bicyclic,
    # not just a second plain monocyclic) ring system is still rejected.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1CC1CC2CCC1C2")


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


def test_heteroaromatic_ring_substituent():
    # A plain pyridine ring hanging off atom 0, attached at the ring carbon
    # two bonds from the nitrogen -- "pyridin-3-yl" (P-29.3.4.1), same
    # `name_branch` path as plain "phenyl" but recognizing the ring's own
    # heteroatom composition via `mol` instead of just its graph shape.
    mol = Chem.MolFromSmiles("Cc1cccnc1")
    graph = adjacency(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    assert name_branch(graph, 1, 0, aromatic_atoms=aromatic_atoms, mol=mol) == ("pyridin-3-yl", True)


def test_ring_substituent_with_own_hydroxyl():
    # A 6-membered ring (1-2-3-4-5-6-1) hanging off atom 0, attached at
    # ring atom 1, with a hydroxyl oxygen (7) on the ring atom directly
    # opposite the attachment point (position 4, unambiguous either way
    # around the ring); cross-checked end-to-end against PubChem CID 21395558 via
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
    # single-hydroxyl case above; cited together with a "di" multiplying prefix. Going one
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


def test_mol_defends_against_unvalidated_heteroatom_branch():
    # PR #487/#488 found the same bug class twice: a caller handed
    # `name_branch` a branch it hadn't fully validated as carbon-plus-
    # `halogens`, so an unrecognized heteroatom (here an amine nitrogen)
    # was silently walked as if it were an ordinary chain carbon, instead
    # of being rejected. `mol` closes this for any not-yet-audited caller.
    mol = Chem.MolFromSmiles("CCCN")
    graph = adjacency(mol)
    assert name_branch(graph, 1, 0, {}) == ("propyl", False)
    with pytest.raises(UnsupportedStructure):
        name_branch(graph, 1, 0, {}, mol=mol)


def test_mol_names_unsaturated_branch():
    # Without `mol` bond orders are invisible, so a double bond would vanish
    # silently; with it the branch is named with its 'ene' ending.
    mol = Chem.MolFromSmiles("CCC=C")
    graph = adjacency(mol)
    assert name_branch(graph, 1, 0, {}) == ("propyl", False)
    assert name_branch(graph, 1, 0, {}, mol=mol) == ("prop-2-en-1-yl", True)
    with pytest.raises(UnsupportedStructure):
        name_branch(graph, 1, 0, {}, mol=mol, unsaturated=False)


def test_halogenated_phenyl_substituent_alphabetical_tiebreak():
    # `halogenated_phenyl_substituent` picks whichever ring-walk direction
    # gives the lowest locant set (P-14.5.2); when a symmetric halogen
    # pattern makes both directions' locant sets identical ({2,4,6} either
    # way here), the direction was previously chosen arbitrarily by
    # ring-neighbor iteration order instead of by the required alphabetical
    # tiebreak (the alphabetically-first-cited substituent -- 'bromo' before
    # 'chloro' -- must get the lower locant). These two SMILES draw the same
    # molecule (one ortho Br, the other ortho Cl, para F) in opposite ring
    # directions and must produce the identical, correct name.
    expected = "1-(2-bromo-6-chloro-4-fluorophenyl)propan-2-one"
    assert smiles_to_iupac("CC(=O)Cc1c(Cl)cc(F)cc1Br") == expected
    assert smiles_to_iupac("CC(=O)Cc1c(Br)cc(F)cc1Cl") == expected


def test_halogenated_phenyl_substituent_non_tied_locant_set():
    # A non-symmetric case (the two walk directions' locant sets differ:
    # {2,3} vs {5,6}) must still resolve by lowest locant set alone,
    # unaffected by the alphabetical tiebreak added above. PubChem CID
    # 83724710.
    assert smiles_to_iupac("CC(=O)Cc1c(Cl)c(Br)ccc1") == "1-(3-bromo-2-chlorophenyl)propan-2-one"
