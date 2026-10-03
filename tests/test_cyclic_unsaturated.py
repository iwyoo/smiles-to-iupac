import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Plain unsubstituted cycloalkenes/dienes are well-known, unambiguous
        # compound names (cyclohexene, cyclohexadiene, cyclopentene,
        # cycloheptene); the substituted variants below are covered by
        # in-line locant-rule citations instead of a per-case lookup.
        ("C1=CCCCC1", "cyclohexene"),
        ("C1=CC=CCC1", "cyclohexa-1,3-diene"),
        ("C1=CCC=CC1", "cyclohexa-1,4-diene"),
        ("C1=CCCC1", "cyclopentene"),
        ("C1=CCCCCC1", "cycloheptene"),
    ],
)
def test_unsubstituted_cyclic_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_single_ene_locant_always_omitted_with_substituent():
    # methyl on the ring carbon three bonds from the double bond: numbering
    # is chosen so the ene is at the (always omitted) position '1' and the
    # substituent gets the lowest achievable locant, '3'.
    assert smiles_to_iupac("C1=CC(C)CCC1") == "3-methylcyclohexene"


def test_substituent_on_ene_carbon_gets_lowest_locant():
    # methyl on a double-bond carbon itself: numbering starts there so both
    # the ene (implicit '1') and methyl land on locant '1', not '2'.
    assert smiles_to_iupac("C1=C(C)CCCC1") == "1-methylcyclohexene"


def test_diene_with_substituent_cites_all_locants():
    assert smiles_to_iupac("C1=CC=C(C)CC1") == "1-methylcyclohexa-1,3-diene"


def test_halogen_substituent():
    assert smiles_to_iupac("C1=CC(Cl)CCC1") == "3-chlorocyclohexene"


def test_unsubstituted_cycloalkyne():
    # cyclooctyne: well-known unambiguous compound name, single ring
    # triple bond so its locant is always omitted (P-14.3.3), same as the
    # single-ene case above.
    assert smiles_to_iupac("C1#CCCCCCC1") == "cyclooctyne"


def test_cycloalkyne_with_substituent():
    # Same locant-tiebreak logic already verified for the ene case
    # (test_single_ene_locant_always_omitted_with_substituent) applied to
    # a triple bond: the yne stays at the omitted '1', methyl gets the
    # lowest achievable locant.
    assert smiles_to_iupac("CC1CCCCCC#C1") == "3-methylcyclooctyne"


def test_mixed_enyne_ring():
    # P-31.1.3.1's own worked example (Blue Book Chapter P-3): the double
    # bond is allocated locant '1' and the triple bond '4' -- lower
    # locants go to the double bond specifically once the combined
    # {1,4}/{1,2}... locant-set choice is tied. No stereo specified in the
    # input here, so no descriptor -- a specified E/Z double bond combined
    # with a triple bond is out of scope (see the E/Z tests below).
    assert smiles_to_iupac("C1=CCC#CCCCCCCCCCC1") == "cyclopentadec-1-en-4-yne"


def test_polycyclic_ring_triple_bond_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1#CC2(CCCC1)CCCCC2")


def test_exocyclic_double_bond_now_resolves():
    # A plain exocyclic '=CH2' substituent on an otherwise saturated,
    # unsubstituted ring is not ring unsaturation -- PubChem CID 14502,
    # matches its own IUPACName exactly.
    assert smiles_to_iupac("C1(=C)CCCCC1") == "methylidenecyclohexane"


def test_exocyclic_triple_bond_substituent_resolves():
    assert smiles_to_iupac("C1(C#CC)CCCCC1") == "(prop-1-yn-1-yl)cyclohexane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real, distinctly registered PubChem CIDs: 638079 ((Z)-cyclooctene),
        # 5463599 ((E)-cyclooctene), 5463153/5463190 (cyclononene Z/E),
        # 5365612/5364362 (cyclodecene Z/E). PubChem's own auto-generated
        # IUPACName for all six is just the bare parent name with no
        # descriptor at all (its namer has the same gap this module did) --
        # correctness is verified directly against P-91.2.2's own worked
        # example ('(Z)-cyclooctene', '(E)-cyclooctene', no locant) instead.
        ("C1CCC/C=C\\CC1", "(Z)-cyclooctene"),
        ("C1CCC/C=C/CC1", "(E)-cyclooctene"),
        ("C1CCCC/C=C\\CC1", "(Z)-cyclononene"),
        ("C1CCCC/C=C/CC1", "(E)-cyclononene"),
        ("C1CCCCC/C=C\\CC1", "(Z)-cyclodecene"),
        ("C1CCCCC/C=C/CC1", "(E)-cyclodecene"),
    ],
)
def test_ring_double_bond_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        # A 3-7-membered ring's C=C is not a genuine stereogenic unit
        # (P-91.2.2: 'Z' is the only physically realizable configuration
        # there) -- RDKit's own `Chem.FindPotentialStereo` already agrees
        # (no element reported), so a stereo marker in the input is
        # silently ignored, same as today.
        "C1CC=CC1",
        "C1CCC=CC1",
        "C1CCC/C=C\\C1",
    ],
)
def test_small_ring_stereo_marker_has_no_effect(smiles):
    unmarked = smiles.replace("/", "").replace("\\", "")
    assert smiles_to_iupac(smiles) == smiles_to_iupac(unmarked)


def test_ring_double_bond_stereo_with_unspecified_ring_stereocenter_raises():
    # The methyl-bearing ring carbon is a genuine (here left unspecified)
    # tetrahedral stereocenter of its own -- a specified double bond
    # alongside an unspecified stereocenter elsewhere is a partially
    # specified molecule, out of scope everywhere in this project (P-92/
    # P-93), not just here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCC/C=C\\CC1")


def test_ring_double_bond_stereo_multiple_bonds_raises():
    # P-91.2.2's multi-bond citation rule (locants become non-redundant
    # once 2+ double bonds are present) needs separate verification --
    # not attempted here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC/C=C\\C/C=C\\C1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A tetrahedral stereocenter on the ring itself (#776) -- a ring
        # double bond breaks the two-ring-walk-direction symmetry a
        # saturated ring has, so the sole substituent-bearing atom is a
        # genuine stereocenter even with no branch stereocenter of its
        # own. Real, distinctly registered PubChem CIDs: 25113480
        # (4-methyl), 96449479 (4-ethyl), 124507122 (4-propyl), 102028993
        # (4-propan-2-yl), 12463333 (4-tert-butyl).
        ("C[C@H]1CCC=CC1", "(4S)-4-methylcyclohexene"),
        ("CC[C@H]1CCC=CC1", "(4S)-4-ethylcyclohexene"),
        ("CCC[C@H]1CCC=CC1", "(4S)-4-propylcyclohexene"),
        ("CC(C)[C@H]1CCC=CC1", "(4S)-4-(propan-2-yl)cyclohexene"),
        ("CC(C)(C)[C@H]1CCC=CC1", "(4S)-4-tert-butylcyclohexene"),
        # The ring attachment atom *and* a stereocenter on its substituent
        # branch, both specified together -- PubChem CID 175915978
        # confirms the citation form combines the ring's own ordinary
        # on-ring "(4S)-" prefix with the branch's existing bracketed
        # descriptor. The three diastereomers/regiomers aren't themselves
        # separately registered CIDs, but exercise the same mechanism
        # with each independent R/S combination.
        ("CC[C@H](C)[C@H]1CCC=CC1", "(4S)-4-[(2S)-butan-2-yl]cyclohexene"),
        ("CC[C@H](C)[C@@H]1CCC=CC1", "(4R)-4-[(2S)-butan-2-yl]cyclohexene"),
        ("CC[C@@H](C)[C@H]1CCC=CC1", "(4S)-4-[(2R)-butan-2-yl]cyclohexene"),
        ("CC[C@@H](C)[C@@H]1CCC=CC1", "(4R)-4-[(2R)-butan-2-yl]cyclohexene"),
        ("C[C@@H](CC)[C@@H]1CCC=CC1", "(4R)-4-[(2S)-butan-2-yl]cyclohexene"),
    ],
)
def test_ring_tetrahedral_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_stereocenter_with_two_substituents_raises():
    # More than one substituent-bearing ring position alongside a
    # specified tetrahedral stereocenter is out of scope for this narrow
    # slice (see module docstring / #776) -- the ring's sole-substituent
    # restriction `ring_and_branch_stereo_display` shares with
    # `ring_branch_stereo_display` still applies.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@H](C)[C@H]1CC(C)C=CC1")
