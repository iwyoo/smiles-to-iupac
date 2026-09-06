import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-14.3.4.2(d): the locant '1' is omitted for unsubstituted two- and
        # three-carbon alkenes/alkynes.
        ("C=C", "ethene"),
        ("C#C", "acetylene"),  # P-31.1.2.1: 'acetylene' is the retained PIN, not 'ethyne'.
        ("CC=C", "propene"),
        ("CC#C", "propyne"),
        # Four carbons and up: the locant is always essential, since more
        # than one position is possible.
        ("CC=CC", "but-2-ene"),
        ("C=CCC", "but-1-ene"),
        ("CC#CC", "but-2-yne"),
        ("C#CCC", "but-1-yne"),
        ("C#CCCC", "pent-1-yne"),
        ("CCCC=C", "pent-1-ene"),
        ("CCC#CC", "pent-2-yne"),
        # P-14.3.3 / P-14.3.4.2(d): once a substituent is present, even a
        # three-carbon chain must cite the otherwise-omittable locant.
        ("C=C(C)C", "2-methylprop-1-ene"),
        ("C=C(C)CC", "2-methylbut-1-ene"),
        ("CC(C)=CC", "2-methylbut-2-ene"),
        ("CC(C)C=C", "3-methylbut-1-ene"),
        # P-14.4(e): lowest locant for the double bond outranks lowest locant
        # for substituents when choosing numbering direction. Numbering left
        # to right gives the double bond locant 2 and the methyl locant 6;
        # numbering right to left would give the methyl locant 2 but pushes
        # the double bond to locant 5. The double bond wins the tie-break.
        ("CC=CCCC(C)C", "6-methylhept-2-ene"),
        # A branched ("compound") substituent (P-29.4) on an unsaturated
        # chain: CCCC(C(C)C)CC=C is a 7-carbon chain (the unique longest
        # chain: the isopropyl branch's own arm is only 2 atoms long,
        # shorter than the 3-atom continuation on the other side, so it can't
        # tie for principal chain) carrying an isopropyl branch (root
        # forking into two methyls -> "propan-2-yl", P-29.3.2.2) plus a
        # double bond at the far end. Numbering from the double-bond end
        # (P-14.4e) gives the double bond locant 1 and the branch locant 4.
        # This is only reachable now because substituents on an unsaturated
        # chain are named via `_substituents.py`'s `name_branch`, the same as
        # for alkanes/cycloalkanes, instead of being rejected outright.
        ("CCCC(C(C)C)CC=C", "4-(propan-2-yl)hept-1-ene"),
        # P-14.3.4.2(b): a two-carbon chain has only one possible bond
        # position, so the -ene/-yne locant is omittable regardless of how
        # many substituents are cited -- with one substituent, its own
        # locant is also omittable (no ambiguity); with two or more, their
        # locants are still needed to distinguish isomers like "1,2-" from
        # "1,1-".
        ("ClC=C", "chloroethene"),
        ("FC#C", "fluoroethyne"),
        ("ClC=CCl", "1,2-dichloroethene"),
        ("FC=CF", "1,2-difluoroethene"),
        ("FC=CCl", "1-chloro-2-fluoroethene"),
        ("BrC#CBr", "1,2-dibromoethyne"),
    ],
)
def test_smiles_to_iupac_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-31.1.1.2: multiplying prefix + comma-separated locants before
        # the (elided) suffix.
        ("C=CC=C", "buta-1,3-diene"),
        ("C=CCC=C", "penta-1,4-diene"),
        ("C=CC=CC", "penta-1,3-diene"),
        ("C#CC#C", "buta-1,3-diyne"),
        # Cross-checked against the Blue Book's own P-31.1.1.2 example,
        # 'nona-1,3,5,7-tetraene (PIN)'.
        ("C=CC=CC=CC=CC", "nona-1,3,5,7-tetraene"),
        ("C=CC=CC=C", "hexa-1,3,5-triene"),
        # P-31.1.1.1 / P-31.1.2.2.1: a double bond and a triple bond combine
        # into a single 'en...yne' suffix, eliding 'ene' to 'en' before the
        # unprefixed 'yne'.
        ("C=CC#C", "but-1-en-3-yne"),
        # 2-methylbuta-1,3-diene (isoprene): cross-checked against the Blue
        # Book's own P-31.1.2.1 note that this is the PIN for isoprene.
        ("C=C(C)C=C", "2-methylbuta-1,3-diene"),
        # A branched ("compound") substituent (P-29.4) together with two
        # double bonds: CCCC(C(C)C)CC=C ('hept-1-ene' single-bond case
        # above) extended with a second double bond at the chain's other
        # end. The 7-carbon chain is still the unique longest chain (the
        # isopropyl branch's own arm is too short to compete either
        # way), and it is symmetric enough that both numbering directions
        # tie on the double-bond locant set {1,6} and on the branch's own
        # locant (4), giving one unambiguous name.
        ("C=CCC(C(C)C)CC=C", "4-(propan-2-yl)hepta-1,6-diene"),
        # P-14.4(e)(ii) / P-31.1.1.1 / P-44.4.1.10.1: lowest locants go to
        # the full set of multiple bonds first ({2,4} either numbering
        # direction gives the same set), then to the double bond
        # specifically when a choice remains. Numbering left-to-right gives
        # the double bond locant 2 (and the triple bond locant 4); the
        # reverse direction would tie on the set {2,4} but swap them, giving
        # the double bond locant 4 instead - left-to-right wins.
        ("CC=CC#CC", "hex-2-en-4-yne"),
        # A multiplying-prefixed 'yne' (diyne/triyne) elides 'ene's
        # trailing 'e' exactly like an unprefixed 'yne' does -- the
        # elision is triggered by the underlying 'yne' word itself, not
        # by whether the final prefixed word happens to start with a
        # vowel. Found via real-data testing (PubChem-verified):
        # 'deca-1,2,3-trien-5,7,9-triyne', not '...triene-5,7,9-triyne'.
        ("C#CC#CC#CC=C=C=C", "deca-1,2,3-trien-5,7,9-triyne"),
    ],
)
def test_smiles_to_iupac_multiply_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_branch_not_on_principal_chain_raises():
    # The double bond sits on a short branch whose own arm can never be as
    # long as the two 6-carbon arms of the main chain, so no longest chain
    # contains it; expressing it would need an alkenyl substituent prefix
    # (P-29.2/P-32.1), which is out of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCC(C=C)CCCCCC")


def test_unsaturated_not_all_multiple_bonds_on_one_chain_raises():
    # Same short branch as above, but the main chain now also carries its
    # own double bond. That bond alone would be fine, but the branch's
    # double bond still can't be expressed without an alkenyl substituent
    # prefix, so the whole molecule stays out of scope (P-44.4.1.1: only a
    # chain carrying *every* multiple bond is eligible as principal chain).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCCCC(C=C)CCCCCC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single, specified C=C double bond (P-91.3/P-93): the primary
        # source's own worked example is "(2Z)-but-2-ene (PIN)" -- P-91.3
        # requires a locant before every stereodescriptor in an acyclic
        # name, even with only one double bond in the molecule (the bare,
        # locant-less "(Z)-" form is reserved for specific ring systems
        # per P-91.2.2, not chains).
        ("C/C=C/C", "(2E)-but-2-ene"),
        ("C/C=C\\C", "(2Z)-but-2-ene"),
        # Combined with a substituent prefix. Note this project's own
        # SMILES for the "E" case canonicalizes to PubChem's CID 5364761
        # structure (verified via Chem.CanonSmiles), not CID 5463022 as an
        # earlier, unverified guess in the task file's own background
        # section assumed -- CIP priority (Cl outranking CH3) flips which
        # double-bond drawing is E vs Z relative to naive left-right
        # geometry, so this is exactly the kind of case this module
        # deliberately delegates to RDKit's rdCIPLabeler rather than
        # guessing.
        ("C/C(Cl)=C\\C", "(2E)-2-chlorobut-2-ene"),
        ("C/C(Cl)=C/C", "(2Z)-2-chlorobut-2-ene"),
    ],
)
def test_smiles_to_iupac_ez_double_bond(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_double_bond_geometry_unaffected():
    # no `/`/`\` at all: not a new rejection case -- named exactly as
    # before (no E/Z prefix), matching this project's prior behavior.
    assert smiles_to_iupac("CC=CC") == "but-2-ene"


def test_non_stereogenic_double_bond_unaffected():
    # one alkene carbon has two identical substituents (both H): no real
    # geometric isomerism exists, so there's nothing to prefix.
    assert smiles_to_iupac("C=C(C)C") == "2-methylprop-1-ene"


def test_specified_double_bond_with_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C/C=C/CC#C")


def test_specified_double_bond_with_second_double_bond_raises():
    # The second double bond (terminal =CH2) is non-stereogenic, so it's
    # left unspecified -- a specified bond alongside an unspecified/
    # non-stereogenic one is still out of scope (every double bond must be
    # specified for the multi-E/Z case below to apply).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C/C=C/C=C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem PUG REST (queried via POST -- GET 400s on a SMILES
        # containing '/'): CID 638071, auto-generated name matches exactly.
        ("C/C=C/C=C/C", "(2E,4E)-hexa-2,4-diene"),
        # CID 5326156.
        ("C/C=C\\C=C/C", "(2Z,4Z)-hexa-2,4-diene"),
        # CID 643786.
        ("C/C=C\\C=C\\C", "(2Z,4E)-hexa-2,4-diene"),
        # CID 5368766 -- confirms the pattern extends to three double bonds.
        ("C/C=C/C=C/C=C/C", "(2E,4E,6E)-octa-2,4,6-triene"),
    ],
)
def test_multi_ez_double_bond(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_partially_specified_diene_raises():
    # Both double bonds are stereogenic, but only one has a slash marker
    # -- a partially-specified molecule is out of scope, same policy as
    # `specified_stereocenters`'s R/S handling.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C/C=C/C=CC")
