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
        # chain: the isopropyl-like fork's own arm is only 2 atoms long,
        # shorter than the 3-atom continuation on the other side, so it can't
        # tie for principal chain) carrying an isopropyl-like branch
        # (root forking into two methyls -> "1-methylethyl", P-46) plus a
        # double bond at the far end. Numbering from the double-bond end
        # (P-14.4e) gives the double bond locant 1 and the branch locant 4.
        # This is only reachable now because substituents on an unsaturated
        # chain are named via `_substituents.py`'s `name_branch`, the same as
        # for alkanes/cycloalkanes, instead of being rejected outright.
        ("CCCC(C(C)C)CC=C", "4-(1-methylethyl)hept-1-ene"),
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
        # isopropyl-like branch's own arm is too short to compete either
        # way), and it is symmetric enough that both numbering directions
        # tie on the double-bond locant set {1,6} and on the branch's own
        # locant (4), giving one unambiguous name.
        ("C=CCC(C(C)C)CC=C", "4-(1-methylethyl)hepta-1,6-diene"),
        # P-14.4(e)(ii) / P-31.1.1.1 / P-44.4.1.10.1: lowest locants go to
        # the full set of multiple bonds first ({2,4} either numbering
        # direction gives the same set), then to the double bond
        # specifically when a choice remains. Numbering left-to-right gives
        # the double bond locant 2 (and the triple bond locant 4); the
        # reverse direction would tie on the set {2,4} but swap them, giving
        # the double bond locant 4 instead - left-to-right wins.
        ("CC=CC#CC", "hex-2-en-4-yne"),
    ],
)
def test_smiles_to_iupac_multiply_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C1=CCCCC1",  # cyclohexene: unsaturation in a ring, out of scope.
    ],
)
def test_out_of_scope_unsaturation_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
