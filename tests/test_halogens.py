import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simple monosubstituted chain (P-35.2.1).
        ("CCCF", "1-fluoropropane"),
        # P-14.3.4.2(b): a monosubstituted, homogeneous two-carbon chain has
        # only one possible structure regardless of numbering direction, so
        # the locant is omittable, analogous to 'ethanol (PIN)'.
        ("CCCl", "chloroethane"),
        # P-14.3.4.2(a): every substituent's locant on a mononuclear (single
        # skeletal atom) parent hydride is always '1' and is never cited, no
        # matter how many substituents there are, e.g. 'chloromethane (PIN)'.
        ("C(Cl)(Cl)(Cl)Cl", "tetrachloromethane"),
        # Multiple identical halogens need a multiplying prefix.
        ("FC(F)CC", "1,1-difluoropropane"),
        # Different halogens alphabetize by their own name (P-14.5.2): bromo,
        # chloro, fluoro, iodo, in that order - here all four at once.
        ("BrCC(Cl)C(F)CI", "1-bromo-2-chloro-3-fluoro-4-iodobutane"),
        # A halogen substituent alongside an alkyl substituent: alphabetical
        # order puts 'chloro' before 'methyl' (P-14.5.2), so it wins the
        # lower locant on this otherwise-symmetric tie.
        ("CC(Cl)CC(C)CC", "2-chloro-4-methylhexane"),
        # A halogen on a ring, reusing _cyclic.py's existing machinery; the
        # only substituent on an otherwise unsubstituted ring, so the locant
        # is omitted (P-14.3.3).
        ("ClC1CCCCC1", "chlorocyclohexane"),
        # A halogen alongside an alkyl substituent on a ring.
        ("CC1CCCCC1Cl", "1-chloro-2-methylcyclohexane"),
        # A halogen on a monospiro ring (P-24.2.1's numbering machinery,
        # reused unchanged from _cyclic.py/_spiro.py).
        ("ClC1CCC2(C1)CCCCC2", "2-chlorospiro[4.5]decane"),
        # A halogen on a bicyclic (von Baeyer) ring, added after merging with
        # _bicyclic.py: same _substituents_for_ring reuse as the spiro case.
        ("ClC1CC2CCC1CC2", "2-chlorobicyclo[2.2.2]octane"),
        # A chloro-substituted alkene: halogen substituent combined with
        # unsaturation.
        ("C=C(Cl)CC", "2-chlorobut-1-ene"),
        # Cross-checked directly against the Blue Book's own P-31.1.2.1
        # example: 'fluoroethyne (PIN)' for fluoroacetylene, HC#C-F.
        ("C#CF", "fluoroethyne"),
        # A compound (branched) substituent whose own chain is a single,
        # mononuclear carbon bearing the halogen (P-14.3.4.2(a) again, this
        # time for a substituent's own internal numbering): 'chloromethyl',
        # not '1-chloromethyl'.
        ("CCCC(CCl)CCC", "4-(chloromethyl)heptane"),
        # A compound substituent whose own two-carbon chain genuinely has two
        # different possible chloro positions (unlike the mononuclear case
        # above), so the locant is essential: '1-chloroethyl'.
        ("CCCC(C(Cl)C)CCC", "4-(1-chloroethyl)heptane"),
        # Regression test for the carbon-only chain search fix: naively
        # running the longest-chain BFS over the full atom graph (including
        # the terminal Cl) finds a 5-atom path C-C-C-C-Cl and would wrongly
        # treat Cl as a 5th skeletal atom, picking the wrong (nonexistent)
        # "pentane-like" principal chain. The true carbon skeleton is only 4
        # atoms long (an isobutyl-shaped chain: (CH3)2CH-CH2-CH2-Cl), giving
        # the correct name below (cross-checked as isoamyl chloride's PIN).
        ("CC(C)CCCl", "1-chloro-3-methylbutane"),
    ],
)
def test_halogen_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_lone_halogen_atom_raises():
    # No carbon atom at all: no hydrocarbon parent hydride to substitute.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("FCl")


def test_non_halogen_heteroatom_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCN")
