import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 7-oxanorbornane: same bicyclo[2.2.1]heptane skeleton as
        # tests/test_bicyclic.py's norbornane case, with the one-atom
        # bridge (always numbered last, P-23.2.3) replaced by O -- the
        # heteroatom has no numbering choice here, it's pinned to locant 7.
        ("C1CC2CCC1O2", "7-oxabicyclo[2.2.1]heptane"),
        # same skeleton, sulfur instead of oxygen.
        ("C1CC2CCC1S2", "7-thiabicyclo[2.2.1]heptane"),
        # same skeleton, nitrogen instead of oxygen.
        ("C1CC2CCC1N2", "7-azabicyclo[2.2.1]heptane"),
        # quinuclidine: bicyclo[2.2.2]octane with a *bridgehead* nitrogen.
        # All three bridges are equal length, so (unlike the case above)
        # there's a real numbering choice between locant 1 and locant 4 for
        # the heteroatom-bearing bridgehead; the heteroatom-lowest-locant
        # tie-break must prefer 1.
        ("C1CN2CCC1CC2", "1-azabicyclo[2.2.2]octane"),
        # a halogen substituent coexists with the ring heteroatom; the
        # nondetachable 'oxa' prefix is cited right before the parent,
        # separated from the detachable 'chloro' prefix by a hyphen.
        ("C1C(Cl)C2CCC1O2", "2-chloro-7-oxabicyclo[2.2.1]heptane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_heteroatom_outside_ring_raises():
    # a plain hydrocarbon bicyclic with an exocyclic -OH substituent: the
    # single heteroatom isn't a *skeletal* ring atom, so this is a
    # different module's territory (deferred, see _alcohol.py).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CC2CCC1C2")


def test_unsupported_heteroatom_element_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CCC1P2")


def test_unsaturated_heteroatom_bicyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2C=CC1O2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 2,7-dioxabicyclo[2.2.1]heptane: same bicyclo[2.2.1]heptane
        # skeleton as the single-heteroatom cases above, with both bridge
        # atoms of the one-atom bridge position replaced -- cross-checked
        # against PubChem's own computed IUPACName for this exact SMILES
        # (C5H8O2; ConnectivitySMILES "C1CC2OCC1O2" independently matches
        # this project's own canonical SMILES for the same structure).
        ("O1CC2CCC1O2", "2,7-dioxabicyclo[2.2.1]heptane"),
        # same skeleton, nitrogen instead of oxygen -- PubChem-verified for
        # this exact SMILES ("2,7-diazabicyclo[2.2.1]heptane", C5H10N2).
        ("N1CC2CCC1N2", "2,7-diazabicyclo[2.2.1]heptane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom_multi(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_tetraaza_multiplying_prefix_elides_vowel():
    # P-16.3.3: 'tetra' + 'aza' -> 'tetraza' (the multiplying prefix's
    # terminal 'a' elides before the 'a'-initial heteroatom prefix), not
    # 'tetraaza' -- PubChem-verified for this exact SMILES
    # ("12-methyl-1,4,7,10-tetrazabicyclo[8.3.2]pentadecane").
    assert (
        smiles_to_iupac("CC1CN2CCNCCNCCN(CC2)C1")
        == "12-methyl-1,4,7,10-tetrazabicyclo[8.3.2]pentadecane"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # bicyclo[2.2.1]heptane, O in a 2-atom bridge and N in the 1-atom
        # bridge: PubChem-verified for this exact SMILES
        # ("2-oxa-7-azabicyclo[2.2.1]heptane", CID 153796098).
        ("O1CC2CCC1N2", "2-oxa-7-azabicyclo[2.2.1]heptane"),
        # same skeleton, N at a *bridgehead* instead of in a bridge --
        # PubChem-verified (CID 23462610): the bridgehead locant (1) is
        # still cited after 'oxa' even though it's numerically lower than
        # oxa's own locant (7), confirming citation order really is by
        # Table 2.8 seniority and not by locant value.
        ("C1CC2CCN1O2", "7-oxa-1-azabicyclo[2.2.1]heptane"),
        # bicyclo[3.2.1]octane, O and N each in a different bridge --
        # PubChem-verified (CID 12069231): 'oxa' is cited before 'aza'
        # even though its own locant (8) is numerically higher than aza's
        # (3), confirming Table 2.8 citation order is by seniority, not by
        # locant value.
        ("C1CC2CNCC1O2", "8-oxa-3-azabicyclo[3.2.1]octane"),
        # same skeleton, O and S -- PubChem-verified (CID 118210548).
        ("C1CC2CSCC1O2", "8-oxa-3-thiabicyclo[3.2.1]octane"),
        # same skeleton, S and N -- PubChem-verified (CID 102218589):
        # confirms S is senior to N too (S > N in Table 2.8), not just
        # O > N.
        ("C1CC2CSCC1N2", "3-thia-8-azabicyclo[3.2.1]octane"),
        # bicyclo[2.2.1]heptane built with O and N in mirror-image
        # positions (both bridging directly off the same bridgehead):
        # isolates the tie-break rule itself, since the overall locant
        # *set* {2,6} is symmetric either way one of the two 2-atom
        # bridges gets numbered first. PubChem's own computed name
        # (CID 154021213) gives O the lower locant, confirming the more
        # senior element -- not the lower atom index -- wins the tie.
        ("C1NC2CC1CO2", "2-oxa-6-azabicyclo[2.2.1]heptane"),
        # same mirror-image construction, S and N instead of O and N --
        # confirms the same tie-break for the S > N pair (CID 137400758).
        ("C1NC2CC1CS2", "2-thia-6-azabicyclo[2.2.1]heptane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom_mixed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_three_mixed_element_heteroatoms_raises():
    # Table 2.8 seniority ordering here only covers exactly one heteroatom
    # of each of two different elements -- a third heteroatom (even a
    # repeat of one already-present element) must still be rejected rather
    # than silently picking an arbitrary order among three.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2NNC1CO2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 2-oxaadamantane: adamantane (tricyclo[3.3.1.1^3,7]decane, see
        # tests/test_polycyclic.py) with a non-bridgehead -CH2- replaced by
        # -O-. Cross-checked against PubChem's own computed IUPACName for
        # this exact SMILES ("2-oxatricyclo[3.3.1.1{3,7}]decane",
        # C9H14O) -- unlike the naphthalene-bridge/phane cases, PubChem's
        # name generator does implement this shape, so it's usable here.
        ("O1C2CC3CC1CC(C2)C3", "2-oxatricyclo[3.3.1.1^3,7]decane"),
        # a halogen substituent coexists with the ring heteroatom, same as
        # the bicyclic case above. PubChem-verified for this exact SMILES
        # ("5-chloro-2-oxatricyclo[3.3.1.1{3,7}]decane", C9H13ClO).
        ("O1C2CC3CC1CC(Cl)(C2)C3", "5-chloro-2-oxatricyclo[3.3.1.1^3,7]decane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom_tricyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_heteroatoms_tricyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O1C2CC3CC1CC(C2)N3")


def test_unsaturated_heteroatom_tricyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O1C2CC3CC1C=C(C2)C3")


def test_tricyclic_hydrocarbon_itself_is_unaffected():
    assert smiles_to_iupac("C1C2CC3CC1CC(C2)C3") == "tricyclo[3.3.1.1^3,7]decane"
