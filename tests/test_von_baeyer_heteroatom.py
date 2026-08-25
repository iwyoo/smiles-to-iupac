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


def test_two_ring_heteroatoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CCN1O2")


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


def test_mixed_element_heteroatoms_raises():
    # P-23.2.1's 'a'-prefix ordering for two different heteroatom kinds
    # (Table 2.8 seniority) is out of scope -- this must not silently pick
    # an arbitrary order.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O1CC2CCC1N2")


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
