import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Norborneol-type bicyclic alcohol (PubChem-confirmed
        # 'bicyclo[2.2.1]heptan-2-ol', CID 19809's own structure).
        ("OC1CC2CCC1C2", "bicyclo[2.2.1]heptan-2-ol"),
        # Same skeleton as test_halogens.py's 'ClC1CC2CCC1CC2' ->
        # '2-chlorobicyclo[2.2.2]octane', with the halogen swapped for
        # -OH (PubChem-confirmed 'bicyclo[2.2.2]octan-2-ol').
        ("OC1CC2CCC1CC2", "bicyclo[2.2.2]octan-2-ol"),
        # A fused (all-bridge-lengths-nonzero-except-one) bicyclic
        # decalinol.
        ("OC1CCC2CCCCC2C1", "bicyclo[4.4.0]decan-3-ol"),
        # A bridgehead hydroxyl on a fused bicyclooctane.
        ("OC12CCCC1CCC2", "bicyclo[3.3.0]octan-1-ol"),
        # A hydroxyl coexisting with a plain alkyl substituent, both on
        # the ring.
        ("OC1(C)CC2CCC1C2", "2-methylbicyclo[2.2.1]heptan-2-ol"),
        # A hydroxyl coexisting with a halogen substituent, both on the
        # ring -- the -OH locant still wins the tie-break ahead of the
        # halogen's (mirrors `_von_baeyer_heteroatom.py`'s own
        # heteroatom-before-substituent rank).
        ("ClC1CC2CCC1C(O)C2", "6-chlorobicyclo[2.2.2]octan-2-ol"),
        # 1-adamantanol (PubChem CID 12178's own structure) -- this
        # project's adamantane already uses the systematic
        # 'tricyclo[3.3.1.1^3,7]decane' name rather than the retained
        # 'adamantane' one (see test_tricyclic.py), so this module's
        # ring_count>=3 path follows the same systematic convention:
        # 'tricyclo[3.3.1.1^3,7]decan-1-ol', not 'adamantan-1-ol'.
        ("OC12CC3CC(CC(C3)C1)C2", "tricyclo[3.3.1.1^3,7]decan-1-ol"),
        # 2-adamantanol (PubChem CID 68159's own structure), same
        # systematic-name divergence as 1-adamantanol above.
        ("OC1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.1^3,7]decan-2-ol"),
    ],
)
def test_von_baeyer_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_hydroxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CC2CCC1C2O")


def test_hydroxyl_on_substituent_branch_raises():
    # The -OH sits on a chain hanging off the bicyclic ring, not on the
    # ring skeleton itself -- a different (still unimplemented) shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC1CC2CCC1C2")


def test_unsaturated_von_baeyer_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CC2C=CC1C2")


def test_disjoint_rings_alcohol_still_raises():
    # Two rings joined only by a single bond (not fused, bridged, or
    # spiro) never matches `find_bicyclic_core`/`find_polycyclic_core`'s
    # own core requirement, so this regression check confirms the new
    # von Baeyer routing doesn't misfire on this still-unsupported shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1(CCCC1)C1CCCC1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed real structures, each PubChem's own IUPACName
        # property matches exactly (#1079, M5 step 2, mirroring #1078's
        # identical ketone extension): the bicyclic suffix engine never
        # learned to cite a specified stereocenter for an alcohol before
        # this, unconditionally rejecting every one of these.
        ("C[C@@]12CC[C@@H](C1(C)C)C[C@H]2O", "(1R,2R,4R)-1,7,7-trimethylbicyclo[2.2.1]heptan-2-ol"),  # isoborneol, CID 6321405
        ("C[C@]12CC[C@H](C1(C)C)C[C@H]2O", "(1S,2R,4S)-1,7,7-trimethylbicyclo[2.2.1]heptan-2-ol"),  # (1S,2R,4S)-borneol, CID 1201518
        ("C[C@@]12CC[C@@H](C1(C)C)C[C@@H]2O", "(1R,2S,4R)-1,7,7-trimethylbicyclo[2.2.1]heptan-2-ol"),  # (1R,2S,4R)-borneol, CID 6552009
    ],
)
def test_von_baeyer_alcohol_specified_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_myrtenol_exocyclic_hydroxyl_still_raises():
    # Superficially similar to the ring-unsaturation shape this module
    # handles, but the hydroxyl sits on an exocyclic CH2 off the ring, not
    # on the ring itself -- ruled out during M6 step 2's own investigation
    # (#1082) by running it through smiles_to_iupac() directly rather than
    # trusting its PubChem name's "...enyl)methanol" phrasing. Still
    # raises today (ring unsaturation, unrelated to this step), unaffected
    # by this stereocenter-citation change.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1(C2CC=C(C1C2)CO)C")
