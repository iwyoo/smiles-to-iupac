import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # adamantane (flagship case): four bridgeheads, each pair connected by
        # a one-carbon bridge (tetrahedral K4 arrangement). Cross-checked
        # against the NIST WebBook / PubChem name for CAS 281-23-2, and
        # independently re-derived by hand from P-23.2.2/P-23.2.3/P-23.2.5.1/
        # P-23.2.5.2 (see _tricyclic.py's module docstring): with main
        # bridgeheads numbered 1 and 5, the other two bridgeheads land at the
        # midpoints of the two 3-atom main-ring segments, positions 3 and 7.
        ("C1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.1^3,7]decane"),
        # twistane: built from scratch (not copied from a SMILES database) to
        # match the descriptor tricyclo[4.4.0.0^3,8]decane reported by
        # Wikipedia/ChemSpider (CAS 253-14-5) -- main bridgeheads directly
        # bonded (c=0) and the secondary bridgeheads (at positions 3 and 8,
        # one on each 4-atom main-ring segment) also directly bonded (d=0).
        ("C1CC2CC3CCC2CC13", "tricyclo[4.4.0.0^3,8]decane"),
        # 1-methyladamantane: methyl on a main bridgehead. All four
        # bridgeheads are equivalent in unsubstituted adamantane, so the
        # lowest-locant tie-break (P-14.4/P-45.2) picks the numbering that
        # puts the substituted one at position 1.
        ("CC12CC3CC(CC(C3)C1)C2", "1-methyltricyclo[3.3.1.1^3,7]decane"),
        # Same, with a halogen substituent -- proves `halogens` threads
        # through find_tricyclic_core/name_tricycloalkane the same way it
        # does through _bicyclic.py (see this repo's history for the exact
        # bug this guards against: halogens dropped by a bicyclic/tricyclic
        # merge that didn't thread the map through).
        ("ClC12CC3CC(CC(C3)C1)C2", "1-chlorotricyclo[3.3.1.1^3,7]decane"),
    ],
)
def test_smiles_to_iupac_tricyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bicyclic_is_not_tricyclic():
    # bicyclo[2.2.2]octane: genuinely bicyclic (cyclomatic number 2), still
    # resolved by _bicyclic.py unaffected by tricyclic support.
    assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"


def test_two_separate_rings_still_raises():
    # two cyclohexane rings joined by a single bond: share no atom at all,
    # not a polycyclic ring system by any of these modules' definitions.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C1CCCCC1")


def test_tetracyclic_raises():
    # cyclomatic number 4 (six branch atoms in the leaf-stripped core, not
    # four): genuinely beyond this module's scope (P-23.2.6, not P-23.2.5).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1C2CC3CC1C1C(C2)C31")


def test_propellane_like_degree_four_raises():
    # [1.1.1]propellane: only two branch atoms, both of degree 4 (directly
    # bonded to each other in addition to the three one-carbon bridges), not
    # four branch atoms of degree 3 -- the exact "secondary bridge
    # reconnects to a main bridgehead" topology this module deliberately
    # excludes (see _tricyclic.py's module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1C23CC12C3")


def test_secondary_bridgeheads_on_same_main_bridge_raises():
    # Built directly from the Blue Book's own P-23.2.5.2 worked example,
    # tricyclo[4.2.2.2^2,5]dodecane: both secondary-bridge attachment points
    # (locants 2 and 5) lie on the *same* 4-atom main-ring segment, so the
    # reduced graph on the four branch atoms has a doubled edge between the
    # two main bridgeheads (via the two length-2 bridges) and another
    # doubled edge between the two secondary bridgeheads (via the segment
    # interior and the secondary bridge itself), rather than the six
    # distinct single bridges (K4) this module's detection requires. A real,
    # named von Baeyer topology that is out of scope for this
    # implementation -- see _tricyclic.py's module docstring -- and
    # correctly rejected rather than mis-named.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CCC1C1CCC2CC1")
