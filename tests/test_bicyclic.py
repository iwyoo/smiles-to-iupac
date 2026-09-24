import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # norbornane: bridges 2, 2, 1 -> bicyclo[2.2.1]heptane; bridgeheads at
        # C1/C4 and the one-carbon bridge at C7 (cross-checked against PubChem
        # CID 9233 / Wikipedia "Norbornane").
        ("C1CC2CCC1C2", "bicyclo[2.2.1]heptane"),
        # bicyclo[2.2.2]octane: three equal-length (2, 2, 2) bridges; RDKit's
        # SSSR reports 3 rings for this SMILES (redundant, non-minimal), not
        # 2 (see _bicyclic.py's find_bicyclic_core docstring), so this is the
        # key regression case for not gating detection on NumRings().
        ("C1CC2CCC1CC2", "bicyclo[2.2.2]octane"),
        # A bare, unsubstituted two-fused-six-membered-ring skeleton
        # (bridges 4, 4, 0) is decahydronaphthalene, not a von Baeyer
        # name at all -- see test_dihydro_aromatic.py's own coverage,
        # routed ahead of this module in core.py.
        # bicyclo[3.3.0]octane: two fused five-membered rings sharing a bond
        # -> bridges 3, 3, 0.
        ("C12CCCC1CCC2", "bicyclo[3.3.0]octane"),
        # methyl on the ring atom next to a bridgehead, in one of the two
        # 2-atom bridges of norbornane (P-23.2.3: numbering starts at the
        # bridgehead adjacent to it and proceeds along that bridge) -> C2.
        ("C1C(C)C2CCC1C2", "2-methylbicyclo[2.2.1]heptane"),
    ],
)
def test_smiles_to_iupac_bicyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_pentagonal_prism_is_not_bicyclic():
    # pentaprismane: a genuinely hexacyclic system (cyclomatic number 6, ten
    # branch atoms of degree 3) -- not mistaken for bicyclic by
    # find_bicyclic_core, and correctly resolved via _polycyclic.py's
    # hexacyclic support instead (see tests/test_hexacyclic.py).
    assert smiles_to_iupac("C12C3C4C1C1C2C2C3C4C12") == "hexacyclo[4.4.0.0^2,5.0^3,9.0^4,8.0^7,10]decane"


def test_two_separate_rings_routes_to_ring_assembly_module():
    # Two disjoint cyclohexane rings joined by a single bond share no
    # atom at all, so this is not a bicyclic system -- it's
    # `_ring_assembly.py`'s own P-28.2.1 ring-assembly shape instead.
    assert smiles_to_iupac("C1CCCCC1C1CCCCC1") == "1,1'-bi(cyclohexane)"




def test_unsaturated_bicyclic_name():
    assert smiles_to_iupac("C1CC2CC=C1C2") == "bicyclo[2.2.1]hept-1-ene"


def test_monospiro_still_resolves_via_spiro_module():
    assert smiles_to_iupac("C1CCCC12CCCCC2") == "spiro[4.5]decane"
