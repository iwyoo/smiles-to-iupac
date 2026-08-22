import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


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
        # decahydronaphthalene (decalin): two fused six-membered rings
        # sharing a bond -> bridges 4, 4, 0.
        ("C1CCC2CCCCC2C1", "bicyclo[4.4.0]decane"),
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


def test_tricyclic_raises():
    # cubane: a genuinely pentacyclic system (cyclomatic number 5, eight
    # branch atoms of degree 3) -- still out of scope even after
    # tricyclic/tetracyclic support (see test_tricyclic.py,
    # tests/test_tetracyclic.py), so this remains unsupported polycyclic
    # territory. (Previously this test used a tetracyclic SMILES, but that
    # topology is now correctly supported by _tetracyclic.py, so it no
    # longer demonstrates "out of scope" -- see _tetracyclic.py.)
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C1C2C3C41")


def test_two_separate_rings_still_raises():
    # two cyclohexane rings joined by a single bond: two rings sharing no
    # atom at all, so not a bicyclic system (P-23.2.2 requires "two or more
    # atoms in common").
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C1CCCCC1")


def test_unsaturated_bicyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CC=C1C2")


def test_heteroatom_bicyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CCC1O2")


def test_monospiro_still_resolves_via_spiro_module():
    assert smiles_to_iupac("C1CCCC12CCCCC2") == "spiro[4.5]decane"
