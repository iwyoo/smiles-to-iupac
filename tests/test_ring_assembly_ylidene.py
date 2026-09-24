import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-28.2.2), spanning three ring sizes.
        ("C1CCC(=C2CCCC2)C1", "1,1'-bi(cyclopentylidene)"),  # CID 549125
        ("C1CCC(=C2CCCCC2)CC1", "1,1'-bi(cyclohexylidene)"),  # CID 138159
        ("C1CC1=C1CC1", "1,1'-bi(cyclopropylidene)"),
        ("C1CCC1=C1CCC1", "1,1'-bi(cyclobutylidene)"),
    ],
)
def test_ring_assembly_ylidene_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituent():
    assert smiles_to_iupac("C1CCC(Cl)C1=C1CCCC1") == "2-chloro-1,1'-bi(cyclopentylidene)"


def test_biphenyl_still_routes_to_single_bond_module():
    assert smiles_to_iupac("c1ccccc1-c1ccccc1") == "1,1'-biphenyl"


def test_different_ring_sizes_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC1=C1CCCC1")


def test_single_ring_exocyclic_double_bond_unaffected():
    # Not a ring assembly at all -- one ring, one open-chain =CH-CH3 tail.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1=CC")


def test_bicyclic_ring_assembly_ylidene_matches_pin_worked_example():
    # The primary source's own literal PIN worked example
    # (tmp/bluebook/P2.txt ~7826): two bicyclo[2.2.1]heptane rings, each
    # attached at position 2, joined by a C=C double bond.
    assert (
        smiles_to_iupac("C1CC2CC1CC2=C1CC2CCC1C2")
        == "2,2'-bi(bicyclo[2.2.1]heptanylidene)"
    )


def test_bicyclic_ring_assembly_ylidene_bridge_carbon_attachment():
    # Attached at the one-carbon bridge (locant 7) instead of locant 2.
    assert (
        smiles_to_iupac("C1CC2CCC1C2=C1C2CCC1CC2")
        == "7,7'-bi(bicyclo[2.2.1]heptanylidene)"
    )


def test_bicyclic_ring_assembly_ylidene_mismatched_bridge_lengths_raises():
    # One side is bicyclo[2.2.1]heptane, the other a different bridge-length
    # bicyclic skeleton -- not "identical cyclic systems" (P-28.1).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2C(=C3CC4CCC3C4)CCC2C1")


def test_bicyclic_ring_assembly_ylidene_halogen_substituent_raises():
    # A halogen substituent on a bicyclic side is deferred (see this
    # module's own docstring for why) -- confirms it's rejected outright,
    # not silently mis-numbered.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CC2CCC1C2=C1CC2CCC1C2")
