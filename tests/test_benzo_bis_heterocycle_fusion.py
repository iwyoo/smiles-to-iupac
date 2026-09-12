import pytest
from rdkit import Chem

from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._benzo_bis_heterocycle_fusion import (
    has_benzo_bis_heterocycle_fusion_name,
    name_benzo_bis_heterocycle_fusion,
)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 11106168 -- "BDT", a well-known organic-
        # semiconductor scaffold.
        ("C1=CSC2=CC3=C(C=C21)SC=C3", "benzo[1,2-b:4,5-b']dithiophene"),
        # PubChem CID 605456 -- the other real benzo-bond-pair isomer.
        ("C1=CC2=C(C=CS2)C3=C1C=CS3", "benzo[1,2-b:3,4-b']dithiophene"),
        # PubChem CID 57083408 -- confirms the mechanism isn't chalcogen-
        # specific.
        ("C1=COC2=CC3=C(C=C21)OC=C3", "benzo[1,2-b:4,5-b']difuran"),
        # PubChem CID 12274424
        ("C1=CC2=C(C=CO2)C3=C1C=CO3", "benzo[1,2-b:3,4-b']difuran"),
    ],
)
def test_benzo_bis_heterocycle_fusion_matches_pubchem(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    assert has_benzo_bis_heterocycle_fusion_name(mol)
    assert name_benzo_bis_heterocycle_fusion(mol) == expected


def test_mixed_parent_rings_out_of_scope():
    # One furan + one thiophene bridged by benzo - needs P-25.3.4.2.1's
    # seniority ordering between two *different* parents, a later step
    # (M5 step 2, #605).
    mol = Chem.MolFromSmiles("C1=CC2=C(C=CO2)C3=C1C=CS3")
    assert not has_benzo_bis_heterocycle_fusion_name(mol)


def test_parents_fused_to_each_other_out_of_scope():
    # dithieno[3,2-b:2',3'-d]thiophene (PubChem CID 137985, "DTT"): three
    # thiophenes in a linear chain, each pair directly ortho-fused to its
    # neighbor - not this module's "two identical parents bridged only by
    # one benzo" shape at all (no benzo ring here, and the middle ring is
    # itself a parent-type thiophene, not benzo).
    mol = Chem.MolFromSmiles("C1=CSC2=C1SC3=C2SC=C3")
    assert not has_benzo_bis_heterocycle_fusion_name(mol)


def test_substituent_out_of_scope():
    mol = Chem.MolFromSmiles("Cc1cc2occc2c(C)c1")
    assert not has_benzo_bis_heterocycle_fusion_name(mol)


def test_wrong_ring_count_out_of_scope():
    mol = Chem.MolFromSmiles("c1ccc2occc2c1")  # benzofuran, 2 rings only
    assert not has_benzo_bis_heterocycle_fusion_name(mol)
