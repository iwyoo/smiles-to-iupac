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
        # PubChem CID 605316 -- a real, structurally distinct isomer from
        # CID 11106168 above (different InChIKey - a different regio-
        # chemistry of the two thiophene S atoms, not just a citation-
        # order variant). Found during #607's PubChem validation pass:
        # an earlier version of this module's citation-order logic sorted
        # by the *attached* benzo ring's own numbering instead of the
        # *parent* thiophene's, which happened to give the same wrong
        # name ("...4,5-b'...") for this isomer as for CID 11106168 -
        # silently conflating two different real compounds into one name.
        ("C1=CSC2=CC3=C(C=CS3)C=C21", "benzo[1,2-b:5,4-b']dithiophene"),
    ],
)
def test_benzo_bis_heterocycle_fusion_matches_pubchem(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    assert has_benzo_bis_heterocycle_fusion_name(mol)
    assert name_benzo_bis_heterocycle_fusion(mol) == expected


def test_mixed_parent_rings_out_of_scope():
    # One furan + one thiophene bridged by benzo - #605 found this isn't
    # even a real multiparent case (the senior parent, furan, would just
    # win outright and use a retained bicyclic base instead - see
    # `thieno[3,2-f][1]benzofuran`, PubChem CID 69035235), so this
    # structure is correctly rejected here rather than misnamed.
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


def test_naphtho_bridge_out_of_scope():
    # naphtho[1,2-b:5,6-b']dithiophene (PubChem CID 58434773), found
    # during #607's validation pass: a real compound in the same spirit
    # (two identical thiophene parents bridged by one carbocycle) but with
    # a naphtho (4-ring total) bridge instead of benzo (3-ring total) -
    # correctly out of this module's exactly-3-rings scope.
    mol = Chem.MolFromSmiles("C1=CC2=C(C=CC3=C2SC=C3)C4=C1C=CS4")
    assert not has_benzo_bis_heterocycle_fusion_name(mol)
