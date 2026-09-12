from rdkit import Chem

from smiles_to_iupac._pyridine_bis_heterocycle_fusion import (
    has_pyridine_bis_heterocycle_fusion_name,
    name_pyridine_bis_heterocycle_fusion,
)


def test_dithieno_pyridine_matches_pubchem():
    # PubChem CID 129867025.
    mol = Chem.MolFromSmiles("C1=CSC2=NC3=C(C=CS3)C=C21")
    assert has_pyridine_bis_heterocycle_fusion_name(mol)
    assert name_pyridine_bis_heterocycle_fusion(mol) == "dithieno[2,3-b:3',2'-e]pyridine"


def test_difuro_pyridine_analog():
    # Constructed by swapping S for O in the same skeleton as the real
    # dithieno compound above - not independently registered on PubChem,
    # but the same mechanism (see module docstring's honest breadth
    # caveat).
    mol = Chem.MolFromSmiles("C1=COC2=NC3=C(C=CO3)C=C21")
    assert name_pyridine_bis_heterocycle_fusion(mol) == "difuro[2,3-b:3',2'-e]pyridine"


def test_mixed_furan_thiophene_out_of_scope():
    mol = Chem.MolFromSmiles("C1=COC2=NC3=C(C=CS3)C=C21")
    assert not has_pyridine_bis_heterocycle_fusion_name(mol)


def test_attached_rings_fused_to_each_other_out_of_scope():
    # dithieno[3,2-b:2',3'-d]thiophene (PubChem CID 137985, "DTT"): no
    # pyridine ring at all, and the two thiophenes are fused to each
    # other via a middle thiophene, not to a shared pyridine parent.
    mol = Chem.MolFromSmiles("C1=CSC2=C1SC3=C2SC=C3")
    assert not has_pyridine_bis_heterocycle_fusion_name(mol)


def test_wrong_ring_count_out_of_scope():
    mol = Chem.MolFromSmiles("C1=CC2=C(C=CS2)N=C1")  # thieno[3,2-b]pyridine, 2 rings
    assert not has_pyridine_bis_heterocycle_fusion_name(mol)


def test_substituent_out_of_scope():
    mol = Chem.MolFromSmiles("Cc1csc2nc3ccsc3cc12")
    assert not has_pyridine_bis_heterocycle_fusion_name(mol)
