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


def test_dithieno_pyridine_second_real_isomer():
    # PubChem CID 10932260 - found during #607's broader PubChem
    # validation pass; PubChem's own name for it lacks primes
    # ("dithieno[3,2-b:2,3-e]pyridine"), but the primed form is the
    # correct PIN per P-25.3.4.1.1 (this project matches other cases
    # where PubChem's synonym is a looser, non-PIN form).
    mol = Chem.MolFromSmiles("C1=CSC2=C1N=C3C=CSC3=C2")
    assert name_pyridine_bis_heterocycle_fusion(mol) == "dithieno[3,2-b:2',3'-e]pyridine"


def test_dithieno_pyridine_third_real_isomer():
    # PubChem CID 901680 - letter 'd' instead of 'e', found in the same
    # validation pass.
    mol = Chem.MolFromSmiles("C1=CSC2=CN=C3C(=C21)C=CS3")
    assert name_pyridine_bis_heterocycle_fusion(mol) == "dithieno[2,3-b:3',2'-d]pyridine"


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
