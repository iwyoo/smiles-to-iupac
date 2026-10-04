import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._pyrimidinedione import has_pyrimidinedione_shape


def test_uracil_name():
    assert smiles_to_iupac("O=c1cc[nH]c(=O)[nH]1") == "pyrimidine-2,4(1H,3H)-dione"


def test_thymine_name():
    # Thymine (5-methyluracil), PubChem CID 1135 isomeric SMILES.
    assert smiles_to_iupac("CC1=CNC(=O)NC1=O") == "5-methylpyrimidine-2,4(1H,3H)-dione"


def test_6_methyl_regiochemistry_gets_locant_6_not_5():
    # A methyl on the ring carbon adjacent to the *other* nitrogen (not
    # adjacent to any oxo carbon) is a distinct regiochemistry from
    # thymine's 5-methyl and must not be conflated with it.
    assert smiles_to_iupac("O=c1[nH]c(=O)cc(C)[nH]1") == "6-methylpyrimidine-2,4(1H,3H)-dione"


def test_pyrimidinone_single_oxo_unaffected():
    assert smiles_to_iupac("O=c1nccc[nH]1") == "pyrimidin-2(1H)-one"


def test_barbituric_acid_names_as_saturated_trione():
    assert smiles_to_iupac("O=C1CC(=O)NC(=O)N1") == "1,3-diazinane-2,4,6-trione"


def test_non_1_3_diazine_dione_not_matched():
    mol = Chem.MolFromSmiles("O=C1C=CC(=O)NN1")
    assert not has_pyrimidinedione_shape(mol)
