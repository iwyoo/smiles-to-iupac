import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._furan_pyran_fusion import has_furan_pyran_fusion_name


def test_detects_furan_pyran_shape():
    mol = Chem.MolFromSmiles("C1=Cc2occc2OC1")
    assert has_furan_pyran_fusion_name(mol)


def test_refuses_to_guess_indicated_hydrogen():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=Cc2occc2OC1")
