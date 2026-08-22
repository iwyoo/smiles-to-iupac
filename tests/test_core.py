import pytest

from chemonym import smiles_to_iupac


def test_smiles_to_iupac_not_implemented():
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("CCO")
