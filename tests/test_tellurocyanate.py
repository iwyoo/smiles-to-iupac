import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_tellurocyanate():
    # PubChem structure match: "methyl tellurocyanate" -- the only
    # tellurocyanate PubChem has registered at all (ethyl/propyl
    # candidates both come back as CID 0).
    assert smiles_to_iupac("C[Te]C#N") == "methyl tellurocyanate"


def test_ethyl_tellurocyanate():
    assert smiles_to_iupac("CC[Te]C#N") == "ethyl tellurocyanate"


def test_propyl_tellurocyanate():
    assert smiles_to_iupac("CCC[Te]C#N") == "propyl tellurocyanate"


def test_branched_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te]C#N")


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Te]C#N")


def test_cyclic_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1[Te]C#N")


def test_two_tellurocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#C[Te]CC[Te]C#N")


def test_selenocyanate_not_confused_with_tellurocyanate():
    assert smiles_to_iupac("C[Se]C#N") == "methyl selenocyanate"
