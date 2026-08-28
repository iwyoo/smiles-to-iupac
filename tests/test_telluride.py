import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_telluride():
    # PubChem PUG REST CID 68977, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Te]C") == "methyltellanylmethane"


def test_diethyl_telluride():
    # PubChem PUG REST CID 69394, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Te]CC") == "ethyltellanylethane"


def test_methyl_ethyl_telluride():
    # The shorter (methyl) side becomes the substituent, the longer
    # (ethyl) side the parent. PubChem PUG REST CID 13981584,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("C[Te]CC") == "methyltellanylethane"


def test_methyl_propyl_telluride():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 15932889,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[Te]C") == "1-methyltellanylpropane"


def test_branched_tellanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te]C(C)C")
