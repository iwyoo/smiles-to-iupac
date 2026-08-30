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


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'tellanyl' outside. Structure verified against PubChem: CID 13975014
    # ("1-propan-2-yltellanylbutane" -- PubChem's own PIN-style name, this
    # project keeps its usual CAS-style substituent name instead).
    assert smiles_to_iupac("CCCC[Te]C(C)C") == "1-(1-methylethyl)tellanylbutane"


def test_both_sides_branched_and_tied_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te]C(C)C")
