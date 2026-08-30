import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_sulfide():
    # 'sulfanyl' substituent-prefix naming mirrors 'oxy' (_ether.py, see
    # test_ether.py, already cross-checked against PubChem) with -O-
    # replaced by -S-, same P-14.3.4.2(b) omitted-locant convention.
    assert smiles_to_iupac("CSC") == "methylsulfanylmethane"


def test_methylsulfanylpropane():
    assert smiles_to_iupac("CSCCC") == "1-methylsulfanylpropane"


def test_diethyl_sulfide():
    assert smiles_to_iupac("CCSCC") == "ethylsulfanylethane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'sulfanyl' outside. Structure verified against PubChem: CID 522478
    # ("1-propan-2-ylsulfanylbutane" -- PubChem's own PIN-style name, this
    # project keeps its usual CAS-style substituent name instead).
    assert smiles_to_iupac("CCCCSC(C)C") == "1-(1-methylethyl)sulfanylbutane"


def test_both_sides_branched_and_tied_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)SC(C)C")
