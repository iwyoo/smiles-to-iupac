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


def test_disulfide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CSSC")


def test_branched_sulfanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)SC(C)C")
