import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_symmetric_alkyl_substituents():
    # PubChem PIN matches: dimethylborinic acid (CID 5326208),
    # diethylborinic acid (CID 545119).
    assert smiles_to_iupac("CB(C)O") == "dimethylborinic acid"
    assert smiles_to_iupac("CCB(CC)O") == "diethylborinic acid"


def test_symmetric_benzene_rings():
    # PubChem PIN match: diphenylborinic acid (CID 17498).
    assert smiles_to_iupac("c1ccccc1B(c1ccccc1)O") == "diphenylborinic acid"


def test_asymmetric_substituents():
    # PubChem PIN matches: ethyl(methyl)borinic acid (CID 22096092),
    # methyl(propan-2-yl)borinic acid (CID 58801635) --
    # `format_mononuclear_prefixes` orders and parenthesizes the two
    # different bare substituents alphabetically (P-16.5.1.3.1), same
    # mechanism as `_phosphinic_acid.py`.
    assert smiles_to_iupac("CCB(C)O") == "ethyl(methyl)borinic acid"
    assert smiles_to_iupac("CC(C)B(C)O") == "methyl(propan-2-yl)borinic acid"


def test_halogenated_substituent():
    # Structure exists on PubChem with the same PubChem-vs-Blue-Book
    # parenthesization gap already documented in
    # `test_phosphinic_acid.py::test_halogenated_substituent`.
    assert smiles_to_iupac("ClCCB(C)O") == "(2-chloroethyl)(methyl)borinic acid"


def test_rejects_second_borinic_acid_group():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OB(C)CCB(C)O")


def test_rejects_unrecognized_heteroatom():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCB(C)O")
