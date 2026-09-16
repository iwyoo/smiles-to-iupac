import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_symmetric_alkyl_substituents():
    # PubChem PIN matches: dimethylphosphinic acid (CID 76777),
    # diethylphosphinic acid (CID 4186629), dipropylphosphinic acid
    # (CID 9124933).
    assert smiles_to_iupac("CP(C)(=O)O") == "dimethylphosphinic acid"
    assert smiles_to_iupac("CCP(CC)(=O)O") == "diethylphosphinic acid"
    assert smiles_to_iupac("CCCP(CCC)(=O)O") == "dipropylphosphinic acid"


def test_symmetric_benzene_rings():
    # PubChem PIN match: diphenylphosphinic acid (CID 15567).
    assert smiles_to_iupac("c1ccccc1P(c1ccccc1)(=O)O") == "diphenylphosphinic acid"


def test_asymmetric_substituents():
    # PubChem PIN matches: ethyl(methyl)phosphinic acid (CID 103893),
    # methyl(propan-2-yl)phosphinic acid (CID 14440203) --
    # `format_mononuclear_prefixes` orders and parenthesizes the two
    # different bare substituents alphabetically (P-16.5.1.3.1), no
    # P-locants used.
    assert smiles_to_iupac("CCP(C)(=O)O") == "ethyl(methyl)phosphinic acid"
    assert smiles_to_iupac("CC(C)P(C)(=O)O") == "methyl(propan-2-yl)phosphinic acid"


def test_halogenated_substituent():
    # Structure match: PubChem CID 56638447 (`CP(=O)(CCCl)O`). PubChem's
    # own auto-generated name omits the parentheses around the first,
    # compound substituent ('2-chloroethyl(methyl)phosphinic acid');
    # `format_mononuclear_prefixes` follows the literal P-16.5.1.3.1 text
    # instead (a compound first substituent is parenthesized too), the
    # same PubChem-vs-Blue-Book correction already documented for
    # `_phosphane.py`/`_phosphanone.py`.
    assert smiles_to_iupac("ClCCP(C)(=O)O") == "(2-chloroethyl)(methyl)phosphinic acid"


def test_rejects_second_phosphinic_acid_group():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OP(C)(=O)CCP(C)(=O)O")


def test_rejects_unrecognized_heteroatom():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCP(C)(=O)O")
