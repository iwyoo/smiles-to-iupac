import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_plain_alkyl_chains():
    # PubChem PIN matches: methylphosphonic acid (CID 13818),
    # ethylphosphonic acid (CID 204482), propylphosphonic acid
    # (CID 204484), butylphosphonic acid (CID 76839), pentylphosphonic
    # acid (CID 302502), hexylphosphonic acid (CID 312552),
    # heptylphosphonic acid (CID 3870222), octylphosphonic acid
    # (CID 78452).
    assert smiles_to_iupac("CP(=O)(O)O") == "methylphosphonic acid"
    assert smiles_to_iupac("CCP(=O)(O)O") == "ethylphosphonic acid"
    assert smiles_to_iupac("CCCP(=O)(O)O") == "propylphosphonic acid"
    assert smiles_to_iupac("CCCCP(=O)(O)O") == "butylphosphonic acid"
    assert smiles_to_iupac("CCCCCP(=O)(O)O") == "pentylphosphonic acid"
    assert smiles_to_iupac("CCCCCCP(=O)(O)O") == "hexylphosphonic acid"
    assert smiles_to_iupac("CCCCCCCP(=O)(O)O") == "heptylphosphonic acid"
    assert smiles_to_iupac("CCCCCCCCP(=O)(O)O") == "octylphosphonic acid"


def test_branched_chains():
    # PubChem PIN matches: propan-2-ylphosphonic acid (CID 197187),
    # tert-butylphosphonic acid (CID 312546).
    assert smiles_to_iupac("CC(C)P(=O)(O)O") == "(propan-2-yl)phosphonic acid"
    assert smiles_to_iupac("CC(C)(C)P(=O)(O)O") == "tert-butylphosphonic acid"


def test_benzene_ring():
    # PubChem PIN match: phenylphosphonic acid (CID 15295) -- unlike the
    # sulfonic/selenonic/telluronic acid modules, a ring substituent uses
    # the exact same '<R-prefix>phosphonic acid' shape as an acyclic one
    # (P-67.1.1.2), no separate benzene-ring-retained-name path needed.
    assert smiles_to_iupac("c1ccccc1P(=O)(O)O") == "phenylphosphonic acid"


def test_saturated_ring():
    # PubChem PIN match: cyclohexylphosphonic acid (CID 70494).
    assert smiles_to_iupac("C1CCCCC1P(=O)(O)O") == "cyclohexylphosphonic acid"


def test_halogenated_chains():
    # PubChem PIN matches: chloromethylphosphonic acid (CID 75723),
    # 2-chloroethylphosphonic acid (CID 27982, "ethephon"),
    # 3-chloropropylphosphonic acid (CID 10920790), 2-bromoethylphosphonic
    # acid (CID 3613627), trichloromethylphosphonic acid (CID 80089).
    assert smiles_to_iupac("ClCP(=O)(O)O") == "(chloromethyl)phosphonic acid"
    assert smiles_to_iupac("ClCCP(=O)(O)O") == "(2-chloroethyl)phosphonic acid"
    assert smiles_to_iupac("ClCCCP(=O)(O)O") == "(3-chloropropyl)phosphonic acid"
    assert smiles_to_iupac("BrCCP(=O)(O)O") == "(2-bromoethyl)phosphonic acid"
    assert smiles_to_iupac("ClC(Cl)(Cl)P(=O)(O)O") == "(trichloromethyl)phosphonic acid"


def test_rejects_second_phosphonic_acid_group():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OP(=O)(O)CCP(=O)(O)O")


def test_rejects_unrecognized_heteroatom():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCP(=O)(O)O")
