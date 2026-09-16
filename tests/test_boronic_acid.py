import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_plain_alkyl_chains():
    # PubChem PIN matches: methylboronic acid (CID 139377), ethylboronic
    # acid (CID 521157), propylboronic acid (CID 351065), butylboronic
    # acid (CID 20479), pentylboronic acid (CID 352103), hexylboronic
    # acid (CID 351064).
    assert smiles_to_iupac("CB(O)O") == "methylboronic acid"
    assert smiles_to_iupac("CCB(O)O") == "ethylboronic acid"
    assert smiles_to_iupac("CCCB(O)O") == "propylboronic acid"
    assert smiles_to_iupac("CCCCB(O)O") == "butylboronic acid"
    assert smiles_to_iupac("CCCCCB(O)O") == "pentylboronic acid"
    assert smiles_to_iupac("CCCCCCB(O)O") == "hexylboronic acid"


def test_branched_chains():
    # PubChem PIN matches: propan-2-ylboronic acid (CID 2734750),
    # tert-butylboronic acid (CID 5156670).
    assert smiles_to_iupac("CC(C)B(O)O") == "(propan-2-yl)boronic acid"
    assert smiles_to_iupac("CC(C)(C)B(O)O") == "tert-butylboronic acid"


def test_benzene_ring():
    # PubChem PIN match: phenylboronic acid (CID 66827) -- like
    # `_phosphonic_acid.py`, a ring substituent uses the exact same
    # '<R-prefix>boronic acid' shape as an acyclic one (P-68.1.4.1).
    assert smiles_to_iupac("c1ccccc1B(O)O") == "phenylboronic acid"


def test_saturated_ring():
    # PubChem PIN match: cyclohexylboronic acid (CID 199578).
    assert smiles_to_iupac("C1CCCCC1B(O)O") == "cyclohexylboronic acid"


def test_halogenated_chains():
    # PubChem PIN matches: chloromethylboronic acid (CID 15828379),
    # 2-chloroethylboronic acid (CID 86052048), 2-bromoethylboronic acid
    # (CID 20826997).
    assert smiles_to_iupac("ClCB(O)O") == "(chloromethyl)boronic acid"
    assert smiles_to_iupac("ClCCB(O)O") == "(2-chloroethyl)boronic acid"
    assert smiles_to_iupac("BrCCB(O)O") == "(2-bromoethyl)boronic acid"


def test_rejects_second_boronic_acid_group():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OB(O)CCB(O)O")


def test_rejects_unrecognized_heteroatom():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCB(O)O")
