import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_guanidine():
    # PubChem structure match: "guanidine".
    assert smiles_to_iupac("NC(=N)N") == "guanidine"


def test_n_methylguanidine():
    # PubChem structure match ("2-methylguanidine", deprecated numeral
    # locant); this project uses the Blue Book's own letter-locant style
    # established by `_urea.py`/`_thiourea.py`.
    assert smiles_to_iupac("CNC(=N)N") == "N-methylguanidine"


def test_n_ethylguanidine():
    assert smiles_to_iupac("CCNC(=N)N") == "N-ethylguanidine"


def test_n_n_dimethylguanidine_same_nitrogen():
    # PubChem structure match: "1,1-dimethylguanidine".
    assert smiles_to_iupac("CN(C)C(=N)N") == "N,N-dimethylguanidine"


def test_n_ethyl_n_methylguanidine_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=N)N") == "N-ethyl-N-methylguanidine"


def test_n_n_prime_dimethylguanidine_different_nitrogens():
    # Blue Book P-66.4.1.2.1.2 worked example: "N,N′-dimethylguanidine
    # (PIN) (not N,N′′-dimethylguanidine)". PubChem structure match:
    # "1,2-dimethylguanidine" (deprecated numeral locants).
    assert smiles_to_iupac("CNC(=N)NC") == "N,N'-dimethylguanidine"


def test_different_substituents_on_different_amino_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNC(=N)NC")


def test_different_substituent_counts_on_amino_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(C)C(=N)NC")


def test_imino_nitrogen_substitution_not_supported():
    # N''-substitution (on the double-bonded nitrogen) is a confirmed Blue
    # Book shape ("N,N,N′,N′-tetramethyl-N′′-phenylguanidine (PIN)") but
    # deferred as its own follow-up.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=NC)N")


def test_branched_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)NC(=N)N")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=N)N")


def test_amidine_not_confused_with_guanidine():
    assert smiles_to_iupac("CC(=N)N") == "ethanimidamide"


def test_urea_not_confused_with_guanidine():
    assert smiles_to_iupac("NC(=O)N") == "urea"
