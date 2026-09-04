import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_thiourea():
    # PubChem structure match: "thiourea".
    assert smiles_to_iupac("NC(=S)N") == "thiourea"


def test_n_methylthiourea():
    # PubChem structure match ("methylthiourea"); this project uses the
    # Blue Book's own letter-locant style established by `_urea.py`, not
    # PubChem's own convention.
    assert smiles_to_iupac("CNC(=S)N") == "N-methylthiourea"


def test_n_n_dimethylthiourea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=S)N") == "N,N-dimethylthiourea"


def test_n_ethyl_n_methylthiourea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=S)N") == "N-ethyl-N-methylthiourea"


def test_n_n_prime_dimethylthiourea_different_nitrogens():
    assert smiles_to_iupac("CNC(=S)NC") == "N,N'-dimethylthiourea"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # One different substituent on each nitrogen -- the alphabetically
        # first substituent name becomes 'N-', the other 'N''-', same rule
        # as `_urea.py` (PubChem structure match, its own numeric-locant
        # style): `CCNC(=S)NC` -> '1-ethyl-3-methylthiourea' (CID
        # 15568242), `CCCNC(=S)NC` -> '1-methyl-3-propylthiourea' (CID
        # 3690223).
        ("CCNC(=S)NC", "N-ethyl-N'-methylthiourea"),
        ("CCCNC(=S)NC", "N-methyl-N'-propylthiourea"),
    ],
)
def test_different_substituents_on_different_nitrogens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_different_substituent_counts_on_different_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCN(C)C(=S)NC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 1711921, 737374.
        ("CC(C)NC(=S)N", "N-(propan-2-yl)thiourea"),
        ("CC(C)(C)NC(=S)N", "N-tert-butylthiourea"),
    ],
)
def test_branched_n_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_identical_branched_substituents_different_nitrogens_parenthesized_when_compound():
    # PubChem CID 2758386.
    assert smiles_to_iupac("CC(C)NC(=S)NC(C)C") == "N,N'-di(propan-2-yl)thiourea"


def test_two_identical_branched_substituents_different_nitrogens_not_parenthesized_when_retained():
    # PubChem CID 2801221.
    assert smiles_to_iupac("CC(C)(C)NC(=S)NC(C)(C)C") == "N,N'-ditert-butylthiourea"


def test_two_different_substituents_alphabetized_ignoring_italic_prefix():
    # PubChem CID 4611181.
    assert smiles_to_iupac("CC(C)(C)NC(=S)NCC") == "N-tert-butyl-N'-ethylthiourea"


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=S)N")


def test_urea_not_confused_with_thiourea():
    assert smiles_to_iupac("NC(=O)N") == "urea"


def test_n_phenylthiourea():
    # PubChem structure match: "phenylthiourea" (CID 676454).
    assert smiles_to_iupac("NC(=S)Nc1ccccc1") == "N-phenylthiourea"


def test_n_methyl_n_prime_phenylthiourea_different_nitrogens():
    # PubChem structure match: "1-methyl-3-phenylthiourea" (CID 698294).
    assert smiles_to_iupac("CNC(=S)Nc1ccccc1") == "N-methyl-N'-phenylthiourea"


def test_n_n_prime_diphenylthiourea():
    # PubChem structure match: "1,3-diphenylthiourea" (CID 700999).
    assert smiles_to_iupac("c1ccc(NC(=S)Nc2ccccc2)cc1") == "N,N'-diphenylthiourea"


def test_substituted_phenyl_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=S)Nc1ccc(C)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=S)N")
