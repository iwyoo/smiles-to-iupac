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


def test_branched_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)NC(=S)N")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=S)N")


def test_urea_not_confused_with_thiourea():
    assert smiles_to_iupac("NC(=O)N") == "urea"
