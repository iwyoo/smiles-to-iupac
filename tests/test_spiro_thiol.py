import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 114818784, 168960680, 123897870
        # respectively). The parent stem's final 'e' is kept throughout
        # (see `test_von_baeyer_thiol.py`'s identical note).
        ("SC1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-thiol"),
        ("SC1CCC2(CC1)CCCCC2", "spiro[5.5]undecane-3-thiol"),
        ("SC1CCC2(CCC2)CC1", "spiro[3.5]nonane-7-thiol"),
    ],
)
def test_spiro_thiol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_spiro_thiols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC2(C1)CCCCC2S")


def test_spiro_thiol_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCC1CCCC12CCCCC2")


def test_unsaturated_spiro_thiol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC12C=CCCC2")
