import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 70075: mononuclear parent (P-14.3.4.2(a), no locant).
        ("CN=O", "nitrosomethane"),
        # PubChem CID 79124: two-carbon chain, locant omitted (same
        # P-14.3.4.2(b) rule `_nitro.py`'s own docstring confirms for
        # 'nitroethane', without the PubChem-vs-PIN mismatch nitro has).
        ("CCN=O", "nitrosoethane"),
        # PubChem CID 21528903: nitroso coexisting with a halogen
        # substituent.
        ("ClCCN=O", "1-chloro-2-nitrosoethane"),
        # PubChem CID 13116308: two nitroso groups (multiplying prefix).
        ("O=NCCN=O", "1,2-dinitrosoethane"),
    ],
)
def test_nitroso(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nitro_not_misnamed_as_nitroso():
    assert smiles_to_iupac("C[N+](=O)[O-]") == "nitromethane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1N=O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CN=O")
