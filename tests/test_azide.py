import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 79079: mononuclear parent (P-14.3.4.2(a), no locant).
        ("CN=[N+]=[N-]", "azidomethane"),
        # PubChem CID 79118: two-carbon chain, locant omitted
        # (P-14.3.4.2(b) -- unlike `_nitro.py`'s own two-carbon case, this
        # one matches PubChem's auto-generated name too, no discrepancy).
        ("CCN=[N+]=[N-]", "azidoethane"),
        # PubChem CID 12405073.
        ("CC(C)N=[N+]=[N-]", "2-azidopropane"),
        # PubChem CID 11469047: azide coexisting with a halogen substituent.
        ("ClCCN=[N+]=[N-]", "1-azido-2-chloroethane"),
        # PubChem CID 17884280: two azide groups (multiplying prefix).
        ("[N-]=[N+]=NCCCN=[N+]=[N-]", "1,3-diazidopropane"),
    ],
)
def test_azide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(N=[N+]=[N-])CC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCN=[N+]=[N-]")


def test_non_azide_triple_nitrogen_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNNN")
