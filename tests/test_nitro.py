import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 6375: mononuclear parent (P-14.3.4.2(a), no locant).
        ("C[N+](=O)[O-]", "nitromethane"),
        # PubChem CID 7903: three-carbon chain, locant required (only one
        # substituent position is lowest, unlike the two-carbon case below).
        ("CCC[N+](=O)[O-]", "1-nitropropane"),
        # PubChem CID 398.
        ("CC(C)[N+](=O)[O-]", "2-nitropropane"),
        # PubChem CID 136442: nitro coexisting with a halogen substituent.
        ("ClCC[N+](=O)[O-]", "1-chloro-2-nitroethane"),
        # PubChem CID 138670: two nitro groups (multiplying prefix).
        ("O=[N+]([O-])CCC[N+](=O)[O-]", "1,3-dinitropropane"),
    ],
)
def test_nitro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_carbon_chain_omits_locant():
    # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly one
    # substituent omits the locant, the same rule already verified for
    # halogens (test_halogens.py's 'chloroethane' for CCCl) -- PubChem's
    # own auto-generated name for this exact structure (CID 6587) is
    # '1-nitroethane', which contradicts this general, unconditional rule
    # (P1.pdf) and every other single-substituent two-carbon-chain case
    # this project has verified; treated here as another instance of the
    # known PubChem-autoname-vs-PIN mismatch documented throughout this
    # project (see e.g. `_hydroxylamine.py`, `_sulfoxide.py`).
    assert smiles_to_iupac("CC[N+](=O)[O-]") == "nitroethane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N+](=O)[O-])CC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N+](=O)[O-]")
