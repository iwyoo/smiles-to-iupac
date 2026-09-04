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


def test_phenyl_nitro_direct_bond():
    # P-44.1.2.2 rule (1): 'nitro' has no suffix form, so the ring is
    # always senior to a chain of the same class -- confirmed by PubChem
    # CID 7416.
    assert smiles_to_iupac("c1ccccc1[N+](=O)[O-]") == "nitrobenzene"


def test_phenyl_nitro_chain():
    # PubChem CID 80208 gives "2-nitroethylbenzene" (no parentheses), but
    # this codebase follows `_azide.py`'s identical, Blue-Book-verified
    # rule instead (the Blue Book's own worked example is
    # '(2-azidoethyl)benzene (PIN)', parenthesized) -- a locant-bearing
    # compound substituent prefix is enclosed regardless of which simple
    # prefix (azido/nitro/...) it carries, so 'nitro' gets the same
    # treatment PubChem's own algorithmic name doesn't apply consistently
    # here.
    assert smiles_to_iupac("c1ccccc1CC[N+](=O)[O-]") == "(2-nitroethyl)benzene"
    # PubChem CID 54143589 gives "7-nitroheptylbenzene" (same
    # no-parentheses inconsistency); ring wins regardless of chain
    # length, same as 'heptylbenzene (PIN)'.
    assert smiles_to_iupac("c1ccccc1CCCCCCC[N+](=O)[O-]") == "(7-nitroheptyl)benzene"


def test_phenyl_nitro_halogen_coexistence():
    assert smiles_to_iupac("c1ccccc1CC(Cl)[N+](=O)[O-]") == "(2-chloro-2-nitroethyl)benzene"


def test_phenyl_nitro_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1C[N+](=O)[O-]")


def test_phenyl_nitro_unsaturation_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1C[N+](=O)[O-]")
