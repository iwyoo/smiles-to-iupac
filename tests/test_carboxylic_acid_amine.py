import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # GABA (gamma-aminobutyric acid): PubChem-verified exactly
        # ('4-aminobutanoic acid', CID 119).
        ("NCCCC(=O)O", "4-aminobutanoic acid"),
        # beta-alanine: PubChem-verified exactly ('3-aminopropanoic acid',
        # CID 239).
        ("C(CN)C(=O)O", "3-aminopropanoic acid"),
        # 6-aminohexanoic acid: PubChem-verified exactly (CID 564).
        ("C(CCC(=O)O)CCN", "6-aminohexanoic acid"),
        # 2-aminoisobutyric acid (AIB): both the amine and a methyl branch
        # sit at C2 -- PubChem-verified exactly
        # ('2-amino-2-methylpropanoic acid', CID 6119), proving the amine
        # and an ordinary substituent group and alphabetize/cite together
        # correctly.
        ("CC(C)(C(=O)O)N", "2-amino-2-methylpropanoic acid"),
        # glycine: the amine and acid on the same (only) chain carbon.
        # PubChem's own PIN ('2-aminoacetic acid', CID 750) uses the
        # retained 'acetic acid' stem; this project's own
        # `_carboxylic_acid.py` always uses the systematic 'ethanoic acid'
        # stem instead (see that module's tests), so this module follows
        # that same pre-existing convention rather than introducing a new
        # divergence (see module docstring).
        ("NCC(=O)O", "2-aminoethanoic acid"),
        # a halogen substituent coexists with both the acid and the amine --
        # same reasoning as glycine above for the acid-stem divergence
        # (PubChem CID 10308266 gives '2-amino-2-chloroacetic acid').
        ("NC(Cl)C(=O)O", "2-amino-2-chloroethanoic acid"),
    ],
)
def test_smiles_to_iupac_carboxylic_acid_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_carbamic_acid_not_this_module():
    # H2N-COOH: the amine nitrogen is bonded directly to the acid carbon
    # itself, not to some other chain carbon -- carbamic acid, a distinct
    # retained functional class (_carbamate.py's territory), not this
    # module's "amine elsewhere on the chain" shape. Pre-existing test
    # (tests/test_carbamate.py::test_free_carbamic_acid_not_supported)
    # confirms it's still rejected overall.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(N)=O")


def test_secondary_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNCC(=O)O")


def test_two_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(N)CC(=O)O")


def test_two_carboxylic_acids_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(C(=O)O)C(=O)O")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCC(C(=O)O)CC1")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC=CC(=O)O")


def test_plain_carboxylic_acid_still_works():
    assert smiles_to_iupac("CCC(=O)O") == "propanoic acid"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"
