import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: 'COCC(C)=O' -> '1-methoxypropan-2-one'.
        ("COCC(C)=O", "1-methoxypropan-2-one"),
        ("COCC(CC)=O", "1-methoxybutan-2-one"),
    ],
)
def test_smiles_to_iupac_ether_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCC(C)=O") == "1-(propan-2-yloxy)propan-2-one"


def test_branched_alkoxy_r_prime_on_main_chain():
    # Found via real-data testing (smiles-to-iupac-realdata-test's pubchem
    # diff): this used to come out as
    # '5-(propan-2-yl)oxy-4,4-dimethylpentan-2-one', 'oxy' wrongly sitting
    # outside the enclosing parenthesis instead of fused inside it -- see
    # `test_ether.py`'s own note on the same bug across all 8 ether
    # modules. PubChem-verified reference:
    # '4,4-dimethyl-5-propan-2-yloxypentan-2-one'.
    assert smiles_to_iupac("CC(=O)CC(C)(C)COC(C)C") == "4,4-dimethyl-5-(propan-2-yloxy)pentan-2-one"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C(C)=O") == "3-chloro-4-methoxybutan-2-one"


def test_ester_still_routes_correctly():
    assert smiles_to_iupac("COC(C)=O") == "methyl ethanoate"


def test_two_ketones():
    assert smiles_to_iupac("O=CC(=O)COC") == "3-methoxy-2-oxopropanal"


def test_two_ethers():
    assert smiles_to_iupac("COCC(OC)C(C)=O") == "3,4-dimethoxybutan-2-one"


def test_ring():
    assert smiles_to_iupac("O=C1CCCCC1COC") == "2-(methoxymethyl)cyclohexan-1-one"


def test_specified_stereocenter():
    assert smiles_to_iupac("COC[C@@H](C)C(C)=O") == "(3R)-4-methoxy-3-methylbutan-2-one"


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_plain_ketone_still_works():
    assert smiles_to_iupac("CCC(C)=O") == "butan-2-one"
