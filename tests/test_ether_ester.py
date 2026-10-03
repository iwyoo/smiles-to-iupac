import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem's own IUPACName is 'methyl 2-methoxyacetate' (retained
        # 'acetate' stem); this project keeps its existing 'ethanoate'
        # stem convention (see plain `_ester.py`), so this differs from
        # PubChem only in that stem choice, not in substance.
        ("COCC(=O)OC", "methyl 2-methoxyethanoate"),
        ("CCOCC(=O)OC", "methyl 2-ethoxyethanoate"),
    ],
)
def test_smiles_to_iupac_ether_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCC(=O)OC") == "methyl 2-(propan-2-yloxy)ethanoate"


def test_halogen_on_acyl_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C(=O)OC") == "methyl 2-chloro-3-methoxypropanoate"


def test_plain_ester_still_routes_correctly():
    assert smiles_to_iupac("COC(=O)C") == "methyl ethanoate"
    assert smiles_to_iupac("CC(=O)OC") == "methyl ethanoate"


def test_ether_on_alcohol_part():
    assert smiles_to_iupac("CC(=O)OCCOC") == "2-methoxyethyl ethanoate"


def test_two_esters_names_polyester():
    assert smiles_to_iupac("COCC(=O)OCOC(=O)C") == "methylene ethanoate 2-methoxyethanoate"


def test_two_ethers():
    assert smiles_to_iupac("COCC(OC)C(=O)OC") == "methyl 2,3-dimethoxypropanoate"


def test_ring():
    assert smiles_to_iupac("O=C(OC)C1CCCCC1COC") == "methyl 2-(methoxymethyl)cyclohexane-1-carboxylate"


def test_specified_stereocenter():
    assert smiles_to_iupac("COC[C@@H](C)C(=O)OC") == "methyl (2R)-3-methoxy-2-methylpropanoate"


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"
