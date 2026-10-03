import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem's own IUPACName is '2-methoxyacetaldehyde' (retained
        # 'acetaldehyde' stem); this project keeps its existing 'ethanal'
        # stem convention (see plain `_aldehyde.py`), so this differs from
        # PubChem only in that stem choice, not in substance.
        ("COCC=O", "2-methoxyethanal"),
        ("CCOCC=O", "2-ethoxyethanal"),
    ],
)
def test_smiles_to_iupac_ether_aldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCC=O") == "2-(propan-2-yloxy)ethanal"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C=O") == "2-chloro-3-methoxypropanal"


def test_ester_still_routes_correctly():
    assert smiles_to_iupac("COC(C)=O") == "methyl ethanoate"
    assert smiles_to_iupac("COC=O") == "methyl methanoate"


def test_two_aldehydes():
    assert smiles_to_iupac("O=CC(COC)C=O") == "2-(methoxymethyl)propanedial"


def test_two_ethers():
    assert smiles_to_iupac("COCC(OC)C=O") == "2,3-dimethoxypropanal"


def test_ring():
    assert smiles_to_iupac("O=CC1CCCCC1COC") == "2-(methoxymethyl)cyclohexane-1-carbaldehyde"


def test_specified_stereocenter():
    assert smiles_to_iupac("COC[C@@H](C)C=O") == "(2R)-3-methoxy-2-methylpropanal"


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_plain_aldehyde_still_works():
    assert smiles_to_iupac("CCC=O") == "propanal"
