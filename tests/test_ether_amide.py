import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem's own IUPACName is '2-methoxyacetamide' (retained
        # 'acetamide' stem); this project keeps its existing 'ethanamide'
        # stem convention (see plain `_amide.py`), so this differs from
        # PubChem only in that stem choice, not in substance.
        ("COCC(N)=O", "2-methoxyethanamide"),
        ("CCOCC(N)=O", "2-ethoxyethanamide"),
    ],
)
def test_smiles_to_iupac_ether_amide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCC(N)=O") == "2-(propan-2-yloxy)ethanamide"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C(N)=O") == "2-chloro-3-methoxypropanamide"


def test_carbamate_still_routes_correctly():
    assert smiles_to_iupac("COC(N)=O") == "methyl carbamate"
    assert smiles_to_iupac("CCOC(N)=O") == "ethyl carbamate"


def test_n_alkyl_amide():
    assert smiles_to_iupac("CNC(=O)COC") == "2-methoxy-N-methylethanamide"


def test_two_amides():
    assert smiles_to_iupac("NC(=O)C(COC)C(N)=O") == "2-(methoxymethyl)propanediamide"


def test_two_ethers():
    assert smiles_to_iupac("COCC(OC)C(N)=O") == "2,3-dimethoxypropanamide"


def test_ring():
    assert smiles_to_iupac("NC(=O)C1CCCCC1COC") == "2-(methoxymethyl)cyclohexane-1-carboxamide"


def test_specified_stereocenter():
    assert smiles_to_iupac("COC[C@@H](C)C(N)=O") == "(2R)-3-methoxy-2-methylpropanamide"


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_plain_amide_still_works():
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"
