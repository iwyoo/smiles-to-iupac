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
    assert smiles_to_iupac("CC(C)OCC(=O)OC") == "methyl 2-(propan-2-yl)oxyethanoate"


def test_halogen_on_acyl_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C(=O)OC") == "methyl 2-chloro-3-methoxypropanoate"


def test_plain_ester_still_routes_correctly():
    assert smiles_to_iupac("COC(=O)C") == "methyl ethanoate"
    assert smiles_to_iupac("CC(=O)OC") == "methyl ethanoate"


def test_ether_on_alcohol_part_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OCCOC")


def test_two_esters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COCC(=O)OCOC(=O)C")


def test_two_ethers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COCC(OC)C(=O)OC")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(OC)C1CCCCC1COC")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC[C@@H](C)C(=O)OC")


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"
