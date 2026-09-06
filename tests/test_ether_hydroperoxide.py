import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem's own IUPACName is "1-hydroperoxy-2-methoxyethane" --
        # naming the molecule as a substituted ethane with both groups as
        # prefixes, not using the 'peroxol' suffix at all. This
        # contradicts P-41/P-43's own text (a hydroperoxide is always
        # cited via the 'peroxol' suffix when it's the senior/only
        # suffix-eligible group present); this project follows the Blue
        # Book primary source directly, the same policy
        # `_hydroperoxide_amine.py` established for this class (PR #426).
        ("COCCOO", "2-methoxyethane-1-peroxol"),
        ("CCOCCOO", "2-ethoxyethane-1-peroxol"),
    ],
)
def test_smiles_to_iupac_ether_hydroperoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCCOO") == "2-(propan-2-yloxy)ethane-1-peroxol"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)COO") == "2-chloro-3-methoxypropane-1-peroxol"


def test_two_hydroperoxides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOCC(OO)COC")


def test_two_ethers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COCC(OC)COO")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOC1CCCCC1COC")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOCC=CCOC")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OO[C@@H](C)COC")


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_plain_hydroperoxide_still_works():
    assert smiles_to_iupac("CCOO") == "ethaneperoxol"
