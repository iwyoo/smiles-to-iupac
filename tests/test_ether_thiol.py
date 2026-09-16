import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID matches confirm the thiol locant is omitted at
        # chain length 2 regardless of the alkoxy substituent (unlike a
        # chain length 3+, e.g. plain `_thiol.py`'s 'propane-1-thiol').
        ("COCCS", "2-methoxyethanethiol"),
        ("CCOCCS", "2-ethoxyethanethiol"),
    ],
)
def test_smiles_to_iupac_ether_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCCS") == "2-(propan-2-yloxy)ethanethiol"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)CS") == "2-chloro-3-methoxypropane-1-thiol"


def test_two_thiols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCC(S)COC")


def test_two_ethers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COCC(OC)CS")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCCC1COC")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCC=CCOC")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S[C@@H](C)COC")


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_plain_thiol_still_works():
    assert smiles_to_iupac("CCCS") == "propane-1-thiol"
