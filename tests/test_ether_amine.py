import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID matches confirm the amine locant is omitted at
        # chain length 2 regardless of the alkoxy substituent (unlike a
        # chain length 3+, e.g. plain `_amine.py`'s 'propan-1-amine').
        ("COCCN", "2-methoxyethanamine"),
        ("CCOCCN", "2-ethoxyethanamine"),
        ("CCCOCCN", "2-propoxyethanamine"),
    ],
)
def test_smiles_to_iupac_ether_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_alkoxy_r_prime():
    # Mirrors `_ether.py`'s own enclosure convention: the whole 'yloxy'
    # group (including 'oxy') is parenthesized when it carries its own
    # internal locant, e.g. '(butan-2-yloxy)' -- P-16.3.3.
    assert smiles_to_iupac("CC(C)OCCN") == "2-(propan-2-yloxy)ethanamine"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)CN") == "2-chloro-3-methoxypropan-1-amine"


def test_two_amines():
    assert smiles_to_iupac("NCC(N)COC") == "3-methoxypropane-1,2-diamine"


def test_two_ethers():
    assert smiles_to_iupac("COCC(OC)CN") == "2,3-dimethoxypropan-1-amine"


def test_ring():
    assert smiles_to_iupac("NC1CCCCC1COC") == "2-(methoxymethyl)cyclohexan-1-amine"


def test_unsaturated_chain():
    assert smiles_to_iupac("NCC=CCOC") == "4-methoxybut-2-en-1-amine"


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](C)COC")


def test_hydroxyl_coexisting():
    assert smiles_to_iupac("NCC(O)COC") == "1-amino-3-methoxypropan-2-ol"


def test_plain_ether_still_works():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"
