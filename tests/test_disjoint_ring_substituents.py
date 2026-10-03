import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCCC1CCC1CCCCC1", "1,1'-(ethane-1,2-diyl)dicyclohexane"),  # multiplicative PIN, P-15.3
        ("C1CCCCC1CCCc1ccccc1", "(3-cyclohexylpropyl)benzene"),  # CID 561990
        ("c1ccccc1CCCc1ccccc1", "1,1'-(propane-1,3-diyl)dibenzene"),  # multiplicative PIN, P-15.3
        ("C1CCCCC1CC1CCCCC1", "1,1'-methylenedicyclohexane"),  # multiplicative PIN, P-15.3
        ("c1ccsc1Cc1ccsc1", "2-[(thiophen-3-yl)methyl]thiophene"),
        ("c1ccsc1CCC1CCCCC1", "2-(2-cyclohexylethyl)thiophene"),
    ],
)
def test_disjoint_ring_pair_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_rings_joined_by_double_bond_use_the_larger_ring_as_parent():
    assert smiles_to_iupac("C1CCC1=C1CCCC1") == "cyclobutylidenecyclopentane"


def test_three_disjoint_rings_multiply_the_two_terminal_rings():
    assert (
        smiles_to_iupac("C1CCCCC1CCC1CCCCC1CCC1CCCCC1")
        == "1,1'-[cyclohexane-1,2-diyldi(ethane-2,1-diyl)]dicyclohexane"
    )


def test_ring_with_extra_substituent():
    assert smiles_to_iupac("CC1CCCCC1CC1CCCCC1") == "1-(cyclohexylmethyl)-2-methylcyclohexane"


def test_unsaturated_ring_in_a_pair_is_not_misnamed_as_saturated_is_named():
    assert smiles_to_iupac("C1=CCCCC1CC1=CCCCC1") == "1-[(cyclohex-2-en-1-yl)methyl]cyclohex-1-ene"
