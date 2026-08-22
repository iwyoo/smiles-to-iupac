import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._substituents import name_branch


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # isopropyl-like branch: root forks into two equal-length methyls, so
        # the longest chain from the fixed root (P-46) is 2 atoms (ethyl) with
        # the other methyl cited at locant 1.
        ("CCCC(C(C)C)CCC", "4-(1-methylethyl)heptane"),
        # tert-butyl-like branch: root forks into three methyls.
        ("CCCC(C(C)(C)C)CCC", "4-(1,1-dimethylethyl)heptane"),
        # sec-butyl-like branch: root forks into a methyl and an ethyl
        # continuation; the longer (ethyl) continuation wins the chain.
        ("CCCCC(C(C)CC)CCCCC", "5-(1-methylpropyl)decane"),
        # a compound substituent alongside a simple one on the same chain;
        # alphanumerical order (P-14.5.2) puts 'methyl' before the compound
        # substituent's own key ('methylethyl'), so it gets the lower locant.
        ("CCCC(C)C(C(C)C)CCC", "4-methyl-5-(1-methylethyl)octane"),
        # a compound substituent on a ring, alongside a simple one.
        ("CC1CCCCC1C(C)CC", "1-methyl-2-(1-methylpropyl)cyclohexane"),
        # two identical compound substituents: 'bis', not 'di' (P-14.2.2).
        (
            "CCCCC(C(C)CC)CCCCC(C(C)CC)CCCCC",
            "5,10-bis(1-methylpropyl)pentadecane",
        ),
    ],
)
def test_compound_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_substituent_raises():
    # Two separate cyclopropane rings joined by a chain: more than one ring
    # overall, so this is rejected before a cyclic substituent could even be
    # considered (see smiles_to_iupac.core's polycyclic check).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC1CCCCC1CC1")


def test_cyclic_substituent_raises():
    # Directly exercise the P-29.3.3 scope-out: a 3-membered ring (1-2-3-1)
    # hanging off atom 0 cannot be named as a compound substituent.
    graph = {0: [1], 1: [0, 2, 3], 2: [1, 3], 3: [1, 2]}
    with pytest.raises(UnsupportedStructure):
        name_branch(graph, 1, 0)
