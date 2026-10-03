import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-28.3), N=3 spanning all three
        # substitution patterns.
        ("C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=CC=C3", "11,21:24,31-terphenyl"),  # CID 7115, p-terphenyl (PIN worked example)
        ("C1=CC=C(C=C1)C2=CC(=CC=C2)C3=CC=CC=C3", "11,21:23,31-terphenyl"),  # CID 7076, m-terphenyl
        ("C1=CC=C(C=C1)C2=CC=CC=C2C3=CC=CC=C3", "11,21:22,31-terphenyl"),  # CID 6766, o-terphenyl
        # Real PubChem structures, N=4-6, all-para chain.
        ("C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=C(C=C3)C4=CC=CC=C4", "11,21:24,31:34,41-quaterphenyl"),  # CID 8677
        (
            "C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=C(C=C3)C4=CC=C(C=C4)C5=CC=CC=C5",
            "11,21:24,31:34,41:44,51-quinquephenyl",
        ),  # CID 137813
        (
            "C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=C(C=C3)C4=CC=C(C=C4)C5=CC=C(C=C5)C6=CC=CC=C6",
            "11,21:24,31:34,41:44,51:54,61-sexiphenyl",
        ),  # CID 78254
        # Real PubChem structure, N=4, both middle rings meta-linked -- the
        # primary source's own tie-break worked example ("not
        # 11,21:23,31:32,41-quaterphenyl").
        (
            "C1=CC=C(C=C1)C2=CC(=CC=C2)C3=CC=CC(=C3)C4=CC=CC=C4",
            "11,21:23,31:33,41-quaterphenyl",
        ),  # CID 14422
    ],
)
def test_ring_assembly_chain_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituent_lowest_locant():
    # Chlorine on a terminal ring, para to the attachment point.
    assert smiles_to_iupac("c1cc(Cl)ccc1-c1ccc(-c2ccccc2)cc1") == "14-chloro-11,21:24,31-terphenyl"


def test_biphenyl_still_routes_to_two_ring_module():
    assert smiles_to_iupac("c1ccccc1-c1ccccc1") == "1,1'-biphenyl"


def test_two_rings_plus_extra_ring_still_out_of_scope():
    # Fused two-ring aromatic (naphthalene) is a different module's shape
    # entirely, unaffected by this one.
    assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


def test_branched_three_ring_assembly_raises():
    # 1,3,5-triphenylbenzene -- three rings all attached to a central ring,
    # not an unbranched chain (P-28.6, out of scope).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc(-c2cc(-c3ccccc3)cc(-c4ccccc4)c2)cc1")


def test_seven_ring_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(
            "c1ccc(-c2ccc(-c3ccc(-c4ccc(-c5ccc(-c6ccc(-c7ccccc7)cc6)cc5)cc4)cc3)cc2)cc1"
        )


def test_non_benzene_ring_chain_routes_to_saturated_case():
    # Cycloalkane chains are `_ring_assembly_chain.py`'s own job too now
    # (see tests/test_ring_assembly_cycloalkane_chain.py) -- this just
    # confirms the aromatic-only branch no longer mis-fails on it.
    assert smiles_to_iupac("C1CC1C1CC1C1CC1") == "11,21:22,31-tercyclopropane"


def test_mixed_ring_kinds_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC1c1ccc(-c2ccccc2)cc1")


def test_mixed_ring_sizes():
    assert smiles_to_iupac("C1CC1C1CCCC1C1CC1") == "1,2-dicyclopropylcyclopentane"
