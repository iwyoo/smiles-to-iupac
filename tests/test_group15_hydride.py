import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book worked examples (P-68/P-69.1, `tmp/bluebook/P6a.txt`
        # lines 7810/8010/8608-8609): ethylarsane and trimethylbismuthane
        # cited directly; stibane's own stem is confirmed by the epic's
        # `bromodi(ethenyl)stibane` example (its unsaturated substituent
        # is separately out of scope, see below).
        ("CC[AsH2]", "ethylarsane"),
        ("C[Bi](C)C", "trimethylbismuthane"),
        ("CC[As](CC)CC", "triethylarsane"),
        ("CC[Sb](CC)CC", "triethylstibane"),
        ("C[As](C)C", "trimethylarsane"),
        ("C[Sb](C)C", "trimethylstibane"),
        # A halogen bonded directly to the metal alongside multiplied
        # alkyl substituents, mirroring `_group13_hydride.py`'s own
        # established convention.
        ("CC[As](CC)Cl", "chlorodi(ethyl)arsane"),
    ],
)
def test_group15_hydride_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bare_metal_hydride_name():
    assert smiles_to_iupac("[AsH3]") == "arsane"


def test_two_metal_atoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[As](CC)CC.CC[Sb](CC)CC")


def test_other_heteroatom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CO[As](C)C")


def test_unsaturated_substituent_raises():
    # The epic's own cited worked example, `BrSb(CH=CH2)2` ->
    # 'bromodi(ethenyl)stibane (PIN)', needs unsaturated-substituent
    # support this shared mechanism doesn't have yet -- a separate,
    # larger follow-up step, not this one.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Br[Sb](C=C)C=C")


def test_group13_hydride_unaffected():
    assert smiles_to_iupac("CC[Al](CC)CC") == "triethylalumane"


def test_group14_hydride_unaffected():
    assert smiles_to_iupac("CC[Ge](CC)(CC)CC") == "tetraethylgermane"


def test_phosphane_unaffected():
    assert smiles_to_iupac("CC[P](CC)CC") == "triethylphosphane"
