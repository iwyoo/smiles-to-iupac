import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PH3, already handled by `_phosphane.py`'s own mononuclear case
        # -- confirmed still unaffected by this module.
        ("P", "phosphane"),
        # Blue Book's own worked example (P6a.pdf, line 7709):
        # "diphosphane (preselected name)". PubChem's own auto-generated
        # name for this exact structure ("phosphanylphosphane", CID
        # 139283) disagrees with the actual PIN -- the primary text is
        # followed instead, same as `_silane_chain.py`'s own precedent.
        ("PP", "diphosphane"),
        ("PPP", "triphosphane"),
        # 'tetra' keeps its terminal 'a' before 'phosphane' (which starts
        # with a consonant, unlike '-ane'): tetraphosphane, not
        # tetrphosphane.
        ("PPPP", "tetraphosphane"),
    ],
)
def test_phosphane_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_phosphane_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("P(P)(P)P")


def test_cyclic_phosphane_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("P1PPPP1")


def test_carbon_phosphorus_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CPP")
