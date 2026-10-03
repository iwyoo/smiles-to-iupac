import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # SiH4, the silicon analogue of methane -- 'silane' itself, no
        # 'mono-' prefix, same pattern as the retained alkane names.
        ("[SiH4]", "silane"),
        # cross-checked against well-known compounds disilane (Si2H6) and
        # trisilane (Si3H8).
        ("[SiH3][SiH3]", "disilane"),
        ("[SiH3][SiH2][SiH3]", "trisilane"),
        # 'tetra'/'penta' keep their terminal 'a' before 'silane' (which
        # starts with a consonant, unlike '-ane'): tetrasilane, not
        # tetrsilane.
        ("[SiH3][SiH2][SiH2][SiH3]", "tetrasilane"),
        ("[SiH3][SiH2][SiH2][SiH2][SiH3]", "pentasilane"),
    ],
)
def test_smiles_to_iupac_silane_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_silane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SiH3][Si]([SiH3])([SiH3])[SiH3]")


def test_cyclic_silane():
    assert smiles_to_iupac("[SiH2]1[SiH2][SiH2][SiH2][SiH2]1") == "pentasilolane"


def test_carbon_silicon_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[SiH2][SiH3]")
