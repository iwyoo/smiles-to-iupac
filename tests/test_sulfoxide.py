import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-63.6's own worked example for ethyl butyl sulfoxide
        # (R=butyl, R'=ethyl): '1-(ethanesulfinyl)butane (PIN)'.
        ("CCCCS(=O)CC", "1-(ethanesulfinyl)butane"),
        # Symmetric case, mononuclear parent (P-14.3.4.2(a), no locant):
        # dimethyl sulfoxide (DMSO). Note this differs from PubChem's own
        # auto-generated (non-PIN) name 'methylsulfinylmethane' -- see
        # module docstring for why the acid-derived acyl prefix
        # 'methanesulfinyl' is used instead of a plain 'methylsulfinyl'.
        ("CS(=O)C", "(methanesulfinyl)methane"),
        # Symmetric two-carbon case (P-14.3.4.2(b): locant omitted despite
        # the enclosed prefix, mirroring the sulfone worked example in
        # `_sulfone.py`'s own test file exactly).
        ("CCS(=O)CC", "(ethanesulfinyl)ethane"),
    ],
)
def test_sulfoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)S(=O)C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S1CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CS(=O)C")


def test_two_sulfoxide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)CS(=O)C")
