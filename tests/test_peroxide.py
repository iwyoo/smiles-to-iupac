import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 'peroxy' substituent-prefix naming mirrors 'oxy' (_ether.py, see
        # test_ether.py: 'methoxymethane'/'1-methoxypropane', already
        # cross-checked against PubChem) with -O- replaced by -O-O-,
        # including the same P-14.3.4.2(b) omitted-locant convention for a
        # homogeneous two-carbon chain vs. the cited locant on a longer one.
        ("COOC", "methylperoxymethane"),
        ("CCOOCC", "ethylperoxyethane"),
        ("COOCC", "methylperoxyethane"),
        ("COOCCC", "1-methylperoxypropane"),
    ],
)
def test_peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_asymmetric_prefers_longer_chain_as_parent():
    assert smiles_to_iupac("CCCCOOC") == "1-methylperoxybutane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'peroxy' outside. Structure verified against PubChem: CID 22572410
    # ("1-propan-2-ylperoxybutane" -- PubChem's own PIN-style name, this
    # project keeps its usual CAS-style substituent name instead).
    assert smiles_to_iupac("CCCCOOC(C)C") == "1-(1-methylethyl)peroxybutane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOC1CCCCC1")


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=COOC")


def test_three_oxygens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOCOC")
