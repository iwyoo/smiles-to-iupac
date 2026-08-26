import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # dispiro[4.1.4^7.2^5]tridecane (see tests/test_polyspiro.py) with
        # the single-atom internal arc (locant 6, bonded directly to both
        # spiro atoms, so its position is structurally forced either way --
        # no numbering choice for the heteroatom rule to even act on)
        # replaced by O. PubChem-verified for this exact SMILES
        # ('6-oxadispiro[4.1.4^7.2^5]tridecane', CID 45098582).
        ("C1CCC2(C1)CCC3(O2)CCCC3", "6-oxadispiro[4.1.4^7.2^5]tridecane"),
        # same position, N instead of O -- not independently found on
        # PubChem, but the O/N/S -> oxa/aza/thia prefix mapping itself is
        # already independently verified elsewhere (_von_baeyer_heteroatom.py,
        # _spiro_heteroatom.py); this only exercises that same mapping on
        # a position whose locant is structurally forced regardless of
        # element.
        ("C1CCC2(C1)CCC3(N2)CCCC3", "6-azadispiro[4.1.4^7.2^5]tridecane"),
        # a halogen substituent coexists with the ring heteroatom, the same
        # substituent-handling machinery (`_substituents_for_ring`) already
        # verified by every other single-heteroatom module in this project.
        ("ClC1CCCC12CCC1(CCCC1)O2", "1-chloro-6-oxadispiro[4.1.4^7.2^5]tridecane"),
    ],
)
def test_smiles_to_iupac_linear_polyspiro_heteroatom(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A terminal ring whose two non-spiro atoms are both bonded
        # directly to the spiro atom (a 3-membered terminal ring): a
        # genuine numbering-direction choice exists, so SP-1.8.1's "low
        # locants to heteroatoms" rule (see module docstring) must pick the
        # direction that gives O locant 1, not 2. PubChem's own computed
        # name for the identical structure (ConnectivitySMILES matching
        # CID 114820437) instead says '2-oxa...', contradicting its own
        # SP-1.8.1 text -- treated as a PubChem generator limitation for
        # this shape (see module docstring) and not followed; this
        # assertion instead reflects a direct, manually reviewed
        # application of the quoted rule rather than an independently
        # confirmed name string.
        ("C1CCC2(C1)CCC3(CC2)CO3", "1-oxadispiro[2.2.4^6.2^3]dodecane"),
        # same shape, both terminal rings 3-membered (a genuine tie for
        # which terminal ring is numbered first too) -- same reasoning;
        # PubChem's name for the identical structure (ConnectivitySMILES
        # matching CID 134451841) says '8-oxa...', which this module
        # deliberately does not follow (see above).
        ("C1CC12CCC3(CC2)CO3", "1-oxadispiro[2.2.2^6.2^3]decane"),
    ],
)
def test_smiles_to_iupac_linear_polyspiro_heteroatom_tie_break(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_heteroatoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2(C1)CCC3(O2)CCCN3")


def test_heteroatom_outside_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCC12CCC3(CC2)CCC3")


def test_unsupported_heteroatom_element_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2(C1)CCC3(P2)CCCC3")


def test_unsaturated_heteroatom_polyspiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2(C1)C=CC3(O2)CCCC3")
