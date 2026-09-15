import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the -SO3H).
        # PubChem CID 20376717.
        ("OS(=O)CS(=O)(=O)O", "sulfinomethanesulfonic acid"),
        # Two-carbon chain: not the P-14.3.4.2(b) no-locant case, since a
        # sulfino substituent is present in addition to the -SO3H suffix
        # (mirrors `_sulfonic_acid_thiol.py`'s identical -1- citation for
        # the same chain shape; PubChem's own computed name for this CID,
        # 13983537, omits the locant, but PubChem's algorithmic names are
        # not always strict PINs -- see that module's own test for the
        # same discrepancy with a thiol in place of the sulfino group).
        ("OS(=O)CCS(=O)(=O)O", "2-sulfinoethanesulfonic acid"),
        # PubChem CID 57310813.
        ("OS(=O)CCCS(=O)(=O)O", "3-sulfinopropane-1-sulfonic acid"),
        # Multiple sulfinic acids (multiplying prefix). PubChem CID
        # 123745263.
        ("OS(=O)CC(S(=O)O)CS(=O)(=O)O", "2,3-disulfinopropane-1-sulfonic acid"),
    ],
)
def test_sulfonic_acid_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_sulfonic_acid_still_routes_normally():
    # No sulfinic acid present -- must still reach `_sulfonic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_plain_sulfinic_acid_still_routes_normally():
    # No sulfonic acid present -- must still reach `_sulfinic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"


def test_multiple_sulfonic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CCS(=O)(=O)O")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(S(=O)O)S(=O)(=O)O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(S(=O)O)S(=O)(=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(S(=O)O)(CC1)S(=O)(=O)O")
