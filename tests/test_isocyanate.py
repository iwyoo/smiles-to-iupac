import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 8022.
        ("CCN=C=O", "isocyanatoethane"),
        # PubChem CID 61277.
        ("CC(C)N=C=O", "2-isocyanatopropane"),
    ],
)
def test_isocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mononuclear_case_follows_blue_book_prefix_rule():
    # P-61.8's own 'isocyanatoborane (PIN)' worked example shows the plain
    # 'isocyanato' prefix applies even to a very small parent hydride, so
    # methyl isocyanate is expected to be 'isocyanatomethane' by the same
    # mechanical rule -- even though PubChem's own auto-generated name for
    # this exact structure (CID 12228) is 'methylimino(oxo)methane', a
    # different, non-prefix parent selection. Treated as another instance
    # of the PubChem-autoname-vs-PIN mismatch documented elsewhere in this
    # project (see module docstring). This is also a regression check for
    # the "own carbon must be excluded from the chain search" fix (module
    # docstring) -- without it, this exact mononuclear case is the one that
    # would misfire.
    assert smiles_to_iupac("CN=C=O") == "isocyanatomethane"


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(N=C=O)CC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCN=C=O")


def test_two_isocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C=NCN=C=O")
