import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 11646: mononuclear parent -- also a regression check
        # for the "own carbon must be excluded from the chain search" fix
        # (module docstring); this exact case is the one that would
        # misfire without it.
        ("C[N+]#[C-]", "isocyanomethane"),
        # PubChem CID 12226.
        ("CC[N+]#[C-]", "isocyanoethane"),
    ],
)
def test_isocyanide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N+]#[C-])CC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N+]#[C-]")


def test_two_isocyanide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[C-]#[N+]C[N+]#[C-]")
