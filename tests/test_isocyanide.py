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


def test_phenyl_isocyanide_direct_bond():
    # P-44.1.2.2 rule (1): 'isocyano' has no suffix form, so the ring is
    # always senior to a chain of the same class, mirroring `_nitro.py`'s/
    # `_azide.py`'s identical rule -- confirmed by PubChem CID 13606.
    assert smiles_to_iupac("c1ccccc1[N+]#[C-]") == "isocyanobenzene"


def test_phenyl_isocyanide_chain():
    # PubChem CID 82558/427710 give "isocyanomethylbenzene"/
    # "2-isocyanoethylbenzene" (no parentheses), but this codebase follows
    # `_nitro.py`'s/`_azide.py`'s identical, Blue-Book-verified rule
    # instead -- a locant-bearing compound substituent prefix is enclosed
    # regardless of which simple prefix it carries.
    assert smiles_to_iupac("c1ccccc1C[N+]#[C-]") == "(isocyanomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CC[N+]#[C-]") == "(2-isocyanoethyl)benzene"


def test_phenyl_isocyanide_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1C[N+]#[C-]")


def test_phenyl_isocyanide_unsaturation_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1C[N+]#[C-]")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N+]#[C-]")


def test_two_isocyanide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[C-]#[N+]C[N+]#[C-]")
