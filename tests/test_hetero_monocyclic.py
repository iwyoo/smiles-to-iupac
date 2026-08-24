import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # oxa series (P-22.2.1 Table 22.1), formulas cross-checked:
        # oxirane C2H4O, oxetane C3H6O, oxolane C4H8O, oxane C5H10O,
        # oxepane C6H12O.
        ("C1CO1", "oxirane"),
        ("C1CCO1", "oxetane"),
        ("C1CCCO1", "oxolane"),
        ("C1CCCCO1", "oxane"),
        ("C1CCCCCO1", "oxepane"),
        # thia series, formulas cross-checked: thiirane C2H4S, thietane
        # C3H6S, thiolane C4H8S, thiane C5H10S, thiepane C6H12S.
        ("C1CS1", "thiirane"),
        ("C1CCS1", "thietane"),
        ("C1CCCS1", "thiolane"),
        ("C1CCCCS1", "thiane"),
        ("C1CCCCCS1", "thiepane"),
        # aza series: aziridine/azetidine coincide with the literal
        # Hantzsch-Widman stem construction; pyrrolidine/piperidine are
        # retained names that are themselves the PIN (not
        # azolidine/azinane) per P-22.2.1 Table 2.3; azepane is the stem
        # form again. Formulas cross-checked: aziridine C2H5N, azetidine
        # C3H7N, pyrrolidine C4H9N, piperidine C5H11N, azepane C6H13N.
        ("C1CN1", "aziridine"),
        ("C1CCN1", "azetidine"),
        ("C1CCCN1", "pyrrolidine"),
        ("C1CCCCN1", "piperidine"),
        ("C1CCCCCN1", "azepane"),
    ],
)
def test_smiles_to_iupac_hetero_monocyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_hetero_monocyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCCCO1")


def test_ring_size_outside_scope_raises():
    # 8-membered ring: out of scope (this module covers 3-7).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCCCO1")


def test_two_heteroatoms_raises():
    # 1,4-dioxane: two ring heteroatoms, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1COCCO1")


def test_unsupported_heteroatom_element_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCP1")
