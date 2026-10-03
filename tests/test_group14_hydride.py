import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book worked example (P-69.1, `tmp/bluebook/P6a.txt` lines
        # 8602-8613): tetraethylplumbane cited directly.
        ("CC[Pb](CC)(CC)CC", "tetraethylplumbane"),
        # Real PubChem-confirmed structures spanning all three Group 14
        # elements (Ge/Sn/Pb), matching PubChem's own PIN-generated names.
        ("CC[Ge](CC)(CC)CC", "tetraethylgermane"),  # CID 11703
        ("CC[Sn](CC)(CC)CC", "tetraethylstannane"),  # CID 11704
        ("C[Ge](C)(C)C", "tetramethylgermane"),
        ("C[Sn](C)(C)C", "tetramethylstannane"),
        ("C[Pb](C)(C)C", "tetramethylplumbane"),
        # Plain phenyl substituents, mirroring `_group13_hydride.py`'s
        # equivalent cases.
        ("c1ccccc1[GeH3]", "phenylgermane"),
        ("c1ccccc1[SnH2]c1ccccc1", "diphenylstannane"),
        # A halogen bonded directly to the metal alongside multiplied
        # phenyl substituents: mirrors `_group13_hydride.py`'s own
        # established convention (`chlorodi(ethyl)alumane`), not
        # PubChem's own generated string for this shape (CID 12540 gives
        # 'chloro(triphenyl)stannane' -- this project deliberately puts
        # the multiplying prefix outside the parentheses instead, per the
        # Blue Book's 'ethyldi(methyl)phosphane (PIN)' worked example).
        ("Cl[Sn](c1ccccc1)(c1ccccc1)c1ccccc1", "chlorotri(phenyl)stannane"),
    ],
)
def test_group14_hydride_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bare_metal_hydride_name():
    assert smiles_to_iupac("[GeH4]") == "germane"


def test_two_metal_atoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Ge](C)(C)C.C[Sn](C)(C)C")


def test_other_heteroatom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CO[Ge](C)(C)C")


def test_group13_hydride_unaffected():
    assert smiles_to_iupac("CC[Al](CC)CC") == "triethylalumane"


def test_silicon_still_unsupported_is_named():
    assert smiles_to_iupac("C[Si](C)(C)C") == "tetramethylsilane"
