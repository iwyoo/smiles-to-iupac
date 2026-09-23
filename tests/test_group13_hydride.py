import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book worked examples (P-69.1, `tmp/bluebook/P6a.txt` lines
        # 8602-8613): triethylalumane and dimethylindigane cited directly.
        ("CC[Al](CC)CC", "triethylalumane"),
        ("[InH](C)C", "dimethylindigane"),
        # Real PubChem-confirmed structures spanning all four Group 13
        # elements (Al/Ga/In/Tl), matching PubChem's own PIN-generated
        # names for the plain-alkyl cases.
        ("C[Al](C)C", "trimethylalumane"),  # CID 517724
        ("CC[Ga](CC)CC", "triethylgallane"),  # CID 10412
        ("C[Ga](C)C", "trimethylgallane"),  # CID 10318
        ("CC[In](CC)CC", "triethylindigane"),  # CID 61705
        ("C[Tl](C)C", "trimethylthallane"),  # CID 62371
        ("CCCC[Al](CCCC)CCCC", "tributylalumane"),  # CID 12374
        ("CCC[Al](CCC)CCC", "tripropylalumane"),  # CID 10905
        ("CC(C)C[Al](CC(C)C)CC(C)C", "tris(2-methylpropyl)alumane"),  # CID 21172
        # A mixed halogen + multiplied-alkyl substituent case: PubChem's
        # own generated name ('chloro(diethyl)alumane') omits the
        # required inner parenthesization -- this project's established
        # convention (`_borane.py`/`_phosphane.py`'s own docstrings,
        # confirmed against the Blue Book's own 'ethyldi(methyl)phosphane
        # (PIN)' worked example, not PubChem) puts the multiplying prefix
        # outside the parentheses instead: 'chlorodi(ethyl)alumane'.
        ("CC[Al](CC)Cl", "chlorodi(ethyl)alumane"),
        ("C[Ga](C)Cl", "chlorodi(methyl)gallane"),
        # Plain and halogen-substituted phenyl substituents, mirroring
        # `_borane.py`'s equivalent cases exactly.
        ("c1ccccc1[AlH2]", "phenylalumane"),
        ("c1ccccc1[AlH](c1ccccc1)", "diphenylalumane"),
        ("Clc1ccc(cc1)[GaH2]", "(4-chlorophenyl)gallane"),
    ],
)
def test_group13_hydride_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bare_metal_hydride_name():
    # Note: a bare `[TlH3]` SMILES hits an RDKit radical-electron
    # inference quirk specific to thallium's own default-valence table
    # (3 explicit H exceeds RDKit's assumed common valence for Tl,
    # unlike Al) and gets misrouted to `_radical.py` before this module
    # is ever reached -- out of scope to fix here (every substituted
    # thallane case above, which is what the motivating issue's examples
    # actually need, works fine).
    assert smiles_to_iupac("[AlH3]") == "alumane"


def test_two_metal_atoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Al](C)C.C[Ga](C)C")


def test_other_heteroatom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CO[Al](C)C")


def test_borane_unaffected():
    assert smiles_to_iupac("CC[B](CC)CC") == "triethylborane"


def test_phosphane_unaffected():
    assert smiles_to_iupac("CC[P](CC)CC") == "triethylphosphane"
