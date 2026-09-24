import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_terpyridine():
    # Real PubChem structure, CID 70848 (2,2':6',2''-terpyridine).
    assert smiles_to_iupac("C1=CC=NC(=C1)C2=NC(=CC=C2)C3=CC=CC=N3") == "12,22:26,32-terpyridine"


def test_quaterpyridine_matches_primary_source_pin_worked_example():
    # Real PubChem structure, CID 10876425 -- exact PIN worked example
    # match (tmp/bluebook/P2.txt ~7918: "12,22:26,32:36,42-quaterpyridine").
    assert (
        smiles_to_iupac("C1=CC=NC(=C1)C2=NC(=CC=C2)C3=CC=CC(=N3)C4=CC=CC=N4")
        == "12,22:26,32:36,42-quaterpyridine"
    )


def test_halogen_substituent():
    assert (
        smiles_to_iupac("C1=CC=NC(=C1Cl)C2=NC(=CC=C2)C3=CC=CC=N3")
        == "13-chloro-12,22:26,32-terpyridine"
    )


def test_two_pyridine_rings_still_out_of_scope():
    # N=2 is out of scope here -- `_ring_assembly.py` (benzo-only) doesn't
    # claim it either, a separate gap not part of this milestone.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc(-c2ccccn2)nc1")


def test_other_heteroaromatic_ring_still_out_of_scope():
    # furan (a different role-sequence/symmetry) is a separate follow-up,
    # not attempted here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc(-c2ccc(-c3ccco3)o2)o1")


def test_mixed_pyridine_and_benzo_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc(-c2ccc(-c3ccccn3)cc2)cc1")


def test_benzo_and_cycloalkane_chains_unaffected():
    assert (
        smiles_to_iupac("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1") == "11,21:24,31-terphenyl"
    )
    assert smiles_to_iupac("C1CC1C1CC1C1CC1") == "11,21:22,31-tercyclopropane"
