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


def test_two_pyridine_rings_routes_to_ring_assembly_module():
    # N=2 is `_ring_assembly.py`'s own job (P-28.2.1, not this module's
    # P-28.3 composite-locant scheme) -- confirms the N=3-6 module here no
    # longer claims (and mis-fails on) the N=2 shape.
    assert smiles_to_iupac("c1ccc(-c2ccccn2)nc1") == "2,2'-bipyridine"


def test_locant_prefixed_hantzsch_widman_ring_still_out_of_scope():
    # 1,3-thiazole matches a _ROLE_SEQUENCES parent but is excluded from
    # `_NON_NH_ROLE_SEQUENCES` (see _ring_assembly_chain.py's own note)
    # since its name starts with a locant and P-28.3.1 has no confirmed
    # worked example for how 'ter' composes with such a name.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1csc(-c2csc(-c3cscn3)n2)n1")


def test_mixed_pyridine_and_benzo():
    assert smiles_to_iupac("c1ccc(-c2ccc(-c3ccccn3)cc2)cc1") == "2-([1,1'-biphenyl]-4-yl)pyridine"


def test_benzo_and_cycloalkane_chains_unaffected():
    assert (
        smiles_to_iupac("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1") == "11,21:24,31-terphenyl"
    )
    assert smiles_to_iupac("C1CC1C1CC1C1CC1") == "11,21:22,31-tercyclopropane"
