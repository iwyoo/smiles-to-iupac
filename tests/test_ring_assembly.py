import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_biphenyl():
    # 'biphenyl' with primed/unprimed ring-locant notation (P-28) is a
    # well-known retained name (PIN '1,1'-biphenyl'); '4,4'-dichloro-1,1'-
    # biphenyl' below is also a real, independently-named compound (a PCB
    # congener).
    assert smiles_to_iupac("c1ccc(cc1)-c1ccccc1") == "1,1'-biphenyl"


def test_biphenyl_alternate_smiles_order():
    assert smiles_to_iupac("c1ccc(-c2ccccc2)cc1") == "1,1'-biphenyl"


def test_4_chlorobiphenyl():
    assert smiles_to_iupac("Clc1ccc(cc1)-c1ccccc1") == "4-chloro-1,1'-biphenyl"


def test_4_4_dichlorobiphenyl():
    assert smiles_to_iupac("Clc1ccc(cc1)-c1ccc(Cl)cc1") == "4,4'-dichloro-1,1'-biphenyl"


def test_2_chlorobiphenyl_lowest_locant():
    assert smiles_to_iupac("Clc1ccccc1-c1ccccc1") == "2-chloro-1,1'-biphenyl"


def test_3_3_dichlorobiphenyl_unprimed_lower_than_primed():
    # Both rings carry one chlorine at the meta position: 3 (unprimed) is
    # preferred over 3' for whichever ring is cited first, so the locant set
    # is {3, 3'} either way, but the naming must still pick a consistent
    # unprimed/primed assignment.
    assert smiles_to_iupac("Clc1cccc(c1)-c1cccc(Cl)c1") == "3,3'-dichloro-1,1'-biphenyl"


def test_naphthalene_still_raises_not_ring_assembly():
    # Fused (shared-atom) two-ring aromatic must still go through
    # _aromatic.py, not be misdetected as a ring assembly.
    assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


def test_terphenyl_routes_to_ring_assembly_chain_module():
    # Three or more rings are now `_ring_assembly_chain.py`'s job (see
    # tests/test_ring_assembly_chain.py) -- this just confirms the N=2
    # module here no longer claims (and mis-fails on) the N=3 shape.
    assert smiles_to_iupac("c1ccc(cc1)-c1ccc(cc1)-c1ccccc1") == "11,21:24,31-terphenyl"
