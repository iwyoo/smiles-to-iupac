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


def test_bicyclopropane():
    # PIN worked example `tmp/bluebook/P2.txt` ~7787: 1,1'-bi(cyclopropane)
    # -- parentheses needed here (unlike every other case below) to avoid
    # misreading as a von Baeyer 'bicyclo...' name.
    assert smiles_to_iupac("C1CC1C1CC1") == "1,1'-bi(cyclopropane)"


def test_bicyclopropane_with_halogen():
    assert smiles_to_iupac("C1CC1(Cl)C1CC1") == "1-chloro-1,1'-bi(cyclopropane)"


def test_bipyridine():
    # Real PubChem structure, CID 1474 (2,2'-bipyridine).
    assert smiles_to_iupac("C1=CC=NC(=C1)C2=CC=CC=N2") == "2,2'-bipyridine"


def test_bifuran():
    # Real PubChem structure, CID 80006 (2,2'-bifuran).
    assert smiles_to_iupac("C1=COC(=C1)C2=CC=CO2") == "2,2'-bifuran"


def test_bithiophene():
    # Real PubChem structure, CID 68120 (2,2'-bithiophene).
    assert smiles_to_iupac("C1=CSC(=C1)C2=CC=CS2") == "2,2'-bithiophene"


def test_biselenophene():
    # Real PubChem structure, CID 5141513 (2,2'-biselenophene).
    assert smiles_to_iupac("C1=C[Se]C(=C1)C2=CC=C[Se]2") == "2,2'-biselenophene"


def test_bipyrimidine():
    # Real PubChem structure, CID 123444 (2,2'-bipyrimidine).
    assert smiles_to_iupac("C1=CN=C(N=C1)C2=NC=CC=N2") == "2,2'-bipyrimidine"


def test_bipyridazine_lowest_attachment_locant():
    # Real PubChem structure, CID 12464241 (3,3'-bipyridazine) -- confirms
    # the attachment point is numbered 3 (the lower of pyridazine's two
    # symmetry-equivalent non-nitrogen alpha positions), not 6.
    assert smiles_to_iupac("C1=CC(=NN=C1)C2=NN=CC=C2") == "3,3'-bipyridazine"


def test_bipyrazine():
    # Real PubChem structure, CID 153669 (2,2'-bipyrazine).
    assert smiles_to_iupac("C1=CN=C(C=N1)C2=NC=CN=C2") == "2,2'-bipyrazine"


def test_halogen_on_bipyridine():
    assert smiles_to_iupac("c1ccc(Cl)c(-c2ccccn2)n1") == "3-chloro-2,2'-bipyridine"


def test_two_thiazole_rings_still_out_of_scope():
    # 1,3-thiazole is excluded by `_NON_NH_ROLE_SEQUENCES` (locant-prefixed
    # Hantzsch-Widman name, see `_ring_assembly_chain.py`'s own note).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1csc(-c2cscn2)n1")


def test_bipyrrole_carbon_attached_indicated_hydrogen():
    # Real PubChem structure, CID 260036 -- joined through carbon, so
    # each ring keeps its own N-H, cited at the front of the name (P-28.2.3).
    assert smiles_to_iupac("C1=CNC(=C1)C2=CC=CN2") == "1H,1'H-2,2'-bipyrrole"


def test_bipyrrole_nitrogen_attached_no_indicated_hydrogen():
    # Joined through the ring nitrogens themselves, so each ring's own N-H
    # position is occupied by the junction and needs no indicated-H prefix.
    assert smiles_to_iupac("c1ccn(-n2cccc2)c1") == "1,1'-bipyrrole"


def test_biimidazole_nitrogen_attached():
    # Real PubChem structure, CID 15034216. Both rings' own N-H position
    # is the junction, so no indicated hydrogen is needed (same mechanism
    # as bipyrrole above).
    assert smiles_to_iupac("c1cn(-n2ccnc2)cn1") == "1,1'-biimidazole"


def test_bipyrazole_nitrogen_attached():
    # Real PubChem structure, CID 21981271.
    assert smiles_to_iupac("c1cnn(-n2cccn2)c1") == "1,1'-bipyrazole"


def test_biimidazole_carbon_attached_still_out_of_scope():
    # A carbon-attached junction leaves a real prototropic-tautomer
    # ambiguity unresolved (which ring nitrogen is "N1"), unlike the
    # N-N-attached case above -- deferred, not yet supported.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1c[nH]c(-c2[nH]ccn2)n1")


def test_mixed_benzo_and_pyridine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc(-c2ccccn2)cc1")
