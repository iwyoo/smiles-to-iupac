import pytest
from rdkit import Chem
from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._fullerene import (
    _C84_D6H_ISOMER_24_SMILES,
    _FULLERENE_C60_SMILES,
    _FULLERENE_C70_SMILES,
    _FULLERENE_C76_SMILES,
)
from smiles_to_iupac._fullerene_spiral import (
    _C84_IPR_ISOMERS,
    _canonical_spiral,
    _pentagon_positions,
    _planar_faces,
)
from smiles_to_iupac._ring_system_seniority import ring_seniority_key


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC12CCC#CCC2", "spiro[2.6]non-6-yne"),
        ("C1CCC2(CC1)CC#CC=C2", "spiro[5.5]undec-1-en-3-yne"),
    ],
)
def test_spiro_triple_bond_and_ene_yne_endings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_monospiro_still_resolves_via_spiro_module():
    assert smiles_to_iupac("C1CCCC12CCCCC2") == "spiro[4.5]decane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC2CC#CC1C2", "bicyclo[3.2.1]oct-2-yne"),
    ],
)
def test_bicyclic_unsaturated_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_epoxycholestane():
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC4C5(C3(CCCC5)C)O4)C") == "5,6-epoxycholestane"
    )


def test_ring_compound_substituent():
    assert smiles_to_iupac("CC(CC)C1CCCCC1") == "(butan-2-yl)cyclohexane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@@H]1CCCC[C@@H]1C", "(1R,2S)-1,2-dimethylcyclohexane"),
    ],
)
def test_smiles_to_iupac_ring_cis_trans(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C@@H]1CCCC[C@H]1CC", "(1R,2R)-1-ethyl-2-methylcyclohexane", id="asymmetric_1_2_disubstituted_stereo"),
        pytest.param("C[C@@H]1C[C@H](C)CCC1", "(1R,3S)-1,3-dimethylcyclohexane", id="1_3_disubstituted_stereo"),
        pytest.param("C[C@H]1CCCCC1CC", "(2S)-1-ethyl-2-methylcyclohexane", id="single_ring_stereocenter_cites_the_specified_elements"),
    ],
)
def test_asymmetric_1_2_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CC=CCC1", "cyclohexa-1,3-diene"),
        ("BrC1C=CCCC1", "3-bromocyclohex-1-ene"),
        ("C=1=C=C=C=C=C=C=C=C=C=C1", "cycloundecaundecaene"),
        ("S=1=S=C=C=C=C=C=C=C=C=C1", "1λ4,2λ4-dithiacycloundecaundecaene"),
        ("C1=C=C=C=C=C=C=C=C=C=[S]=1", "1λ4-thiacycloundecaundecaene"),
        ("C1#CC(C)=CC=C1", "3-methyl-1,2-didehydrobenzene"),
    ],
)
def test_unsubstituted_cyclic_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=C1C=CC(=C)C=C1", "3,6-dimethylidenecyclohexa-1,4-diene"),
        ("C=C1C=CC=CC1=C", "5,6-dimethylidenecyclohexa-1,3-diene"),
    ],
)
def test_quinodimethane_is_a_cyclohexadiene_not_a_benzene_with_two_exocyclic_double_bonds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_enyne_ring():
    assert smiles_to_iupac("C1=CCC#CCCCCCCCCCC1") == "cyclopentadec-1-en-4-yne"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC/C=C/CC1", "(E)-cyclooctene"),
        ("OC1CCC/C=C/CC1", "(4E)-cyclooct-4-en-1-ol"),
        ("OC1CCCCC/C=C\\1", "(2Z)-cyclooct-2-en-1-ol"),
        ("C1CCC/C=C\\CC1C(=O)O", "(3Z)-cyclooct-3-ene-1-carboxylic acid"),
        ("O=C(OC)C1CCC/C=C/CC1", "methyl (4E)-cyclooct-4-ene-1-carboxylate"),
        ("OC1CCC/C=C/CC1O", "(4E)-cyclooct-4-ene-1,2-diol"),
        ("O[C@H]1CCC/C=C/CC[C@@H]1C", "(1S,2S,5E)-2-methylcyclonon-5-en-1-ol"),
        ("OC1C/C=C/CCC/C=C/1", "(2E,7E)-cyclonona-2,7-dien-1-ol"),
        ("C1=CCC[C@H]2O[C@H]2CC1", "(1R,8S)-9-oxabicyclo[6.1.0]non-4-ene"),
    ],
)
def test_ring_double_bond_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1cc2cc(c1)COC2", "3-oxabicyclo[3.3.1]nona-1(9),5,7-triene", id="heteroatom_inside_the_bridge"),
        pytest.param("c1cc2cc(c1)CCO2", "2-oxabicyclo[3.3.1]nona-1(9),5,7-triene", id="other_kekule_structure_gives_fewer_compound_locants"),
        pytest.param("c1cc2cc(c1)OCO2", "2,4-dioxabicyclo[3.3.1]nona-1(9),5,7-triene", id="two_heteroatoms_of_one_element"),
        pytest.param("c1cc2cc(c1)CNCC2", "3-azabicyclo[4.3.1]deca-1(10),6,8-triene", id="four_atom_bridge"),
        pytest.param("c1cc2cc(c1)CSCO2", "2-oxa-4-thiabicyclo[4.3.1]deca-1(10),6,8-triene", id="two_heteroatoms_of_two_elements"),
    ],
)
def test_meta_bridged_benzene_with_a_skeletal_heteroatom_is_a_von_baeyer_triene(smiles, expected):
    # P-31.1.4.2, P-23.3.1
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC/C=C/C=C\\CCCC1", "(1Z,3E)-cycloundeca-1,3-diene"),
        ("C1CC/C=C/CCC#CCCCC1", "(1E)-cyclotridec-1-en-5-yne"),
        ("C1=CC=CC=CC=CC=CC=CC=C1", "cyclotetradeca-1,3,5,7,9,11,13-heptaene"),
        ("C1=C\\C=C/C=C\\C=C/C=C\\1", "(1Z,3Z,5Z,7Z,9E)-cyclodeca-1,3,5,7,9-pentaene"),
        pytest.param("C1CC/C=C\\CC=CCCCC1", "(1Z)-cyclododeca-1,4-diene", id="unspecified_ring_double_bond_gets_no_descriptor"),
        pytest.param("C1=CCCCC/C=C/CC1", "(1E)-cyclodeca-1,5-diene", id="only_the_specified_bond_is_cited"),
    ],
)
def test_large_unsaturated_carbocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H]1CCC=CC1", "(4S)-4-ethylcyclohex-1-ene"),
        ("CC[C@H](C)[C@H]1CCC=CC1", "(4S)-4-[(2S)-butan-2-yl]cyclohex-1-ene"),
    ],
)
def test_ring_tetrahedral_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[C@@H]1[C@H]2CCCC[C@@H]12", "(1R,6S,7r)-7-methylbicyclo[4.1.0]heptane", id="cis_fusion_makes_the_substituted_carbon_pseudoasymmetric"),
        pytest.param("C[C@@H]1[C@H]2CCCC[C@H]12", "(1R,6R)-7-methylbicyclo[4.1.0]heptane", id="trans_fusion_leaves_the_substituted_carbon_without_descriptor"),
        pytest.param("C[C@@H]1[C@@H]2CCC[C@@H]12", "(1S,5S)-6-methylbicyclo[3.1.0]hexane", id="trans_fused_bicyclo_hexane"),
        pytest.param("C[C@@H]1[C@H]2CC[C@H]12", "(1R,4R)-5-methylbicyclo[2.1.0]pentane", id="trans_fused_bicyclo_pentane"),
        pytest.param("OC[C@@H]1[C@H]2CCCC[C@H]12", "[(1S,6S)-bicyclo[4.1.0]heptan-7-yl]methanol", id="trans_fused_ring_as_a_substituent"),
        pytest.param("CC1(C)[C@@H]2CCCC[C@H]21", "(1R,6R)-7,7-dimethylbicyclo[4.1.0]heptane", id="identical_groups_on_the_cyclopropane_carbon"),
    ],
)
def test_bridgehead_descriptors_of_a_bicyclo_n_1_0_alkane_with_one_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_stereocenter_with_two_substituents():
    assert smiles_to_iupac("CC[C@H](C)[C@H]1CC(C)C=CC1") == "(5R)-5-[(2S)-butan-2-yl]-3-methylcyclohex-1-ene"


def test_ortho_fused_tetrabenzenacyclooctaphane_is_named_with_a_warning():
    with pytest.warns(NonPreferredNameWarning, match="phane"):
        name = smiles_to_iupac("Cc1cccc2c1CC1=CC=CC(=C1)CC1=CC=CC(=C1)CC1=CC=CC(=C1)C2")
    assert name == "1³-methyl-1(1,2),3,5,7(1,3)-tetrabenzenacyclooctaphane"


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "c1ncc2cc1CCCc1cc3cc(c1)-c1cc(cc(c1)CC2)CCc1ccc(nc1)CCC3",
            "4(5,2),12(3,5)-dipyridina-1,8(1,3,5)-dibenzenabicyclo[6.6.0]tetradecaphane",
        ),
        (
            "c1cc2cc(c1)Cc1cc(c3c4ccc(c3c1)CCc1cc3ccccc3c3cc(ccc13)CC4)CCCCC2",
            "3(3,10)-phenanthrena-6(8,5,3,1)-naphthalena-8(1,3)-benzenaspiro[5.7]tridecaphane",
        ),
    ],
)
def test_phane_von_baeyer_and_spiro_skeletons(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_metacyclophane_is_named():
    assert (
        smiles_to_iupac("CC1CC2=CC(=CC=C2)CCC3=CC=CC1=C3")
        == "2-methyl-1,4(1,3)-dibenzenacyclohexaphane"
    )


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Oc1cc2ccc1CCc1ccc(cc1)CC2", "1,4(1,4)-dibenzenacyclohexaphan-1²-ol"),
    ],
)
def test_phane_with_principal_group_takes_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_separate_paracyclophane_units_not_misread_as_one_n4_phane():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2.C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2")


_C20_DODECAHEDRON = "c12c3c4c5c1c1c6c2c2c3c3c4c4c5c1c1c6c2c3c41"
_C24_TRUNCATED_OCTAHEDRON = "c12c3c4c1c1c5c6c2c2c7c3c3c8c4c4c1c5c1c(c26)c(c73)c1c84"
_C24_D6D = "c12c3c4c5c6c7c8c9c%10c%11c(c1c1c3c5c7c9c%111)c1c2c4c6c8c%101"
_C28_TD = "c12c3c4c5c1c1c6c2c2c3c3c7c2c2c8c7c7c3c4c3c5c4c(c37)c8c(c14)c62"
_C30_D5H = "c12c3c4c5c1c1c6c2c2c3c3c4c4c5c1c1c5c7c8c9c%10c7c1c4c%10c3c9c2c8c56"
_C30_C2V = "c12c3c4c5c1c1c6c7c8c9c%10c%11c%12c9c7c1c1c5c5c(c%11c7c%10c8c(c26)c3c7c45)c%121"
_C34_CS = "c12c3c4c5c1c1c6c7c8c9c%10c7c7c%11c%10c%10c9c9c8c6c5c5c4c4c6c3c(c2c71)c%11c6c%10c4c59"
_C36_D3H = "c12c3c4c5c1c1c6c7c8c9c(c75)c5c4c4c3c3c2c1c1c2c7c%10c(c23)c4c2c%10c3c7c1c6c8c3c9c52"
_C60_WITH_A_HEPTAGON = (
    "c12c3c4c5c1c1c6c7c8c9c%10c%11c%12c%13c%14c%15c%16c%17c%18c%15c%15c%19c%20c%21c(c%14%19)c%12c%12c%11c9c7c7c6c2c2c3c3c6c(c%20c9c3c4c3c5c4c(c%18c3c%159)c%17c(c%10c%16%13)c8c14)c%21c%12c7c26"
)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(_C20_DODECAHEDRON, "(C20-Ih)[5]fullerene", id="icosahedral_cage_of_pentagons_only"),
        pytest.param(_FULLERENE_C60_SMILES, "(C60-Ih)[5,6]fullerene", id="cubic_family_has_no_principal_axis"),
        pytest.param(_C24_TRUNCATED_OCTAHEDRON, "(C24-Oh)[4,6]fullerene", id="four_membered_rings_in_the_bracket"),
        pytest.param(_C28_TD, "(C28-Td)[5,6]fullerene", id="tetrahedral_with_mirror_planes"),
        pytest.param(_FULLERENE_C76_SMILES, "(C76-D2)[5,6]fullerene", id="dihedral_without_a_principal_axis"),
        pytest.param(_FULLERENE_C70_SMILES, "(C70-D5h(6))[5,6]fullerene", id="axial_pentagon_surrounded_by_hexagons"),
        pytest.param(_C30_D5H, "(C30-D5h(5))[5,6]fullerene", id="axial_pentagon_surrounded_by_pentagons"),
        pytest.param(_C24_D6D, "(C24-D6d(5))[5,6]fullerene", id="dihedral_with_diagonal_mirror_planes"),
        pytest.param(_C30_C2V, "(C30-C2v)[5,6]fullerene", id="cyclic_with_vertical_mirror_planes"),
        pytest.param(_C34_CS, "(C34-Cs)[5,6]fullerene", id="single_mirror_plane"),
        pytest.param(_C60_WITH_A_HEPTAGON, "(C60-C1)[5,6,7]fullerene", id="seven_membered_ring_without_symmetry"),
    ],
)
def test_fullerene_cage_named_from_its_symmetry_and_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_fullerene_whose_axis_ring_has_mixed_neighbours_is_flagged():
    with pytest.warns(NonPreferredNameWarning, match="ring-surround"):
        assert smiles_to_iupac(_C36_D3H) == "(C36-D3h)[5,6]fullerene"


def test_c84_isomer_24_via_spiral_match():
    assert smiles_to_iupac(_C84_D6H_ISOMER_24_SMILES) == "(C84-D6h(24))[5,6]fullerene"


def _pentagons_of(smiles):
    mol = Chem.MolFromSmiles(smiles)
    faces = _planar_faces(mol)
    return _pentagon_positions(_canonical_spiral(faces))


def test_canonical_spiral_matches_published_c60_code():
    assert _pentagons_of(_FULLERENE_C60_SMILES) == (1, 7, 9, 11, 13, 15, 18, 20, 22, 24, 26, 32)


def test_c84_isomer_table_has_24_distinct_entries():
    assert len(_C84_IPR_ISOMERS) == 24
    indices = sorted(index for index, _point_group in _C84_IPR_ISOMERS.values())
    assert indices == list(range(1, 25))
    for pentagons in _C84_IPR_ISOMERS:
        assert len(pentagons) == 12
        assert all(1 <= p <= 44 for p in pentagons)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1CCC2(CC1)CCC2=O", "7-hydroxyspiro[3.5]nonan-1-one"),
    ],
)
def test_hydroxy_monospiro_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1C2CC3CC(C2)(CC1C3=O)O", "5-hydroxyadamantan-2-one"),
    ],
)
def test_hydroxy_von_baeyer_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_biphenyl_assembly_stays_assembly():
    assert smiles_to_iupac("c1ccc(-c2ccc(-c3ccccn3)cc2)cc1") == "2-([1,1'-biphenyl]-4-yl)pyridine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)C=C1c2cccc(n2)Cc2cccc(n2)Cc2cccc(n2)Cc2cccc1n2",
            "[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-2-ylidene]acetic acid",
        ),
        pytest.param("OC(=O)C=C1C=C2Cc3cccc(n3)Cc3cccc(n3)Cc3cccc(n3)CC(=C1)N2",
            "[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-1⁴(1¹H)-ylidene]acetic acid", marks=pytest.mark.slow),
    ],
)
def test_pyridine_amplificants_and_ylidene_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_methyl_on_the_bridge_atom_of_a_bridged_fused_system():
    assert smiles_to_iupac("C1=CC2C(C)C1c1ccccc12") == "9-methyl-1,4-dihydro-1,4-methanonaphthalene"


def test_fused_system_to_be_bridged_has_the_most_skeletal_atoms_and_a_saturated_name_cites_no_hydro_locants():
    assert smiles_to_iupac("C1CCCC2C1CCCC1CCCC(CCC2)C1") == "hexadecahydro-1H-8,12-methanobenzo[13]annulene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC1C2CC3CC1C1C3C21", "octahydro-2,3-methanocyclopropa[cd]pentalen-1-ol", id="suffix_on_a_saturated_cage"),
        pytest.param("CC1C2CC3CC1C1C3C21", "1-methyloctahydro-2,3-methanocyclopropa[cd]pentalene", id="substituent_locant_of_the_cage"),
        pytest.param(
            "CC1C2C=CC3CC1C1C3C21",
            "1-methyl-2,2a,2a1,2b,3,5a-hexahydro-1H-2,3-methanocyclopropa[cd]indene",
            id="divalent_bridge_kept_beside_a_substituent",
        ),
        pytest.param(
            "CC1CC2CCC1c1ccccc12",
            "2-methyl-1,2,3,4-tetrahydro-1,4-ethanonaphthalene",
            id="prefix_on_a_benzo_bridged_system",
        ),
        pytest.param(
            "C=C1CC2CC3CC1CC2C3",
            "6-methylideneoctahydro-1H-2,5-methanoindene",
            id="ylidene_prefix_on_a_saturated_bridged_system",
        ),
    ],
)
def test_saturated_bridged_cage_is_numbered_by_the_bridged_fusion_rules(smiles, expected):
    # P-25.4.3: the bridge atom belongs to the parent name, so it takes no hydro prefix or added hydrogen
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC12CCC3(CC2)CCC3", "dispiro[2.2.3^6.2^3]undecane"),
        ("C1CCCC12CC3(CC2)CCCC3", "dispiro[4.1.4^7.2^5]tridecane"),
    ],
)
def test_smiles_to_iupac_linear_polyspiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_linear_polyspiro():
    assert smiles_to_iupac("C1=CC12CCC3(CC2)CCC3") == "dispiro[2.2.3^6.2^3]undec-1-ene"


def test_monospiro_two_rings_is_unaffected():
    assert smiles_to_iupac("C1CCCC12CCCCC2") == "spiro[4.5]decane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C12(CC2)CCC3(CC3)CCC4(CC4)CC1", "trispiro[2.2.2^6.2.2^11.2^3]pentadecane"),
    ],
)
def test_smiles_to_iupac_branched_polyspiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC2(C1)CCC3(CC2)CO3", "1-oxadispiro[2.2.4^6.2^3]dodecane"),
    ],
)
def test_smiles_to_iupac_linear_polyspiro_heteroatom_tie_break(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_heteroatom_polyspiro():
    assert smiles_to_iupac("C1CCC2(C1)C=CC3(O2)CCCC3") == "6-oxadispiro[4.1.4^7.2^5]tridec-12-ene"


def test_oxadispiro_off_spiro_stereocenter():
    assert (
        smiles_to_iupac("Cl[C@@H]1CC2(CCCC2)OC12CCCC2")
        == "(12R)-12-chloro-6-oxadispiro[4.1.4^7.2^5]tridecane"
    )


@pytest.mark.slow
def test_branched_polyspiro_off_spiro_stereocenter():
    assert (
        smiles_to_iupac("Cl[C@H]1CCC12CCC1(CC1)CCC1(CC1)CC2")
        == "(12S)-12-chlorotrispiro[2.2.2^6.2.3^11.2^3]hexadecane"
    )


def test_substituent_branch_stereocenter_on_polyspiro():
    assert (
        smiles_to_iupac("F[C@@H](Cl)C1CCC12CCC1(CCC1)CC2")
        == "1-[(S)-chloro(fluoro)methyl]dispiro[3.2.3^7.2^4]dodecane"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)C1CC2(C1)CC1(C2)CCCC1", "dispiro[3.1.4^6.1^4]undecane-2-carboxylic acid"),
        ("OC1CCC2(CC1)CC1(C2)CCCC1", "dispiro[4.1.5^7.1^5]tridecan-10-ol"),
        ("O=C1CC2(C1)CC1(C2)CCCC1", "dispiro[3.1.4^6.1^4]undecan-2-one"),
        ("CC1CC2(C1)CC1(C2)CCCC1C", "2,7-dimethyldispiro[3.1.4^6.1^4]undecane"),
        ("NC1CC2(CC1)CC1(C2)CCOC1", "2-oxadispiro[4.1.4^7.1^5]dodecan-9-amine"),
        ("OC1C2(CCC2)CC2(CCC2)CC13CCC3", "trispiro[3.1.3^6.1.3^11.1^4]pentadecan-5-ol"),
    ],
)
def test_polyspiro_ring_bearing_a_suffix_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_fused_bicyclic_is_not_spiro():
    assert smiles_to_iupac("C1CCC2CCCCC2C1") == "decahydronaphthalene"


def test_spiro_hydroxyl_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("OCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanol"


def test_spiro_amine_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("NCC1CCCC12CCCCC2") == "1-(spiro[4.5]decan-1-yl)methanamine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[CH2+]C1CCC2(CC1)CCCC2[CH2+]", "(spiro[4.5]decane-1,8-diyl)bis(methylium)", id="spiro_linker_between_two_carbenium_units"),
    ],
)
def test_spiro_linker_between_cationic_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C@@H]1CCCC12CCOCC2", "(1R)-1-methyl-8-oxaspiro[4.5]decane", id="oxaspiro_other_configuration"),
        pytest.param("C1CCCC12CCOCC2", "8-oxaspiro[4.5]decane", id="plain_oxaspiro_without_stereo_still_resolves"),
    ],
)
def test_oxaspiro_other_configuration_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O1CCC[C@]12OCCC2", "(5R)-1,6-dioxaspiro[4.4]nonane", id="dioxa_one_enantiomer"),
        pytest.param("O1CCC[C@@]12OCCC2", "(5S)-1,6-dioxaspiro[4.4]nonane", id="dioxa_other_enantiomer"),
        pytest.param("O=C1CCN[C@@]12NCCC2=O", "(5R)-1,6-diazaspiro[4.4]nonane-4,9-dione", id="diaza_dione"),
        pytest.param("O=C1CCC[C@]12CCCC2=O", "(5S)-spiro[4.4]nonane-1,6-dione", id="carbocyclic_dione"),
        pytest.param("O=C1CCC[C@@]12CCCC2=O", "(5R)-spiro[4.4]nonane-1,6-dione", id="carbocyclic_dione_other_enantiomer"),
        pytest.param("C1CC[C@]2(C1)CCCC2", "spiro[4.4]nonane", id="no_descriptor_without_stereogenicity"),
    ],
)
def test_spiro_atoms_of_the_type_xabab_take_a_descriptor_from_the_ring_pairing(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_radical_centres_on_substituents_of_a_spiro_ring():
    assert smiles_to_iupac("[CH2]C1CCC2(CC1)CCCC2[CH2]") == "(spiro[4.5]decane-1,8-diyl)dimethyl"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2]C1CCC2(CC1)CCCC2", "(spiro[4.5]decan-8-yl)methyl", id="spiro_radical_on_substituent_branch"),
        pytest.param("[CH]1CCC2(C=CC2)CC1", "spiro[3.5]non-1-en-7-yl", id="unsaturated_spiro_radical"),
    ],
)
def test_spiro_radical_as_substituent_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH]1CCC2(CC1)CCC[C@@H](C)C2", "(8R)-8-methylspiro[5.5]undecan-3-yl", id="stereo_spiro_radical_resolves"),
        pytest.param("C[C@H]1CCC[C@]2(C1)CCCC2", "(7S)-7-methylspiro[4.5]decane", id="single_off_spiro_stereocenter"),
    ],
)
def test_stereo_spiro_radical_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_spiro_thiol_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("SCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanethiol"


def test_bicyclic_is_not_tetracyclic():
    assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1C23CC12C3", "tricyclo[1.1.1.0^1,3]pentane"),
    ],
)
def test_smiles_to_iupac_propellane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_myrtenol_is_a_ring_prefix_on_methanol():
    assert (
        smiles_to_iupac("CC1(C2CC=C(C1C2)CO)C")
        == "(6,6-dimethylbicyclo[3.1.1]hept-2-en-2-yl)methanol"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NC12CC3CC(CC(C3)C1)C2", "adamantan-1-amine"),
    ],
)
def test_von_baeyer_amine_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_amine_ring_unsaturation():
    assert smiles_to_iupac("C1C2CC(C1C=C2)N") == "bicyclo[2.2.1]hept-5-en-2-amine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[C+]12CC3CC(CC(C3)C1)C2", "adamantan-1-ylium"),
    ],
)
def test_von_baeyer_carbenium_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[CH2+]C1CC2CCC1C2[CH2+]", "(bicyclo[2.2.1]heptane-2,7-diyl)bis(methylium)", id="bridged_linker_between_two_carbenium_units"),
    ],
)
def test_bridged_linker_between_cationic_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # same skeleton, nitrogen instead of oxygen.
        ("C1CC2CCC1N2", "7-azabicyclo[2.2.1]heptane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_tetraaza_multiplying_prefix_elides_vowel():
    assert (
        smiles_to_iupac("CC1CN2CCNCCNCCN(CC2)C1")
        == "12-methyl-1,4,7,10-tetrazabicyclo[8.3.2]pentadecane"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC2CCN1O2", "7-oxa-1-azabicyclo[2.2.1]heptane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom_mixed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O1C2CC3CC1CC(Cl)(C2)C3", "5-chloro-2-oxaadamantane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom_tricyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O1C2CC3CC1CC(C2)N3", "2-oxa-6-azaadamantane", id="different_ring_heteroatoms_tricyclic"),
        pytest.param("O1C2CC3OC1CC(C2)C3", "2,4-dioxaadamantane", id="same_element_ring_heteroatoms_tricyclic"),
    ],
)
def test_two_different_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_von_baeyer_ketone_now_resolves():
    assert smiles_to_iupac("O=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-en-2-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[C]12CC3CC(CC(C3)C1)C2", "adamantan-1-yl"),
    ],
)
def test_von_baeyer_radical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_radical_on_substituent_branch():
    assert smiles_to_iupac("[CH2]C1CC2CCC1C2") == "(bicyclo[2.2.1]heptan-2-yl)methyl"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[SeH]C1C2CC3CC1CC(C2)C3", "adamantane-2-selenol"),
    ],
)
def test_von_baeyer_selenol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[SeH]C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-selenol"),
    ],
)
def test_von_baeyer_spiro_selenol_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Se]=C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-selone"),
    ],
)
def test_von_baeyer_spiro_selone_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_selone_ring_unsaturation():
    assert smiles_to_iupac("[Se]=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-selone"


def test_von_baeyer_sulfinamide_name():
    assert smiles_to_iupac("NS(=O)C1CC2CCC1C2") == "bicyclo[2.2.1]heptane-2-sulfinamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NS(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfinamide"),
    ],
)
def test_von_baeyer_spiro_sulfinamide_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CNS(=O)C1CC2CCC1C2", "N-methylbicyclo[2.2.1]heptane-2-sulfinamide"),
        ("NS(=O)C1CCCC2(C1)C=CCCC2", "spiro[5.5]undec-7-ene-2-sulfinamide"),
    ],
)
def test_sulfinamide_on_polycyclic_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OS(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfinic acid"),
    ],
)
def test_von_baeyer_spiro_sulfinic_acid_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_sulfinic_acid_ring_unsaturation():
    assert smiles_to_iupac("OS(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfinic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NS(=O)(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfonamide"),
    ],
)
def test_von_baeyer_spiro_sulfonamide_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_sulfonamide_ring_unsaturation():
    assert smiles_to_iupac("NS(=O)(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfonamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OS(=O)(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfonic acid"),
    ],
)
def test_von_baeyer_spiro_sulfonic_acid_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_sulfonic_acid_ring_unsaturation():
    assert smiles_to_iupac("OS(=O)(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfonic acid"


def test_von_baeyer_tellone_ring_unsaturation():
    assert smiles_to_iupac("[Te]=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-tellone"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[TeH]C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-tellurol"),
    ],
)
def test_von_baeyer_spiro_tellurol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_tellurol_ring_unsaturation():
    assert smiles_to_iupac("[TeH]C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-tellurol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("S=C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-thione"),
    ],
)
def test_von_baeyer_spiro_thione_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_thione_ring_unsaturation():
    assert smiles_to_iupac("S=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-thione"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("SC12CC3CC(CC(C3)C1)C2", "adamantane-1-thiol"),
    ],
)
def test_von_baeyer_thiol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_thiol_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("SCC1CC2CCC1C2") == "(bicyclo[2.2.1]heptan-2-yl)methanethiol"


def _fullerene_derivative(cage_smiles, sites, numbered=True):
    """The cage with the atoms at systematic locants `sites` (atom index + 1 when not `numbered`) saturated: a
    substituent SMILES, None for H, or a charge."""
    import networkx as nx
    from smiles_to_iupac._fullerene_numbering import fullerene_numberings

    cage = Chem.MolFromSmiles(cage_smiles)
    graph = nx.Graph((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in cage.GetBonds())
    atom_of = (
        {locant: atom for atom, locant in fullerene_numberings(graph)[0].items()}
        if numbered
        else {atom + 1: atom for atom in graph}
    )
    rest = graph.subgraph(set(graph) - {atom_of[locant] for locant in sites})
    double = {frozenset(pair) for pair in nx.max_weight_matching(rest, maxcardinality=True)}
    mol = Chem.RWMol()
    for _ in range(cage.GetNumAtoms()):
        mol.AddAtom(Chem.Atom(6))
    for u, v in graph.edges:
        mol.AddBond(u, v, Chem.BondType.DOUBLE if frozenset((u, v)) in double else Chem.BondType.SINGLE)
    for locant, tail in sites.items():
        if tail == "-":
            mol.GetAtomWithIdx(atom_of[locant]).SetFormalCharge(-1)
        elif tail:
            offset = mol.GetNumAtoms()
            mol = Chem.RWMol(Chem.CombineMols(mol, Chem.MolFromSmiles(tail)))
            mol.AddBond(atom_of[locant], offset, Chem.BondType.SINGLE)
    derivative = mol.GetMol()
    Chem.SanitizeMol(derivative)
    return Chem.MolToSmiles(derivative)


_ACETIC = "CC(=O)O"


@pytest.mark.parametrize(
    "cage,sites,expected",
    [
        # P-29.3.4.1 / P-71.2: a fullerene yl group takes the added hydrogen (1(9H)), further pairs take hydro prefixes
        (_FULLERENE_C60_SMILES, {1: _ACETIC, 9: None}, "[(C60-Ih)[5,6]fulleren-1(9H)-yl]acetic acid"),
        (_FULLERENE_C60_SMILES, {1: _ACETIC, 2: None}, "[(C60-Ih)[5,6]fulleren-1(2H)-yl]acetic acid"),
        (_FULLERENE_C60_SMILES, {1: _ACETIC, 9: "C"}, "[9-methyl(C60-Ih)[5,6]fulleren-1(9H)-yl]acetic acid"),
        (_FULLERENE_C60_SMILES, {1: _ACETIC, 9: None, 52: None, 60: None}, "[52,60-dihydro(C60-Ih)[5,6]fulleren-1(9H)-yl]acetic acid"),
        (_FULLERENE_C70_SMILES, {7: _ACETIC, 22: None}, "{(C70-D5h(6))[5,6]fulleren-7(22H)-yl}acetic acid"),
        (_FULLERENE_C60_SMILES, {1: _ACETIC, 9: _ACETIC}, "2,2'-[(C60-Ih)[5,6]fullerene-1,9-diyl]diacetic acid"),
        # P-31.1.4 / P-6: the cage as parent, substituent prefixes before the hydro prefixes (the Blue Book examples)
        (_FULLERENE_C60_SMILES, {1: "C(F)(F)F", 9: None}, "1-(trifluoromethyl)-1,9-dihydro(C60-Ih)[5,6]fullerene"),
        (
            _FULLERENE_C60_SMILES,
            {1: "F", 9: "F", 52: "F", 60: "F"},
            "1,9,52,60-tetrafluoro-1,9,52,60-tetrahydro(C60-Ih)[5,6]fullerene",
        ),
        pytest.param(_FULLERENE_C60_SMILES, {1: "C(C)(C)C", 7: "c1ccccc1"}, "1-tert-butyl-7-phenyl-1,7-dihydro(C60-Ih)[5,6]fullerene", marks=pytest.mark.slow),
        (_FULLERENE_C60_SMILES, {1: "C", 23: "C"}, "1,23-dimethyl-1,23-dihydro(C60-Ih)[5,6]fullerene"),
        # P-72: the carbanion with its added hydrogen
        (_FULLERENE_C60_SMILES, {1: "-", 9: None}, "(C60-Ih)[5,6]fulleren-1(9H)-ide"),
    ],
)
def test_fullerene_derivative_names(cage, sites, expected):
    assert smiles_to_iupac(_fullerene_derivative(cage, sites)) == expected


@pytest.mark.parametrize(
    "cage,pairs",
    [
        # 1-9 in C60 and 8-25 in C70 (the C70 PCBM fusion locants) are [6,6] bonds; the C70 pathway ends on the
        # [6,6] bond bisected by its C2 axis (Fu-3.3)
        (_FULLERENE_C60_SMILES, {(1, 9): (6, 6), (1, 2): (5, 6)}),
        (_FULLERENE_C70_SMILES, {(8, 25): (6, 6), (7, 8): (5, 6), (69, 70): (6, 6)}),
    ],
)
def test_fullerene_numbering_follows_the_spiral_pathway(cage, pairs):
    import networkx as nx
    from smiles_to_iupac._fullerene_numbering import fullerene_numberings

    mol = Chem.MolFromSmiles(cage)
    graph = nx.Graph((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds())
    rings = mol.GetRingInfo().AtomRings()
    for numbering in fullerene_numberings(graph):
        atom = {locant: idx for idx, locant in numbering.items()}
        for (first, second), fusion in pairs.items():
            assert graph.has_edge(atom[first], atom[second])
            assert tuple(sorted(len(r) for r in rings if atom[first] in r and atom[second] in r)) == fusion


def test_fullerene_derivative_outside_the_numbered_cages_stops():
    from smiles_to_iupac._fullerene import _FULLERENE_C76_SMILES

    # atoms 0 and 1 of the C76 reference SMILES are bonded
    with pytest.raises(UnsupportedStructure, match="Fu-3.1"):
        smiles_to_iupac(_fullerene_derivative(_FULLERENE_C76_SMILES, {1: _ACETIC, 2: None}, numbered=False))


def test_caged_polycycle_beyond_a_von_baeyer_or_bridged_fused_search_is_refused():
    homo_dihydro_c60 = (
        "C1C2c3c4c5c6c7c8c(c9c%10c2c2c3c3c%11c4c4c5c5c7c7c%12c8c8c9c9c%10c%10c2c2c3c3c%11c%11c4c4c5c7c5c7c%12c8c8c9c9c%10"
        "c2c2c3c3c%11c4c5c4c7c8c9c2c34)C16"
    )
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(homo_dihydro_c60)


def test_fullerene_with_characteristic_group_is_named_as_the_parent_with_added_hydrogen():
    # P-6 (P-58.2): the principal characteristic group takes the added hydrogen, not hydro prefixes
    assert smiles_to_iupac(_fullerene_derivative(_FULLERENE_C60_SMILES, {1: "O", 9: None})) == (
        "(C60-Ih)[5,6]fulleren-1(9H)-ol"
    )


def test_bridged_fused_naming_does_not_write_rdkit_logs_to_stderr(capfd):
    smiles_to_iupac(_fullerene_derivative(_FULLERENE_C60_SMILES, {1: _ACETIC, 9: None}))
    assert capfd.readouterr().err == ""


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C1CCC2(CC1)OCCO2", "1,4-dioxaspiro[4.5]decan-8-one", id="heteroatom_spiro_ketone"),
        pytest.param("O=C1OC2CCC1C2", "2-oxabicyclo[2.2.1]heptan-3-one", id="bridged_lactone"),
        pytest.param("OC1CCCC12C=CCCC2", "spiro[4.5]dec-6-en-1-ol", id="unsaturated_spiro_alcohol"),
        pytest.param("CNS(=O)(=O)C1CC2CCC1C2", "N-methylbicyclo[2.2.1]heptane-2-sulfonamide", id="n_substituted_bridged_sulfonamide"),
        pytest.param("OC1CC2CCC1N2C", "7-methyl-7-azabicyclo[2.2.1]heptan-2-ol", id="ring_substituent_on_hetero_bridge_atom"),
    ],
)
def test_principal_groups_on_hetero_bridged_and_spiro_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)(C)C1CCCC12CC2", "4-tert-butylspiro[2.4]heptane", id="spiro"),
        pytest.param("CC(C)(C)C1CC2CCC1N2", "2-tert-butyl-7-azabicyclo[2.2.1]heptane", id="hetero_bridged"),
    ],
)
def test_tert_butyl_substituent_on_spiro_and_hetero_bridged_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C=C1CC2CCC1C2", "2-methylidenebicyclo[2.2.1]heptane"),
        ("CC=C1CC2CC1C=C2", "5-ethylidenebicyclo[2.2.1]hept-2-ene"),
        ("CC=C1CC2CC1C1=C2C(C)C1", "7-ethylidene-3-methyltricyclo[4.2.1.0^2,5]non-2(5)-ene"),
    ],
)
def test_von_baeyer_exocyclic_double_bond_is_an_ylidene_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CCC2(CC1)C=Cc1ccccc12", "spiro[cyclohexane-1,1′-indene]", id="carbocycle_and_fused_component_without_indicated_hydrogen"),
        pytest.param("C1CNCCC12c3ccccc3Oc4ccccc24", "spiro[piperidine-4,9′-xanthene]", id="heteromonocycle_component"),
        pytest.param("C1CC2(CCNCC2)c2ccccc12", "2,3-dihydrospiro[indene-1,4′-piperidine]", id="hydro_prefixes_of_the_union"),
        pytest.param("Cc1ccc2c(c1)C=CC21CCCCC1", "5′-methylspiro[cyclohexane-1,1′-indene]", id="primed_substituent_locant"),
        pytest.param("C1CC=CC2(C1)Cc1ccccc1C2", "1′,3′-dihydrospiro[cyclohexane-1,2′-inden]-2-ene", id="double_bond_cited_as_an_ending_of_the_union"),
        pytest.param("C1NCC2(N1)C=Nc1ccccc1N2", "1′H-spiro[imidazolidine-4,2′-quinoxaline]", id="indicated_hydrogen_of_the_second_component"),
        pytest.param("C12(CC3=CC=CC=C3C=C1)CC4=CC=CC=C4C=C2", "1H,1′H-2,2′-spirobi[naphthalene]", id="identical_components_with_indicated_hydrogen"),
        pytest.param("C12(Oc3ccccc3C=C1)C=COc4ccccc42", "2,4′-spirobi[[1]benzopyran]", id="locants_of_a_benzo_component_in_brackets"),
        pytest.param("C12(CC=C3C=CC=CC=C13)CC4=CC=CC=CC4=C2", "1′H,2H-1,2′-spirobi[azulene]", id="lower_spiro_locant_is_unprimed"),
        pytest.param("C12(SSc3ccccc31)Sc4ccccc4S2", "spiro[[1,2]benzodithiole-3,2′-[1,3]benzodithiole]", id="components_ordered_by_heteroatom_locants"),
        pytest.param("[S]12(Oc3ccccc3O1)Oc3ccccc3O2", "2λ4,2′-spirobi[[1,3,2]benzodioxathiole]", id="lambda_spiro_atom_without_double_bond"),
        pytest.param("P12(=Nc3ccccc3C=N1)NC4=CC=CC=C4C=N2", "1H-2λ5,2′-spirobi[[1,3,2]benzodiazaphosphinine]", id="lambda_spiro_atom_taking_part_in_the_mancude_system"),
        pytest.param("CP12(Oc3ccccc3O1)Oc1ccccc1O2", "2-methyl-2H-2λ5,2′-spirobi[[1,3,2]benzodioxaphosphole]", id="lambda_spiro_atom_substituent_stands_for_the_indicated_hydrogen"),
        pytest.param("N1=P2(N=C[PH2]=C1)OC1=C(N2)C=CC=C1", "3H-2λ5,5′λ5-spiro[[1,3,2]benzoxazaphosphole-2,2′-[1,3,2,5]diazadiphosphinine]", id="second_lambda_atom_in_a_mancude_component"),
        pytest.param("N[C@H]1C[C@]12CCc1ccccc12", "(1S,2S)-2′,3′-dihydrospiro[cyclopropane-1,1′-inden]-2-amine", id="stereocentres_with_a_suffix_group"),
        pytest.param("N[C@H]1C[C@]12CC[C@@H](C)c1ccccc12", "(1S,2S,4′R)-4′-methyl-3′,4′-dihydro-2′H-spiro[cyclopropane-1,1′-naphthalen]-2-amine", id="stereocentre_at_a_primed_locant"),
        pytest.param("C1[N+]2(C=CC3=CC=CC=C13)C=CN4C=CC=CC4=C2", "1H-2λ5-spiro[isoquinoline-2,2′-pyrido[1,2-a]pyrazin]-2-ylium", id="cationic_spiro_atom_of_different_components"),
        pytest.param("C1=CC=C2C(=C1)C[N+]3(C2)CC4=CC=CC=C4C3", "1,1′,3,3′-tetrahydro-2λ5,2′-spirobi[isoindol]-2-ylium", id="cationic_spiro_atom_of_identical_components"),
    ],
)
def test_monospiro_union_with_a_polycyclic_component(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C12CCC(CC1)SC23c1ccccc1-c1ccccc13", "3-thiaspiro[bicyclo[2.2.2]octane-2,9′-fluorene]", id="replacement_prefix_before_spiro_with_low_spiro_locant_first"),
        pytest.param("C1CC2CC1CC23OCCCCCCCCCO3", "2′,12′-dioxaspiro[bicyclo[2.2.1]heptane-2,1′-cyclododecane]", id="monocycle_beyond_hantzsch_widman_named_by_replacement"),
        pytest.param("C1CC2CCC1CC21CC2CCC(C2)C1", "spiro[bicyclo[2.2.2]octane-2,3′-bicyclo[3.2.1]octane]", id="descriptor_numbers_order_components"),
        pytest.param("C1OCC2CC1CC21OC2CCC1CC2", "3,3′-dioxaspiro[bicyclo[2.2.2]octane-2,6′-bicyclo[3.2.1]octane]", id="spiro_locants_before_heteroatom_locants"),
        pytest.param("C1CC2SCC1CC21CC2CCC1CS2", "5,6′-dithia-2,2′-spirobi[bicyclo[2.2.2]octane]", id="identical_von_baeyer_components_with_heteroatoms"),
        pytest.param("C1CC2CC1CC[Si]21CC2CCC(C2)C1", "2-sila-2,3′-spirobi[bicyclo[3.2.1]octane]", id="standard_valence_heteroatom_at_the_spiro_atom"),
        pytest.param("C1CC[N+]2(C1)CCC1CCC(C1)C2", "1′λ5-spiro[3-azabicyclo[4.2.1]nonane-3,1′-pyrrolidin]-1′-ylium", id="lowest_spiro_locant_cited_with_lambda_and_ylium"),
        pytest.param("C1CCC2(C1)SC1CC2C2CC12", "7′-thiaspiro[cyclopentane-1,6′-tricyclo[3.2.1.0^2,4]octane]", id="polycyclic_von_baeyer_component"),
        pytest.param("C1CCC2(C1)C1CC3CC(C1)CC2C3", "spiro[adamantane-2,1′-cyclopentane]", id="adamantane_retained_component"),
        pytest.param("C12C3(C4OCC(C3)C=C4)CC(OC1)C=C2", "5,6′-dioxa-2,2′-spirobi[bicyclo[2.2.2]octane]-7,7′-diene", id="ene_endings_follow_the_last_bracket"),
        pytest.param("C1=CC2CC1CC23CC4CCC3C4", "2,2′-spirobi[bicyclo[2.2.1]heptan]-5-ene", id="single_ene_elides_the_final_e"),
    ],
)
def test_monospiro_union_with_a_von_baeyer_component(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "C[C@H]1C[C@H]2CC[C@@H]1C21CCCC1",
            "(1S,2S,4R)-2-methylspiro[bicyclo[2.2.1]heptane-7,1\u2032-cyclopentane]",
            id="bridgeheads_and_a_substituted_atom_of_the_von_baeyer_component",
        ),
        pytest.param(
            "C[C@@H]1C[C@@H]2CC[C@H]1C21CCCC1",
            "(1R,2R,4S)-2-methylspiro[bicyclo[2.2.1]heptane-7,1\u2032-cyclopentane]",
            id="enantiomer_inverts_every_descriptor",
        ),
        pytest.param(
            "CC1(C)[C@H]2CCC3(CCCC3)[C@H]21",
            "(1S,5S)-6,6-dimethylspiro[bicyclo[3.1.0]hexane-2,1\u2032-cyclopentane]",
            id="bicyclo_hexane_bridgeheads",
        ),
        pytest.param(
            "C[C@@H]1CC[C@H]2CCC3(CCCC3)[C@@H]12",
            "(3\u2032aS,6\u2032R,6\u2032aS)-6\u2032-methyl-3\u2032,3\u2032a,4\u2032,5\u2032,6\u2032,6\u2032a-hexahydro-2\u2032H-spiro[cyclopentane-1,1\u2032-pentalene]",
            id="primed_locants_of_a_fused_component",
        ),
        pytest.param(
            "Cl[C@@H]1C2CCCOC2C12CCCC2",
            "(8R)-8-chloro-5-oxaspiro[bicyclo[4.2.0]octane-7,1\u2032-cyclopentane]",
            id="only_the_specified_centre_is_cited_beside_a_replacement_prefix",
        ),
        pytest.param(
            "C1C[C@@H]2C[C@H]1C[C@]23CCCC[C@H]3C",
            "(1R,1\u2032S,2\u2032R,4S)-2\u2032-methylspiro[bicyclo[2.2.1]heptane-2,1\u2032-cyclohexane]",
            id="unprimed_locants_before_primed_ones_in_the_citation",
        ),
    ],
)
def test_monospiro_union_with_a_polycyclic_component_cites_the_configuration(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CC2(C=CNC=C2)C2=CC=CC=C12", "2,3-dihydro-1′H-spiro[indene-1,4′-pyridine]", id="mancude_retained_monocycle_with_indicated_hydrogen"),
        pytest.param("C1CC2(CCN(C=C2)C)c2ccccc12", "1′-methyl-2,2′,3,3′-tetrahydro-1′H-spiro[indene-1,4′-pyridine]", id="hydro_prefixes_of_a_partly_saturated_monocycle"),
        pytest.param("C12C[N+]3(CC(CC1)CC2)COC=C3", "2′H-3λ5-spiro[3-azabicyclo[3.2.2]nonane-3,3′-[1,3]oxazol]-3-ylium", id="hantzsch_widman_mancude_monocycle_with_cationic_spiro_atom"),
        pytest.param("P124(OC3CC(O1)CC(O2)C3)OC=NC4", "3H-1′λ5-spiro[[1,4,2]oxazaphosphole-2,1′-[2,8,9]trioxa[1]phosphaadamantane]", id="replaced_adamantane_component_brackets_its_prefix_locants"),
    ],
)
def test_monospiro_union_with_an_unsaturated_hetero_monocycle(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1=CC2(C=Cc3ccccc3C2)C=C2CC3(C=Cc4ccccc4C3)CC=C12", "1H,1′H,1′′H,3′H-2,2′:7′,2′′-dispiroter[naphthalene]", id="three_identical_fused_components_with_indicated_hydrogen"),
        pytest.param("C1CCC2C(C1)C21C2CCC3(CCCC4OC43)CC21", "7-oxa-2,3′:7′,7′′-dispiroter[bicyclo[4.1.0]heptane]", id="three_identical_von_baeyer_components_with_replacement_prefix"),
        pytest.param("C12C3(CCC(OC1)O2)OC2CCC1(C3O2)OC2CCCC1O2", "6,6′,6′′,8,8′,8′′-hexaoxa-2,7′:2′,7′′-dispiroter[bicyclo[3.2.1]octane]", id="replacement_heteroatoms_take_low_locants_before_the_spiro_atoms"),
        pytest.param("C12C3(CCCC2O1)CC1C2(C1CC3)C3C=CCOC32", "2′′,7-dioxa-2,3′:7′,7′′-dispiroter[bicyclo[4.1.0]heptan]-4′′-ene", id="equal_hetero_numbers_leave_the_spiro_citation_order_to_decide"),
        pytest.param("O1C2(C=CC=C1)OCC1(OC2)C2OCC(O1)C2", "3,6-dioxadispiro[bicyclo[2.2.1]heptane-2,2′-[1,4]dioxane-5′,2′′-pyran]", id="different_components_take_low_spiro_locants_before_the_heteroatoms"),
        pytest.param("c1ccc2c(c1)OS13(O2)(Oc2ccccc2O1)Oc1ccccc1O3", "2λ6,2′,2′′-spiroter[[1,3,2]benzodioxathiole]", id="three_identical_components_on_one_nonstandard_atom"),
        pytest.param("c1ccc2c(c1)OS13(O2)(Oc2ccccc2O1)Oc1ccccc1S3", "2λ6-spiro[bis([1,3,2]benzodioxathiole)-2,2′′:2′,2′′-[1,2,3]benzoxadithiole]", id="two_identical_components_and_a_third_on_one_atom"),
        pytest.param("c1ccc2c(c1)OS13(O2)(Oc2ccccc2S1)c1ccccc1-c1ccccc13", "2λ6-spiro[[1,3,2]benzodioxathiole-2,2′-([1,2,3]benzoxadithiole)-2,5′′-dibenzo[b,d]thiophene]", id="three_different_components_on_one_atom"),
        pytest.param("C1=CS23(C=CC14CCCC4)(Oc1ccccc1O2)Oc1ccccc1O3", "1′′λ6-dispiro[bis([1,3,2]benzodioxathiole)-2,1′′:2′,1′′-thiopyran-4′′,1′′′-cyclopentane]", id="central_component_with_terminals_sharing_a_nonstandard_atom"),
        pytest.param("C1C2(C3C4C3C3C4C32)C2(C3C4C3C3C4C32)C12OC21C2C3C2C2C3C21", "trispiro{1-oxaspiro[2.3]hexane-2,3′:4,3′′:5,3′′′-tris(tetracyclo[3.2.0.0^2,7.0^4,6]heptane)}", id="monocyclic_unit_as_central_component_of_identical_terminals"),
        pytest.param("C1CCC2(CC1)CC1(OC3=C(S1)C1(CCCCC1)OC31CCCCC1)C1(CCCCC1)O2", "trispiro{bis(cyclohexane)-1,4′:6′,1′′-furo[3,4-d][1,3]oxathiole-2′,14′′′-[7]oxadispiro[5.1.5^8.2^6]pentadecane}", id="monocyclic_unit_as_terminal_component_beside_repeated_terminals"),
    ],
)
def test_polyspiro_union_with_polycyclic_components(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1C2(Cc3ccccc3C24CCCC4)CC3(Cc4ccccc4C3)CC15Cc6ccccc6C5", "1′′′,1′′′ ′,3′′′,3′′′ ′-tetrahydro-3′′H-spiro{cyclopentane-1,1′′-trispiro[cyclohexane-1,2′:3,2′′:5,2′′′-tris(indene)]}", id="largest_branched_system_named_first_and_used_as_a_unit"),
        pytest.param("c1ccc2c(c1)CC1(C2)OP23(S1)SP14(OC5(Cc6ccccc6C5)S1)SP1(OC5(Cc6ccccc6C5)S1)(S2)SP1(OC2(Cc5ccccc5C2)S1)(S3)S4", "1′′′ ′′,1′′′ ′′′,1′′′ ′′′ ′,1′′′ ′′′ ′′,3′′′ ′′,3′′′ ′′′,3′′′ ′′′ ′,3′′′ ′′′ ′′-octahydro-1λ5,2′′λ5,2′′′λ5,2′′′ ′λ5-tetraspiro{tetraspiro[[2,4,6,8,9,10]hexathia[1,3,5,7]tetraphosphaadamantane-1,2′:3,2′′:5,2′′′:7,2′′′ ′-tetrakis([1,3,2]oxathiaphosphetane)]-4′,2′′′ ′′:4′′,2′′′ ′′′:4′′′,2′′′ ′′′ ′:4′′′ ′,2′′′ ′′′ ′′-tetrakis(indene)}", id="heteroadamantane_unit_with_nonstandard_spiro_atoms"),
        pytest.param("C1C2(Cc3ccccc3C24CCCC4)CC3(Cc4ccccc4C3)CC15Cc6ccccc6C5=O", "1′′′ ′,3′′′ ′-dihydro-3′′H-spiro{cyclopentane-1,1′′-trispiro[cyclohexane-1,2′:3,2′′:5,2′′′-tris(indene)]}-1′′′(3′′′H)-one", id="suffix_group_on_a_nested_system"),
    ],
)
def test_polyspiro_union_with_a_nested_spiro_system(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nested_spiro_systems_of_equal_size_have_no_preferred_name():
    with pytest.warns(NonPreferredNameWarning, match="P-24.7.4"):
        name = smiles_to_iupac("C%10%11(C2(CCCC2)c3ccccc3C4(CCC4)%10)C5(CCCCC5)c6ccccc6C7(CCCC7)%11")
    assert name.startswith("dispiro{cyclohexane-1,")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C1CC2(CC3CCC2C3)c2ccccc12", "spiro[bicyclo[2.2.1]heptane-2,1′-inden]-3′(2′H)-one", id="ketone_with_added_hydrogen_on_a_union_with_a_von_baeyer_component"),
        pytest.param("C1=CC2(CCC1)CC1=C(C=CCC1)O2", "4,5-dihydro-3H-spiro[[1]benzofuran-2,1′-cyclohexan]-2′-ene", id="double_bond_of_the_second_component_cited_as_an_ending"),
        pytest.param("O=C1C=CC2(CC1)CC=Cc1ccccc12", "2′H-spiro[cyclohexane-1,1′-naphthalen]-2-en-4-one", id="double_bond_ending_before_a_suffix"),
        pytest.param("O=C(O)C1CC2(CCCCC2)c2ccccc12", "2′,3′-dihydrospiro[cyclohexane-1,1′-indene]-3′-carboxylic acid", id="suffix_without_elision_after_hydro_prefixes"),
        pytest.param("CC1CC2(CCCCC2)c2ccccc12", "3′-methyl-2′,3′-dihydrospiro[cyclohexane-1,1′-indene]", id="primed_substituent_locant_with_hydro_prefixes"),
    ],
)
def test_spiro_union_with_suffix_groups_and_double_bonds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCCC12[PH3]CCCC2", "6λ5-phosphaspiro[4.5]decane"),
        ("C1CC2CC[PH3]C(C1)C2", "2λ5-phosphabicyclo[3.3.1]nonane"),
        ("C1CC2CCC1[SH2]2", "7λ4-thiabicyclo[2.2.1]heptane"),
        ("C1CCC[PH3]1", "1λ5-phospholane"),
        ("C1CCCCCCCCCCC[SH2]1", "1λ4-thiacyclotridecane"),
        ("S1CCCCCCCCCCC1", "thiacyclododecane"),
        ("N1=CO1", "oxazirene"),
        ("OC(=O)C1CCCCCCCCCCCO1", "1-oxacyclotridecane-2-carboxylic acid"),
        ("C1CC2CC[SH2]C2C1", "hexahydro-2H-1λ4-cyclopenta[b]thiophene"),
        ("C1CC2(C1)CCSCC2", "7-thiaspiro[3.5]nonane"),
    ],
)
def test_skeletal_heteroatom_with_a_nonstandard_bonding_number_cites_lambda(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C12C3CCCCCCC(C(CCCCCC1)CCCCCC2)C2CCC(C3)C2", "tetracyclo[8.6.6.5^2,9.1^23,26]octacosane"),
        ("C12CC3CCCCC4CCCCC(CC(CCC5CC5CC1)CCCC2)CC(C4)C3", "pentacyclo[13.7.4.3^3,8.0^18,20.1^13,28]triacontane"),
    ],
)
def test_largest_main_bridge_with_dependent_secondary_bridges(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "S1C2(CCCCC2)CC3(C=C4C=CC=CC4=C3)SC5(CCCCC5)CC1",
            "trispiro[bis(cyclohexane)-1,2′:6′,1′′-[1,5]dithiocane-4′,2′′′-indene]",
            id="second_copy_of_the_first_cited_component_follows_the_central_one",
        ),
        pytest.param(
            "C=1C2(C=C3C=CC=CC13)CC1(SCCC3(CCCC3)S2)CCCCC1",
            "trispiro[cyclohexane-1,2′-[1,5]dithiocane-6′,1′′-cyclopentane-4′,2′′′-indene]",
            id="different_terminal_components",
        ),
    ],
)
def test_branched_spiro_union_pair_order(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CC2(C1)CC[PH3]CC2", "7λ5-phosphaspiro[3.5]nonane", id="lambda5_phosphorus"),
        pytest.param("C1[SiH2]CC12CC[PH3]CC2", "7λ5-phospha-2-silaspiro[3.5]nonane", id="lambda5_phosphorus_with_silicon"),
        pytest.param("C1CCS12CCCCC2", "4λ4-thiaspiro[3.5]nonane", id="lambda4_spiro_sulfur"),
        pytest.param("C1[SH4]C[SH2]CC12CCCCC2", "2λ6,4λ4-dithiaspiro[5.5]undecane", id="higher_bonding_number_gets_the_lower_locant"),
        pytest.param(
            "S%11%12%13(OCCO%11)(OCCO%12)OCCO%13",
            "1,4,6,9,10,13-hexaoxa-5λ6-thiaspiro[4.4^5.4^5]tridecane",
            id="three_rings_on_a_lambda6_spiro_atom",
        ),
        pytest.param(
            "S%11%12%13(CC%11)(CCCC%12)CCCCC%13", "3λ6-thiaspiro[2.4^3.5^3]dodecane", id="smaller_ring_numbered_first"
        ),
    ],
)
def test_spiro_systems_with_atoms_of_nonstandard_bonding_number(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[SiH]12O[SiH2]O[SiH](O[SiH2]O1)O2", "bicyclo[3.3.1]tetrasiloxane", id="bicyclic_siloxane"),
        pytest.param("[SiH]12O[SiH]3O[SiH](O1)O[SiH](O2)O3", "tricyclo[3.3.1.1^3,7]tetrasiloxane", id="adamantane_shaped_siloxane_is_the_alternating_cage"),
        pytest.param("O1[SiH]2O[SiH]3O[SiH]1O[SiH]1O[SiH](O[SiH](O1)O3)O2", "tetracyclo[5.5.1.1^3,11.1^5,9]hexasiloxane", id="hexasilasesquioxane"),
        pytest.param("[Si]12(O[SiH2]O[SiH2]O1)O[SiH2]O[SiH2]O2", "spiro[5.5]pentasiloxane", id="spiro_siloxane"),
        pytest.param("[SiH]12[SiH2][SiH2][SiH2][SiH2][SiH]([SiH2][SiH2]1)[SiH2]2", "bicyclo[4.2.1]nonasilane", id="same_element_von_baeyer"),
        pytest.param("[SiH2]1[SiH2][SiH2][Si]2([SiH2]1)[SiH2][SiH2][SiH2][SiH2]2", "spiro[4.4]nonasilane", id="same_element_spiro"),
        pytest.param("P12PPP(P1)PP2", "bicyclo[2.2.1]heptaphosphane", id="same_element_cage_perceived_aromatic"),
        pytest.param("N12[SiH]3[NH][SiH]([NH]3)N([SiH2]2)[SiH2][NH][SiH2]1", "1N-tricyclo[3.3.1.1^2,4]pentasilazane", id="nitrogen_at_locant_one"),
        pytest.param("N1[SiH2]N[SiH]2N[SiH]1N1[SiH2]N2[SiH2]1", "1Si-tricyclo[3.3.1.1^2,4]pentasilazane", id="silicon_at_locant_one"),
    ],
)
def test_alternating_heteroatom_cages_take_the_preselected_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O1CCCC12CC1(CCCN1)CC2", "1-oxa-8-azadispiro[4.1.4^7.2^5]tridecane"),
        ("O1CCC[C@@]12C[C@]1(OCCC1)CC2", "(5R,7R)-1,8-dioxadispiro[4.1.4^7.2^5]tridecane"),
        ("O1[C@]2(N=CC1)C=CC1=CC=CC=C12", "(1R)-5′H-spiro[indene-1,2′-[1,3]oxazole]"),
    ],
)
def test_heteroatom_polyspiro_systems_and_stereogenic_spiro_atoms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C123CCCC(C(CCCC1)C2)C3", "tricyclo[4.4.1.1^1,5]dodecane"),
        ("C12C3C4CCC(CC(CCC1)CCC2)(C4)CC3", "tetracyclo[6.3.3.2^2,6.1^3,6]heptadecane"),
        ("C12CC34CCCC(CC(CCCC1)CCC2)(CC3)C4", "tetracyclo[7.4.3.2^3,7.1^3,7]nonadecane"),
        ("C12C3C4C1C15C4C14C3C2C54", "heptacyclo[4.4.0.0^1,5.0^2,5.0^3,9.0^4,8.0^7,10]decane"),
        ("C1CC12CCC1(CCC3(CC3)CCC3(CC3)CC1)CC2", "tetraspiro[2.2.2.2^9.2.2^14.2^6.2^3]icosane"),
        ("C1CC12C1(CC1)CC1(CC3(CC3)C3(CC3)C1)CC2", "pentaspiro[2.0.2^4.1.1.2^10.0.2^13.1^8.2^3]octadecane"),
        (
            "C1CC12C1(OOC3(C4(CC4)CCC34CC4)OO1)C1(CC1)CC2",
            "5,6,16,17-tetraoxahexaspiro[2.0.2.0.2^8.2.2^13.0^7.2^4.0.2^18.2^3]docosane",
        ),
        ("C12CC#C\\C=C/C=C/CCC(=CCC1)C2", "(5Z,7E)-bicyclo[9.3.1]pentadeca-5,7,11-trien-3-yne"),
        ("C1=2CCCCCCCCCCCCCCCCC(=CC=C1)C2", "bicyclo[16.3.1]docosa-1(22),18,20-triene"),
    ],
)
def test_general_von_baeyer_and_branched_polyspiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def _senior_to(senior, junior, assembly=False):
    a, b = Chem.MolFromSmiles(senior), Chem.MolFromSmiles(junior)
    return ring_seniority_key(a, range(a.GetNumAtoms()), assembly) < ring_seniority_key(b, range(b.GetNumAtoms()), assembly)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CCC2(C1)CCN(CC2)Cc1ccc2ccccc2n1", "8-[(quinolin-2-yl)methyl]-8-azaspiro[4.5]decane", id="spiro_before_fused"),
        pytest.param(
            "c1ccc2ccccc2c1CC1CCC2CCC1CC2",
            "1-[(bicyclo[3.2.2]nonan-2-yl)methyl]naphthalene",
            id="fused_before_bridged_nonfused_with_rings_counted_by_cyclomatic_number",
        ),
        pytest.param(
            "C1CC2CCC1CC2CC1CC2CCC(C1)C2",
            "2-[(bicyclo[3.2.1]octan-3-yl)methyl]bicyclo[2.2.2]octane",
            id="lower_ring_size_descriptor",
        ),
    ],
)
def test_parent_of_two_polycyclic_systems_follows_the_order_of_p44_2_2_2(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "senior,junior",
    [
        pytest.param("C1CCC2(C1)NCC3(CCCC3)CC2", "C1CCC2(C1)NCCc1ccccc12", id="more_spiro_fusions"),
        pytest.param(
            "C1CCC2(C1)CCCC3(C2)CCCC3", "C1CCC2(C1)CCC3(CC2)CCCC3", id="lower_spiro_locants_of_saturated_monocycles"
        ),
        pytest.param(
            "C1CCC2(C1)CCCC3(C2)CCCC3", "C1=CCC2(C1)CCCC3(C2)CCCC3", id="saturated_monocycles_before_unsaturated"
        ),
        pytest.param("C1CCC2(C1)NCCc1ccccc12", "C1CCC2(C1)CCNCC2", id="senior_component"),
        pytest.param(
            "C1COC2(C1)CN1CCC2CC12OCCO2",
            "C1COC2(CCOC23CN2CCC3CC2)O1",
            id="senior_component_cited_first_with_replacement_heteroatoms_as_carbon",
        ),
        pytest.param("C1CCC2(C1)NCCc1ccccc12", "C1CCC2(C1)CNCc1ccccc12", id="lower_spiro_locants_of_a_component"),
        pytest.param("C1CC2CCC1CC2", "C1CC2CCC(C1)C2", id="von_baeyer_lower_ring_size_descriptor"),
        pytest.param("C1CC2CCC3(CCC1C3)C2", "C1CC2CC1CC1CCC2C1", id="von_baeyer_lower_superscripts_as_a_set"),
        pytest.param(
            "C1CC2CC(C1)C1CCC2CC2CCC1CC2",
            "C1CC2CCC(C1)C1CCC2C2CCCC1C2",
            id="von_baeyer_lower_superscripts_in_order_of_citation",
        ),
    ],
)
def test_senior_spiro_and_von_baeyer_system(senior, junior):
    assert _senior_to(senior, junior)


def _modified_fullerene(cage_smiles, deleted=(), cleaved=(), formed=(), inserted=(), replaced=None, saturated=(), inserted_element=6, attached=()):
    """The cage at systematic locants with atoms `deleted`, bonds `cleaved`, bonds `formed`, atoms `inserted` into bonds
    and atoms `replaced` ({locant: atomic number}); the atoms at `saturated` locants carry the extra hydrogen. `attached`
    are (cage locants, SMILES) of rings and bridges, written with [Xe:1], [Xe:2], ... for those atoms (a locant 'insertedN'
    is the Nth atom inserted); a bond between two of them is the cage bond itself."""
    import networkx as nx
    from smiles_to_iupac._fullerene_numbering import fullerene_numberings

    cage = Chem.MolFromSmiles(cage_smiles)
    graph = nx.Graph((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in cage.GetBonds())
    locants = fullerene_numberings(graph)[0]
    bonds = {frozenset((locants[u], locants[v])) for u, v in graph.edges}
    bonds -= {frozenset(pair) for pair in cleaved} | {frozenset(pair) for pair in inserted}
    bonds |= {frozenset(pair) for pair in formed}
    nodes = [locant for locant in sorted(locants.values()) if locant not in deleted]
    bonds = {bond for bond in bonds if not bond & set(deleted)}
    inserted_names = [f"inserted{i}" for i in range(len(inserted))]
    for name, pair in zip(inserted_names, inserted):
        bonds |= {frozenset((pair[0], name)), frozenset((pair[1], name))}
    nodes += inserted_names
    replaced = replaced or {}
    parts = [(places, Chem.MolFromSmiles(smiles)) for places, smiles in attached]
    cage_bonds = {
        frozenset((places[b.GetBeginAtom().GetAtomMapNum() - 1], places[b.GetEndAtom().GetAtomMapNum() - 1]))
        for places, part in parts
        for b in part.GetBonds()
        if b.GetBeginAtom().GetAtomicNum() == b.GetEndAtom().GetAtomicNum() == 54
    }
    cage_degree = {n: sum(1 for b in bonds if n in b) for n in nodes}
    sp3 = {n for places, _ in parts for n in places if cage_degree[n] == 3}
    pi = [n for n in nodes if n not in inserted_names and replaced.get(n, 6) in (6, 14) and n not in set(saturated) | sp3]
    pi_graph = nx.Graph()
    pi_graph.add_nodes_from(pi)
    pi_graph.add_edges_from(tuple(b) for b in bonds if b <= set(pi) and b not in cage_bonds)
    double = {frozenset(pair) for pair in nx.max_weight_matching(pi_graph, maxcardinality=True)}
    mol = Chem.RWMol()
    index = {}
    for n in nodes:
        index[n] = mol.AddAtom(Chem.Atom(inserted_element if n in inserted_names else replaced.get(n, 6)))
    for bond in bonds:
        a, b = tuple(bond)
        mol.AddBond(index[a], index[b], Chem.BondType.DOUBLE if bond in double else Chem.BondType.SINGLE)
    for places, part in parts:
        placed = {}
        for part_atom in part.GetAtoms():
            if part_atom.GetAtomicNum() == 54:
                placed[part_atom.GetIdx()] = index[places[part_atom.GetAtomMapNum() - 1]]
            else:
                added = Chem.Atom(part_atom.GetAtomicNum())
                added.SetNoImplicit(True)
                added.SetNumExplicitHs(part_atom.GetTotalNumHs())
                placed[part_atom.GetIdx()] = mol.AddAtom(added)
        for part_bond in part.GetBonds():
            u, v = part_bond.GetBeginAtomIdx(), part_bond.GetEndAtomIdx()
            if not (part.GetAtomWithIdx(u).GetAtomicNum() == 54 and part.GetAtomWithIdx(v).GetAtomicNum() == 54):
                mol.AddBond(placed[u], placed[v], part_bond.GetBondType())
    for n in nodes:
        atom = mol.GetAtomWithIdx(index[n])
        atom.SetNoImplicit(True)
        if n in inserted_names:
            atom.SetNumExplicitHs({6: 4, 7: 3, 8: 2}[inserted_element] - atom.GetDegree())
        elif atom.GetAtomicNum() in (6, 14):
            atom.SetNumExplicitHs(max(0, (4 if n in saturated else 3) - atom.GetDegree()))
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    return Chem.MolToSmiles(result)


@pytest.mark.parametrize(
    "operations,expected",
    [
        pytest.param({"deleted": (1, 9)}, "1,9-dinor(C60-Ih)[5,6]fullerene", id="nor_even_number_of_atoms"),
        pytest.param({"deleted": (1,), "saturated": (2,)}, "2H-1-nor(C60-Ih)[5,6]fullerene", id="nor_odd_number_gives_indicated_hydrogen"),
        pytest.param({"cleaved": ((1, 9),)}, "1,9-seco(C60-Ih)[5,6]fullerene", id="seco"),
        pytest.param({"cleaved": ((1, 9), (2, 3), (10, 11))}, "1,9:2,3:10,11-triseco(C60-Ih)[5,6]fullerene", id="several_seco_pairs"),
        pytest.param({"inserted": ((1, 9),)}, "1(9)aH-1(9)a-homo(C60-Ih)[5,6]fullerene", id="homo"),
        pytest.param({"inserted": ((1, 9), (2, 3))}, "1(2)aH,3(15)aH-1(2)a,3(15)a-dihomo(C60-Ih)[5,6]fullerene", id="two_homo_atoms_each_with_indicated_hydrogen"),
        pytest.param({"deleted": (1,), "formed": ((2, 9),), "saturated": (2,)}, "2H-2,9-cyclo-1-nor(C60-Ih)[5,6]fullerene", id="cyclo_with_nor"),
        pytest.param({"replaced": {1: 7}, "saturated": (2,)}, "2H-1-aza(C60-Ih)[5,6]fullerene", id="trivalent_heteroatom_gives_indicated_hydrogen"),
        pytest.param({"replaced": {1: 7, 2: 7}}, "1,2-diaza(C60-Ih)[5,6]fullerene", id="two_replacements"),
        pytest.param({"replaced": {1: 7, 9: 14}, "saturated": (2,)}, "2H-1-aza-9-sila(C60-Ih)[5,6]fullerene", id="replacement_prefixes_in_order_of_seniority"),
        pytest.param({"inserted": ((1, 9),), "inserted_element": 7}, "1aH-1a-aza-1(9)a-homo(C60-Ih)[5,6]fullerene", id="replaced_homo_atom"),
        pytest.param({"inserted": ((1, 9),), "inserted_element": 8}, "1a-oxa-1(9)a-homo(C60-Ih)[5,6]fullerene", id="homo_atom_without_hydrogen"),
        pytest.param({"deleted": (1,), "saturated": (2,), "cage": 70}, "2H-1-nor(C70-D5h(6))[5,6]fullerene", id="c70_cage"),
        pytest.param({"deleted": (1, 9, 32, 33, 38, 39)}, "1,9,32,33,38,39-hexanor(C60-Ih)[5,6]fullerene", id="fragment_with_more_than_half_the_atoms"),
        pytest.param({"deleted": (1, 2, 3, 9, 12, 15), "cleaved": ((6, 7),)}, "6,7-seco-1,2,3,9,12,15-hexanor(C60-Ih)[5,6]fullerene", id="seco_and_nor_in_order_of_citation"),
        pytest.param(
            {"cleaved": ((1, 2), (5, 6), (11, 12), (20, 21), (29, 30), (40, 41), (49, 50), (59, 60)), "cage": 70},
            "1,2:5,6:11,12:20,21:29,30:40,41:49,50:59,60-octaseco(C70-D5h(6))[5,6]fullerene",
            id="seco_pairs_of_the_c70_cage",
        ),
    ],
)
def test_structurally_modified_and_replaced_fullerenes(operations, expected):
    operations = dict(operations)
    cage = _FULLERENE_C70_SMILES if operations.pop("cage", 60) == 70 else _FULLERENE_C60_SMILES
    assert smiles_to_iupac(_modified_fullerene(cage, **operations)) == expected


_CYCLOPROPANE, _OXIRENE, _DIOXOLE = "C1[Xe:1][Xe:2]1", "O1[Xe:1][Xe:2]1", "C1O[Xe:1][Xe:2]O1"


@pytest.mark.parametrize(
    "operations,expected",
    [
        # P-27.6.1: every bond of a cage is named by its lowest locants, here a bond that is not numbered 1,9
        pytest.param({"attached": [((3, 15), _CYCLOPROPANE)]}, "3′H-cyclopropa[1,9](C60-Ih)[5,6]fullerene", id="hydrocarbon_ring_omits_its_locants"),
        pytest.param({"attached": [((1, 9), _DIOXOLE)]}, "2′H-[1,3]dioxolo[4′,5′:1,9](C60-Ih)[5,6]fullerene", id="heterocycle_cites_primed_locants"),
        pytest.param(
            {"attached": [((1, 9), _OXIRENE), ((2, 12), _OXIRENE), ((7, 8), _OXIRENE)]},
            "tris(oxireno)[2′,3′:1,9;2″,3″:2,12;2‴,3‴:7,8](C60-Ih)[5,6]fullerene",
            id="identical_rings_in_the_order_of_the_lower_fullerene_locant",
        ),
        pytest.param(
            {"attached": [((7, 22), _CYCLOPROPANE), ((33, 34), _CYCLOPROPANE), ((46, 47), _CYCLOPROPANE)], "cage": 70},
            "3′H,3″H,3‴H-tricyclopropa[7,22:33,34:46,47](C70-D5h(6))[5,6]fullerene",
            id="hydrocarbon_rings_separated_by_colons",
        ),
        pytest.param(
            {"attached": [((2, 12), _OXIRENE), ((1, 9), _CYCLOPROPANE)]},
            "3′H-cyclopropa[1,9]oxireno[2″,3″:2,12](C60-Ih)[5,6]fullerene",
            id="different_rings_primed_in_alphanumerical_order",
        ),
        pytest.param(
            {"attached": [((8, 25), "C1=CS[Xe:2][Xe:1]1")], "cage": 70},
            "thieno[3′,2′:8,25](C70-D5h(6))[5,6]fullerene",
            id="ring_locants_follow_the_direction_of_the_fullerene_locants",
        ),
        pytest.param({"attached": [((1, 9), "C1CC[Xe:1][Xe:2]1")]}, "4′,5′-dihydro-3′H-cyclopenta[1,9](C60-Ih)[5,6]fullerene", id="hydro_prefixes_with_indicated_hydrogen"),
        pytest.param({"attached": [((1, 9), "C1C=CC[Xe:1][Xe:2]1")]}, "3′,6′-dihydrobenzo[1,9](C60-Ih)[5,6]fullerene", id="hydro_prefixes_without_indicated_hydrogen"),
        pytest.param({"attached": [((1, 9), "C1=C2C=CC=CC2=C[Xe:1][Xe:2]1")]}, "naphtho[2′,3′:1,9](C60-Ih)[5,6]fullerene", id="ring_system_of_two_rings"),
        pytest.param({"attached": [((1, 9), "N1[Xe:2][Xe:1]C2=CC=CC=C12")]}, "1′H-indolo[2′,3′:1,9](C60-Ih)[5,6]fullerene", id="heterocyclic_ring_system"),
        # P-27.6.1: a ring system with no single fusion component name keeps its own fusion name and numbering as one prefix
        pytest.param({"attached": [((1, 9), "O1[Xe:1][Xe:2]c2c1ccc1ccccc21")]}, "naphtho[2,1-b]furo[1′,2′:1,9](C60-Ih)[5,6]fullerene", id="composite_ring_system_is_one_prefix"),
        pytest.param({"attached": [((1, 9), "C1[Xe:1][Xe:2]c2c1ccc1ccccc21")]}, "3′H-cyclopenta[a]naphtho[1′,2′:1,9](C60-Ih)[5,6]fullerene", id="indicated_hydrogen_in_the_numbering_of_the_whole_system"),
        pytest.param(
            {"deleted": (1, 2, 3, 4, 5, 6, 9, 12, 15, 18), "attached": [((7, 8), "C1C=C[Xe:1][Xe:2]1")]},
            "3′H-cyclopenta[7,8]-1,2,3,4,5,6,9,12,15,18-decanor(C60-Ih)[5,6]fullerene",
            id="ring_on_a_nor_fragment_instead_of_cyclo_and_seco",
        ),
        pytest.param(
            {"deleted": (1,), "saturated": (2,), "attached": [((7, 8), _CYCLOPROPANE)]},
            "2H,3′H-cyclopropa[7,8]-1-nor(C60-Ih)[5,6]fullerene",
            id="indicated_hydrogen_of_the_cage_and_of_the_ring",
        ),
        # P-27.6.1, P-31.1.4.2.4: substituents of a fused ring carry its primed locants; a substituted sp3 position still cites its hydrogen
        pytest.param({"attached": [((1, 9), "C1(C)[Xe:1][Xe:2]1")]}, "3′-methyl-3′H-cyclopropa[1,9](C60-Ih)[5,6]fullerene", id="substituent_on_a_hydrocarbon_ring"),
        pytest.param({"attached": [((1, 9), "CC1(C)O[Xe:1][Xe:2]O1")]}, "2′,2′-dimethyl-2′H-[1,3]dioxolo[4′,5′:1,9](C60-Ih)[5,6]fullerene", id="two_substituents_on_one_position"),
        pytest.param({"attached": [((1, 9), "ClC1=CS[Xe:2][Xe:1]1")]}, "4′-chlorothieno[2′,3′:1,9](C60-Ih)[5,6]fullerene", id="lowest_locant_among_the_numberings_of_the_ring"),
        # P-27.6.1, P-25.3.1.3: an alicyclic bridged ring system stays saturated and is cited by its von Baeyer name ending in 'ano'
        pytest.param({"attached": [((1, 9), "[Xe:1]1[Xe:2]C2CCC1C2")]}, "bicyclo[2.2.1]heptano[2′,3′:1,9](C60-Ih)[5,6]fullerene", id="bridged_alicyclic_ring_system"),
        pytest.param({"attached": [((1, 9), "[Xe:1]1[Xe:2]C2CCC1O2")]}, "7-oxabicyclo[2.2.1]heptano[2′,3′:1,9](C60-Ih)[5,6]fullerene", id="bridged_ring_system_with_a_replacement_prefix"),
        pytest.param({"attached": [((1, 9), "[Xe:1]1[Xe:2]C2C=CC1C2")]}, "bicyclo[2.2.1]hept-2-eno[5′,6′:1,9](C60-Ih)[5,6]fullerene", id="double_bond_numbered_before_the_fusion_locants"),
        pytest.param(
            {"attached": [((1, 9), "[Xe:1]1[Xe:2]C2CCC1O2"), ((2, 12), "[Xe:1]1[Xe:2]C2CCC1C2")]},
            "bicyclo[2.2.1]heptano[2′,3′:1,9]-7-oxabicyclo[2.2.1]heptano[2″,3″:2,12](C60-Ih)[5,6]fullerene",
            id="prefixes_in_alphanumerical_order_ignoring_locants",
        ),
    ],
)
def test_fullerenes_with_fused_rings(operations, expected):
    operations = dict(operations)
    cage = _FULLERENE_C70_SMILES if operations.pop("cage", 60) == 70 else _FULLERENE_C60_SMILES
    assert smiles_to_iupac(_modified_fullerene(cage, **operations)) == expected


def test_bridged_alicyclic_ring_on_a_fullerene_has_no_book_prefix_form():
    with pytest.warns(NonPreferredNameWarning, match="P-27.6.1"):
        smiles_to_iupac(_modified_fullerene(_FULLERENE_C60_SMILES, attached=[((1, 9), "[Xe:1]1[Xe:2]C2CCC1C2")]))


_ETHANO, _ETHENO, _EPOXY = "[Xe:1]CC[Xe:2]", "[Xe:1]C=C[Xe:2]", "[Xe:1]O[Xe:2]"


@pytest.mark.parametrize(
    "operations,expected",
    [
        pytest.param({"attached": [((1, 4), _ETHANO)], "cage": 70}, "1,4-ethano(C70-D5h(6))[5,6]fullerene", id="bridge_between_nonadjacent_atoms"),
        pytest.param(
            {"attached": [((1, 4), "[Xe:1]C1=CC=CC=C1[Xe:2]")], "cage": 70},
            "1,4-[1,2]benzeno(C70-D5h(6))[5,6]fullerene",
            id="bridge_that_is_a_ring_with_its_locants",
        ),
        pytest.param(
            {"attached": [((1, 7), _ETHANO), ((3, 13), _ETHANO)]},
            "1,7:3,13-diethano(C60-Ih)[5,6]fullerene",
            id="identical_bridges_multiplied",
        ),
        pytest.param(
            {"attached": [((1, 7), _ETHANO), ((3, 13), _ETHENO), ((5, 8), _EPOXY)]},
            "1,11-epoxy-2,10-ethano-4,16-etheno(C60-Ih)[5,6]fullerene",
            id="different_bridges_in_alphabetical_order_with_the_lowest_locant_set",
        ),
        pytest.param(
            {"attached": [((1, 9, 32, 33), "C1(COCC2=CC=C(C=C2)COCC3[Xe:3][Xe:4]3)[Xe:1][Xe:2]1")]},
            "3′H,3″H-3′,3″-(methanooxymethano[1,4]benzenomethanooxymethano)dicyclopropa[1,9:32,33](C60-Ih)[5,6]fullerene",
            id="bridge_between_two_fused_rings",
        ),
        pytest.param(
            {"attached": [((1, 9, 31), "C1(CC[Xe:3])[Xe:1][Xe:2]1"), ((2, 12, 32), "C1(CC[Xe:3])[Xe:1][Xe:2]1")]},
            "3′H,3″H-3′,23:3″,24-diethanodicyclopropa[1,9:2,12](C60-Ih)[5,6]fullerene",
            id="bridges_between_fused_rings_and_the_cage",
        ),
    ],
)
def test_bridged_fullerenes(operations, expected):
    operations = dict(operations)
    cage = _FULLERENE_C70_SMILES if operations.pop("cage", 60) == 70 else _FULLERENE_C60_SMILES
    assert smiles_to_iupac(_modified_fullerene(cage, **operations)) == expected


@pytest.mark.parametrize(
    "operations,expected",
    [
        pytest.param(
            {"inserted": ((1, 9),), "attached": [(("inserted0",), "[Xe:1]1CCCCC1")]},
            "spiro[cyclohexane-1,1′a-[1(9)a]homo(C60-Ih)[5,6]fullerene]",
            id="homofullerene_atom_is_the_spiro_atom",
        ),
        pytest.param(
            {"attached": [((1, 9), "C1(CCCCC1)([Xe:1])[Xe:2]")]},
            "spiro[cyclohexane-1,3′′-cyclopropa[1,9](C60-Ih)[5,6]fullerene]",
            id="atom_of_a_fused_ring_is_the_spiro_atom",
        ),
    ],
)
def test_spiro_fullerenes(operations, expected):
    assert smiles_to_iupac(_modified_fullerene(_FULLERENE_C60_SMILES, **operations)) == expected


def test_spiro_atom_of_a_fused_ring_with_other_saturated_positions_is_not_named():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(_modified_fullerene(_FULLERENE_C60_SMILES, attached=[((1, 9), "[Xe:1]CCC2(CCCCC2)[Xe:2]")]))


def test_bridge_with_a_substituent_is_not_named():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(_modified_fullerene(_FULLERENE_C60_SMILES, attached=[((1, 7), "[Xe:1]C(Cl)C[Xe:2]")]))
