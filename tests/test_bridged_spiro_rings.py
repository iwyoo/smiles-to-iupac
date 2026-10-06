import pytest
from rdkit import Chem
from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._fullerene import (
    _C84_D6H_ISOMER_24_SMILES,
    _FULLERENE_C60_SMILES,
    _FULLERENE_C70_SMILES,
)
from smiles_to_iupac._fullerene_spiral import (
    _C84_IPR_ISOMERS,
    _canonical_spiral,
    _pentagon_positions,
    _planar_faces,
)


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
        ("C[C@@H]1CCCC[C@@H]1C", "cis-1,2-dimethylcyclohexane"),
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
    ],
)
def test_unsubstituted_cyclic_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_enyne_ring():
    assert smiles_to_iupac("C1=CCC#CCCCCCCCCCC1") == "cyclopentadec-1-en-4-yne"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC/C=C/CC1", "(E)-cyclooctene"),
    ],
)
def test_ring_double_bond_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_double_bond_stereo_multiple_bonds_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC/C=C\\C/C=C\\C1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H]1CCC=CC1", "(4S)-4-ethylcyclohexene"),
        ("CC[C@H](C)[C@H]1CCC=CC1", "(4S)-4-[(2S)-butan-2-yl]cyclohexene"),
    ],
)
def test_ring_tetrahedral_stereocenter(smiles, expected):
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


def test_c70_fullerene():
    assert smiles_to_iupac(_FULLERENE_C70_SMILES) == "(C70-D5h(6))[5,6]fullerene"


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


def test_heptacyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C15C4C14C3C2C54")


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
        ("C1C2CC3CC(C2)(CC1C3=O)O", "5-hydroxytricyclo[3.3.1.1^3,7]decan-2-one"),
    ],
)
def test_hydroxy_von_baeyer_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_pentacyclic_propellane_like_degree_four_branch_atom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C15C2C31C45")


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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC12CCC3(CC2)CCC3", "dispiro[2.2.3^6.2^3]undecane"),
        ("C1CCCC12CC3(CC2)CCCC3", "dispiro[4.1.4^7.2^5]tridecane"),
    ],
)
def test_smiles_to_iupac_linear_polyspiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_linear_polyspiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC12CCC3(CC2)CCC3")


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


def test_unsaturated_heteroatom_polyspiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2(C1)C=CC3(O2)CCCC3")


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


def test_substituent_branch_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("F[C@@H](Cl)C1CCC12CCC1(CCC1)CC2")


def test_fused_bicyclic_is_not_spiro():
    assert smiles_to_iupac("C1CCC2CCCCC2C1") == "decahydronaphthalene"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C1CC12CCCC=C2", id="unsaturated_spiro_raises"),
        pytest.param("OC1CCCC12CCCCC2O", id="multiple_spiro_hydroxyls_raises"),
    ],
)
def test_unsaturated_spiro_raises_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_spiro_hydroxyl_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("OCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanol"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OC1CCCC12C=CCCC2", id="unsaturated_spiro_alcohol_raises"),
        pytest.param("O[C@H]1CCCC12CCCCC2", id="stereocenter_alongside_spiro_alcohol_raises"),
        pytest.param("NC1CCCC2(C1)CCCCC2N", id="multiple_spiro_amines_raises"),
    ],
)
def test_unsaturated_spiro_alcohol_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_spiro_amine_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("NCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanamine"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NC1CCCC12C=CCCC2", id="unsaturated_spiro_amine_raises"),
        pytest.param("[CH2+]C1CCC2(CC1)CCCC2[CH2+]", id="multiple_spiro_carbenium_centers_raises"),
        pytest.param("[CH2+]C1CCC2(CC1)CCCC2", id="spiro_carbenium_on_substituent_branch_raises"),
        pytest.param("C1CCC2(OC=CCC2)C1", id="unsaturated_heteroatom_spiro_raises"),
    ],
)
def test_unsaturated_spiro_amine_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    "smiles",
    [
        pytest.param("O=C1CCCC2(C1)CCCC(=O)C2", id="multiple_spiro_ketones_raises"),
        pytest.param("O=C1CCCC12C=CCCC2", id="unsaturated_spiro_ketone_raises"),
        pytest.param("[CH2]C1CCC2(CC1)CCCC2[CH2]", id="multiple_spiro_radical_centers_raises"),
        pytest.param("[CH2]C1CCC2(CC1)CCCC2", id="spiro_radical_on_substituent_branch_raises"),
        pytest.param("[CH]1CCC2(C=CC2)CC1", id="unsaturated_spiro_radical_raises"),
    ],
)
def test_multiple_spiro_ketones_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH]1CCC2(CC1)CCC[C@@H](C)C2", "(8R)-8-methylspiro[5.5]undecan-3-yl", id="stereo_spiro_radical_resolves"),
        pytest.param("C[C@H]1CCC[C@]2(C1)CCCC2", "(7S)-7-methylspiro[4.5]decane", id="single_off_spiro_stereocenter"),
    ],
)
def test_stereo_spiro_radical_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_spiro_thiols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC2(C1)CCCCC2S")


def test_spiro_thiol_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("SCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanethiol"


def test_unsaturated_spiro_thiol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC12C=CCCC2")


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


def test_multiple_ring_hydroxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CC2CCC1C2O")


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
    "smiles",
    [
        pytest.param("[CH2+]C1CC2CCC1C2[CH2+]", id="multiple_ring_carbenium_centers_raises"),
        pytest.param("[CH2+]C1CC2CCC1C2", id="carbenium_on_substituent_branch_raises"),
    ],
)
def test_multiple_ring_carbenium_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


def test_three_mixed_element_heteroatoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2NNC1CO2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O1C2CC3CC1CC(Cl)(C2)C3", "5-chloro-2-oxatricyclo[3.3.1.1^3,7]decane"),
    ],
)
def test_smiles_to_iupac_von_baeyer_heteroatom_tricyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O1C2CC3CC1CC(C2)N3", "2-oxa-6-azatricyclo[3.3.1.1^3,7]decane", id="different_ring_heteroatoms_tricyclic"),
        pytest.param("O1C2CC3OC1CC(C2)C3", "2,4-dioxatricyclo[3.3.1.1^3,7]decane", id="same_element_ring_heteroatoms_tricyclic"),
    ],
)
def test_two_different_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_disjoint_ring_systems_joined_by_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1(C)C2CC[C@@]1(C)CN(CCC1C3CC4CC(C3)CC1C4)C2")


def test_unsaturated_von_baeyer_ketone_now_resolves():
    assert smiles_to_iupac("O=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-en-2-one"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("O=C1C[C@H]2CC[C@H]1[C@@H]2O", id="stereocenter_with_coexisting_hydroxyl_still_raises"),
        pytest.param("O=C1CC2C=CC1C2O", id="unsaturation_with_coexisting_hydroxyl_still_raises"),
    ],
)
def test_von_baeyer_ketone_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[C]12CC3CC(CC(C3)C1)C2", "adamantan-1-yl"),
    ],
)
def test_von_baeyer_radical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_radical_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C1CC2CCC1C2")


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
    "smiles",
    [
        pytest.param("[SeH]C1CC2CCC1C([SeH])C2", id="multiple_ring_selenols_raises"),
        pytest.param("[SeH]C1CCCC2(C1)C=CCCC2", id="unsaturated_monospiro_selenol_still_raises"),
    ],
)
def test_multiple_ring_selenols_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Se]=C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-selone"),
    ],
)
def test_von_baeyer_spiro_selone_reviewed(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_selones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CC2CCC1C(=[Se])C2")


def test_von_baeyer_selone_ring_unsaturation():
    assert smiles_to_iupac("[Se]=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-selone"


def test_unsaturated_monospiro_selone_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCCC2(C1)C=CCCC2")


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
    "smiles",
    [
        pytest.param("CNS(=O)C1CC2CCC1C2", id="n_substituted_sulfinamide_on_polycyclic_ring_raises"),
        pytest.param("NS(=O)C1CCCC2(C1)C=CCCC2", id="unsaturated_monospiro_sulfinamide_still_raises"),
    ],
)
def test_n_substituted_sulfinamide_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


def test_unsaturated_monospiro_sulfinic_acid_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CCCC2(C1)C=CCCC2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NS(=O)(=O)C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-sulfonamide"),
    ],
)
def test_von_baeyer_spiro_sulfonamide_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CNS(=O)(=O)C1CC2CCC1C2", id="n_substituted_sulfonamide_on_polycyclic_ring_raises"),
        pytest.param("NS(=O)(=O)C1CC2CCC1C2S(N)(=O)=O", id="multiple_ring_sulfonamides_raises"),
    ],
)
def test_n_substituted_sulfonamide_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_von_baeyer_sulfonamide_ring_unsaturation():
    assert smiles_to_iupac("NS(=O)(=O)C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-sulfonamide"


def test_unsaturated_monospiro_sulfonamide_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)C1CCCC2(C1)C=CCCC2")


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


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OS(=O)(=O)C1CCCC2(C1)C=CCCC2", id="unsaturated_monospiro_sulfonic_acid_still_raises"),
        pytest.param("[Te]=C1CC2CCC1C(=[Te])C2", id="multiple_ring_tellones_raises"),
    ],
)
def test_unsaturated_monospiro_sulfonic_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_von_baeyer_tellone_ring_unsaturation():
    assert smiles_to_iupac("[Te]=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-tellone"


def test_unsaturated_monospiro_tellone_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC2(C1)C=CCCC2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[TeH]C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-tellurol"),
    ],
)
def test_von_baeyer_spiro_tellurol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_tellurols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CC2CCC1C2[TeH]")


def test_von_baeyer_tellurol_ring_unsaturation():
    assert smiles_to_iupac("[TeH]C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-tellurol"


def test_unsaturated_monospiro_tellurol_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CCCC2(C1)C=CCCC2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("S=C1CCC2(CC1)CCCC2", "spiro[4.5]decane-8-thione"),
    ],
)
def test_von_baeyer_spiro_thione_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_thiones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CC2CCC1C(=S)C2")


def test_von_baeyer_thione_ring_unsaturation():
    assert smiles_to_iupac("S=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-ene-2-thione"


def test_unsaturated_monospiro_thione_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CCCC2(C1)C=CCCC2")


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


def test_fullerene_with_characteristic_group_is_named_as_the_parent_with_added_hydrogen():
    # P-6 (P-58.2): the principal characteristic group takes the added hydrogen, not hydro prefixes
    assert smiles_to_iupac(_fullerene_derivative(_FULLERENE_C60_SMILES, {1: "O", 9: None})) == (
        "(C60-Ih)[5,6]fulleren-1(9H)-ol"
    )


def test_bridged_fused_naming_does_not_write_rdkit_logs_to_stderr(capfd):
    smiles_to_iupac(_fullerene_derivative(_FULLERENE_C60_SMILES, {1: _ACETIC, 9: None}))
    assert capfd.readouterr().err == ""
