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
        ("C1#CC(C)=CC=C1", "3-methyl-1,2-didehydrobenzene"),
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC/C=C/C=C\\CCCC1", "(1Z,3E)-cycloundeca-1,3-diene"),
        ("C1CC/C=C/CCC#CCCCC1", "(1E)-cyclotridec-1-en-5-yne"),
        ("C1=CC=CC=CC=CC=CC=CC=C1", "cyclotetradeca-1,3,5,7,9,11,13-heptaene"),
        ("C1=C\\C=C/C=C\\C=C/C=C\\1", "(1Z,3Z,5Z,7Z,9E)-cyclodeca-1,3,5,7,9-pentaene"),
    ],
)
def test_large_unsaturated_carbocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_double_bond_stereo_partly_unspecified_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC/C=C\\CC=CCCCC1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H]1CCC=CC1", "(4S)-4-ethylcyclohex-1-ene"),
        ("CC[C@H](C)[C@H]1CCC=CC1", "(4S)-4-[(2S)-butan-2-yl]cyclohex-1-ene"),
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
        ("C1C2CC3CC(C2)(CC1C3=O)O", "5-hydroxyadamantan-2-one"),
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


def test_spiro_hydroxyl_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("OCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanol"


def test_spiro_amine_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("NCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanamine"


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


def test_radical_centres_on_substituents_of_a_spiro_ring():
    assert smiles_to_iupac("[CH2]C1CCC2(CC1)CCCC2[CH2]") == "[1-(ylomethyl)spiro[4.5]decan-8-yl]methyl"


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
        ("C=C1CC2CCC1C2", "2-methylidenebicyclo[2.2.1]heptane"),
        ("CC=C1CC2CC1C=C2", "5-ethylidenebicyclo[2.2.1]hept-2-ene"),
        ("C=C1CC2CC3CC1CC2C3", "5-methylidenetricyclo[4.3.1.0^3,8]decane"),
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
    ],
)
def test_monospiro_union_with_a_von_baeyer_component(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CC2(C=CNC=C2)C2=CC=CC=C12", "2,3-dihydro-1′H-spiro[indene-1,4′-pyridine]", id="mancude_retained_monocycle_with_indicated_hydrogen"),
        pytest.param("C1CC2(CCN(C=C2)C)c2ccccc12", "1′-methyl-2,2′,3,3′-tetrahydro-1′H-spiro[indene-1,4′-pyridine]", id="hydro_prefixes_of_a_partly_saturated_monocycle"),
        pytest.param("C12C[N+]3(CC(CC1)CC2)COC=C3", "2′H-3λ5-spiro[3-azabicyclo[3.2.2]nonane-3,3′-[1,3]oxazol]-3-ylium", id="hantzsch_widman_mancude_monocycle_with_cationic_spiro_atom"),
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1C2(Cc3ccccc3C24CCCC4)CC3(Cc4ccccc4C3)CC15Cc6ccccc6C5", "1′′′,1′′′′,3′′′,3′′′′-tetrahydro-3′′H-spiro{cyclopentane-1,1′′-trispiro[cyclohexane-1,2′:3,2′′:5,2′′′-tris(indene)]}", id="largest_branched_system_named_first_and_used_as_a_unit"),
        pytest.param("c1ccc2c(c1)CC1(C2)OP23(S1)SP14(OC5(Cc6ccccc6C5)S1)SP1(OC5(Cc6ccccc6C5)S1)(S2)SP1(OC2(Cc5ccccc5C2)S1)(S3)S4", "1′′′′′,1′′′′′′,1′′′′′′′,1′′′′′′′′,3′′′′′,3′′′′′′,3′′′′′′′,3′′′′′′′′-octahydro-1λ5,2′′λ5,2′′′λ5,2′′′′λ5-tetraspiro{tetraspiro[2,4,6,8,9,10-hexathia-1,3,5,7-tetraphosphaadamantane-1,2′:3,2′′:5,2′′′:7,2′′′′-tetrakis([1,3,2]oxathiaphosphetane)]-4′,2′′′′′:4′′,2′′′′′′:4′′′,2′′′′′′′:4′′′′,2′′′′′′′′-tetrakis(indene)}", id="heteroadamantane_unit_with_nonstandard_spiro_atoms"),
        pytest.param("C1C2(Cc3ccccc3C24CCCC4)CC3(Cc4ccccc4C3)CC15Cc6ccccc6C5=O", "1′′′′,3′′′′-dihydro-3′′H-spiro{cyclopentane-1,1′′-trispiro[cyclohexane-1,2′:3,2′′:5,2′′′-tris(indene)]}-1′′′(3′′′H)-one", id="suffix_group_on_a_nested_system"),
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
        pytest.param("[SiH]12O[SiH]3O[SiH](O1)O[SiH](O2)O3", "tricyclo[3.3.1.1^3,7]tetrasiloxane", id="adamantane_shaped_siloxane"),
        pytest.param("O1[SiH]2O[SiH]3O[SiH]1O[SiH]1O[SiH](O[SiH](O1)O3)O2", "tetracyclo[5.5.1.1^3,11.1^5,9]hexasiloxane", id="hexasilasesquioxane"),
        pytest.param("[Si]12(O[SiH2]O[SiH2]O1)O[SiH2]O[SiH2]O2", "spiro[5.5]pentasiloxane", id="spiro_siloxane"),
        pytest.param("N12[SiH2]N3[SiH2]N([SiH2]1)[SiH2]N([SiH2]2)[SiH2]3", "1N-tricyclo[3.3.1.1^3,7]hexasilazane", id="bridgehead_nitrogen_named_by_prefix"),
    ],
)
def test_alternating_heteroatom_cages_take_the_preselected_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
