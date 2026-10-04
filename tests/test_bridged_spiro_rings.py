import pytest
from rdkit import Chem
from smiles_to_iupac import smiles_to_iupac
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


def test_asymmetric_1_2_disubstituted_stereo():
    assert smiles_to_iupac("C[C@@H]1CCCC[C@H]1CC") == "(1R,2R)-1-ethyl-2-methylcyclohexane"


def test_1_3_disubstituted_stereo():
    assert smiles_to_iupac("C[C@@H]1C[C@H](C)CCC1") == "(1R,3S)-1,3-dimethylcyclohexane"


def test_single_ring_stereocenter_cites_the_specified_elements():
    assert smiles_to_iupac("C[C@H]1CCCCC1CC") == "(2S)-1-ethyl-2-methylcyclohexane"


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


def test_substituted_tetrabenzenacyclooctaphane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cccc2c1CC1=CC=CC(=C1)CC1=CC=CC(=C1)CC1=CC=CC(=C1)C2")


def test_substituted_metacyclophane_is_named():
    assert (
        smiles_to_iupac("CC1CC2=CC(=CC=C2)CCC3=CC=CC1=C3")
        == "2-methyl-1,4(1,3)-dibenzenacyclohexaphane"
    )


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Oc1cc2ccc1CCc1ccc(cc1)CC2", "1,4(1,4)-dibenzenacyclohexaphan-12-ol"),
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
            "2-[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-2-ylidene]ethanoic acid",
        ),
        (
            "OC(=O)C=C1C=C2Cc3cccc(n3)Cc3cccc(n3)Cc3cccc(n3)CC(=C1)N2",
            "2-[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-14(11H)-ylidene]ethanoic acid",
        ),
    ],
)
def test_pyridine_amplificants_and_ylidene_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_extra_unsaturation_outside_aromatic_ring_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC2C(C)C1c1ccccc12")


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


def test_unsaturated_spiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC12CCCC=C2")


def test_multiple_spiro_hydroxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCC12CCCCC2O")


def test_spiro_hydroxyl_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("OCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanol"


def test_unsaturated_spiro_alcohol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCC12C=CCCC2")


def test_stereocenter_alongside_spiro_alcohol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[C@H]1CCCC12CCCCC2")


def test_multiple_spiro_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC2(C1)CCCCC2N")


def test_spiro_amine_on_substituent_branch_is_a_prefix():
    assert smiles_to_iupac("NCC1CCCC12CCCCC2") == "(spiro[4.5]decan-1-yl)methanamine"


def test_unsaturated_spiro_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC12C=CCCC2")


def test_multiple_spiro_carbenium_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CCC2(CC1)CCCC2[CH2+]")


def test_spiro_carbenium_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CCC2(CC1)CCCC2")


def test_unsaturated_heteroatom_spiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2(OC=CCC2)C1")


def test_oxaspiro_other_configuration():
    assert smiles_to_iupac("C[C@@H]1CCCC12CCOCC2") == "(1R)-1-methyl-8-oxaspiro[4.5]decane"


def test_plain_oxaspiro_without_stereo_still_resolves():
    assert smiles_to_iupac("C1CCCC12CCOCC2") == "8-oxaspiro[4.5]decane"


def test_multiple_spiro_ketones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC2(C1)CCCC(=O)C2")


def test_unsaturated_spiro_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC12C=CCCC2")


def test_multiple_spiro_radical_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C1CCC2(CC1)CCCC2[CH2]")


def test_spiro_radical_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C1CCC2(CC1)CCCC2")


def test_unsaturated_spiro_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH]1CCC2(C=CC2)CC1")


def test_stereo_spiro_radical_resolves():
    assert smiles_to_iupac("[CH]1CCC2(CC1)CCC[C@@H](C)C2") == "(8R)-8-methylspiro[5.5]undecan-3-yl"


def test_single_off_spiro_stereocenter():
    assert smiles_to_iupac("C[C@H]1CCC[C@]2(C1)CCCC2") == "(7S)-7-methylspiro[4.5]decane"


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


def test_multiple_ring_carbenium_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CC2CCC1C2[CH2+]")


def test_carbenium_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CC2CCC1C2")


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


def test_two_different_ring_heteroatoms_tricyclic():
    assert smiles_to_iupac("O1C2CC3CC1CC(C2)N3") == "2-oxa-6-azatricyclo[3.3.1.1^3,7]decane"


def test_two_same_element_ring_heteroatoms_tricyclic():
    assert smiles_to_iupac("O1C2CC3OC1CC(C2)C3") == "2,4-dioxatricyclo[3.3.1.1^3,7]decane"


def test_disjoint_ring_systems_joined_by_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1(C)C2CC[C@@]1(C)CN(CCC1C3CC4CC(C3)CC1C4)C2")


def test_unsaturated_von_baeyer_ketone_now_resolves():
    assert smiles_to_iupac("O=C1CC2C=CC1C2") == "bicyclo[2.2.1]hept-5-en-2-one"


def test_von_baeyer_ketone_stereocenter_with_coexisting_hydroxyl_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1C[C@H]2CC[C@H]1[C@@H]2O")


def test_von_baeyer_ketone_unsaturation_with_coexisting_hydroxyl_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2C=CC1C2O")


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


def test_multiple_ring_selenols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CC2CCC1C([SeH])C2")


def test_unsaturated_monospiro_selenol_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CCCC2(C1)C=CCCC2")


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


def test_n_substituted_sulfinamide_on_polycyclic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNS(=O)C1CC2CCC1C2")


def test_unsaturated_monospiro_sulfinamide_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)C1CCCC2(C1)C=CCCC2")


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


def test_n_substituted_sulfonamide_on_polycyclic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNS(=O)(=O)C1CC2CCC1C2")


def test_multiple_ring_sulfonamides_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)C1CC2CCC1C2S(N)(=O)=O")


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


def test_unsaturated_monospiro_sulfonic_acid_still_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC2(C1)C=CCCC2")


def test_multiple_ring_tellones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CC2CCC1C(=[Te])C2")


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


def test_substituted_fullerene_stops_at_the_numbering_the_blue_book_leaves_open():
    from rdkit import Chem
    from smiles_to_iupac._common import UnsupportedStructure
    from smiles_to_iupac._fullerene import _FULLERENE_C60_SMILES

    cage = Chem.MolFromSmiles(_FULLERENE_C60_SMILES)
    Chem.Kekulize(cage, clearAromaticFlags=True)
    editable = Chem.RWMol(cage)
    next(b for b in editable.GetAtomWithIdx(0).GetBonds() if b.GetBondTypeAsDouble() == 2.0).SetBondType(Chem.BondType.SINGLE)
    carbon = editable.AddAtom(Chem.Atom(6))
    editable.AddBond(0, carbon, Chem.BondType.SINGLE)
    derivative = editable.GetMol()
    Chem.SanitizeMol(derivative)
    with pytest.raises(UnsupportedStructure, match="P-27.3"):
        smiles_to_iupac(Chem.MolToSmiles(derivative))
