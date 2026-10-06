import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClP", "chlorophosphane"),
    ],
)
def test_smiles_to_iupac_simple_phosphane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituted_alkyl_chain():
    assert smiles_to_iupac("ClCCP") == "(2-chloroethyl)phosphane"


def test_multiplied_compound_substituent_with_different_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCP(C(C)C)C(C)C")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=CP", "ethenylphosphane", id="unsaturated_substituent"),
        pytest.param("C1CCCCC1P", "cyclohexylphosphane", id="non_aromatic_ring"),
        pytest.param("CP(c1ccc(Cl)cc1)", "(4-chlorophenyl)(methyl)phosphane", id="halogenated_phenyl_mixed_with_alkyl"),
    ],
)
def test_unsaturated_substituent_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("PP", "diphosphane"),
    ],
)
def test_phosphane_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_phosphane_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("P(P)(P)P")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("P1PPPP1", "2,3,4,5-tetrahydro-1H-pentaphosphole", id="cyclic_phosphane_chain"),
        pytest.param("CPP", "methyldiphosphane", id="carbon_phosphorus_mix_is_a_diphosphane"),
    ],
)
def test_cyclic_phosphane_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_zero_substituents_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=P")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=Pc1ccccc1", "phenylphosphanone", id="single_phenyl_substituent"),
        pytest.param("CP(=O)c1ccccc1", "methyl(phenyl)-λ5-phosphanone", id="mixed_alkyl_and_phenyl_substituents"),
    ],
)
def test_single_phenyl_substituent_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("O=Pc1ccccc1C", id="substituted_phenyl_raises"),
        pytest.param("CC(C)P=O", id="branched_substituent_raises"),
    ],
)
def test_substituted_phenyl_raises_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_ring_phosphorus_oxide_is_a_lambda5_heterone():
    assert smiles_to_iupac("O=P1CCCCC1") == "1λ5-phosphinan-1-one"


@pytest.mark.parametrize(
    "smiles, name",
    [
        ("O=P1(O)OCCO1", "2-hydroxy-1,3,2λ5-dioxaphospholan-2-one"),
        ("O=P1(OC)OC(C)CO1", "2-methoxy-4-methyl-1,3,2λ5-dioxaphospholan-2-one"),
        ("OP1(=O)OCCCO1", "2-hydroxy-1,3,2λ5-dioxaphosphinan-2-one"),
        ("S=P1(O)OCCO1", "2-hydroxy-1,3,2λ5-dioxaphospholane-2-thione"),
    ],
)
def test_cyclic_phosphate_diester_is_a_heterocycle(smiles, name):
    assert smiles_to_iupac(smiles) == name


def test_salt_of_partial_ester():
    assert smiles_to_iupac("COP(=O)(O)[O-].[Na+]") == "sodium methyl hydrogen phosphate"
    assert smiles_to_iupac("CCOP(=O)(O)[O-].[Na+]") == "sodium ethyl hydrogen phosphate"
    assert smiles_to_iupac("COP(=O)(O)[O-].[K+]") == "potassium methyl hydrogen phosphate"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("COP(=O)(O)[O-].[Ca+2]", id="salt_of_partial_ester_multivalent_cation_raises"),
        pytest.param("OP(=O)(O)O", id="phosphoric_acid_itself_unaffected"),
    ],
)
def test_salt_of_partial_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_phosphindole():
    assert smiles_to_iupac("C1C=C2C=CC=CC2=P1") == "2H-phosphindole"


def test_asymmetric_substituents():
    assert smiles_to_iupac("CCP(C)(=O)O") == "ethyl(methyl)phosphinic acid"
    assert smiles_to_iupac("CC(C)P(C)(=O)O") == "methyl(propan-2-yl)phosphinic acid"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OP(C)(=O)CCP(C)(=O)O", id="second_phosphinic_acid_group"),
        pytest.param("NCP(C)(=O)O", id="unrecognized_heteroatom"),
    ],
)
def test_rejects_second_phosphinic_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1=CC=PC=C1", "phosphinine", id="phosphinine"),
        pytest.param("CC1=CC=CC=P1", "2-methylphosphinine", id="methylphosphinine"),
    ],
)
def test_phosphinine_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COP(OC)OC", "trimethyl phosphite"),
        ("COP(OC)O", "dimethyl hydrogen phosphite"),
    ],
)
def test_phosphite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COP(=O)(OC)OC", "trimethyl phosphate", id="phosphate_ester_unaffected"),
        pytest.param("c1ccccc1P(=O)(O)O", "phenylphosphonic acid", id="benzene_ring"),
    ],
)
def test_phosphate_ester_unaffected_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OP(=O)(O)CCP(=O)(O)O", id="second_phosphonic_acid_group"),
        pytest.param("NCP(=O)(O)O", id="unrecognized_heteroatom__phosphonic_acid"),
    ],
)
def test_rejects_second_phosphonic_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[PH4+]", "phosphanium", id="phosphanium"),
        pytest.param("C[P+](C)(C)C", "tetramethylphosphanium", id="tetramethylphosphanium"),
    ],
)
def test_phosphanium_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C[P+](C)(C)C(C)C", id="branched_quaternary_phosphonium_raises"),
        pytest.param("C[P+]1(C)CCCC1", id="ring_quaternary_phosphonium_raises"),
        pytest.param("C[P+](C)(C)Cl", id="halogen_substituted_quaternary_phosphonium_raises"),
    ],
)
def test_branched_quaternary_phosphonium_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_phosphonium_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P@H+](CC)CCC")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P@@+](CC)(CCC)CCCC")


def test_tetraphenylphosphanium():
    assert smiles_to_iupac("c1ccccc1[P+](c1ccccc1)(c1ccccc1)c1ccccc1") == "tetraphenylphosphanium"


def test_substituted_phenyl_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[P+](C)(C)C")
