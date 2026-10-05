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


def test_unsaturated_substituent():
    assert smiles_to_iupac("C=CP") == "ethenylphosphane"


def test_non_aromatic_ring():
    assert smiles_to_iupac("C1CCCCC1P") == "cyclohexylphosphane"


def test_halogenated_phenyl_mixed_with_alkyl():
    assert smiles_to_iupac("CP(c1ccc(Cl)cc1)") == "(4-chlorophenyl)(methyl)phosphane"


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


def test_cyclic_phosphane_chain():
    assert smiles_to_iupac("P1PPPP1") == "2,3,4,5-tetrahydro-1H-pentaphosphole"


def test_carbon_phosphorus_mix_is_a_diphosphane():
    assert smiles_to_iupac("CPP") == "methyldiphosphane"


def test_zero_substituents_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=P")


def test_single_phenyl_substituent():
    assert smiles_to_iupac("O=Pc1ccccc1") == "phenylphosphanone"


def test_mixed_alkyl_and_phenyl_substituents():
    assert smiles_to_iupac("CP(=O)c1ccccc1") == "methyl(phenyl)-λ5-phosphanone"


def test_substituted_phenyl_raises():
    # A substituted ring is not the Blue Book's plain 'phenyl' shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=Pc1ccccc1C")


def test_branched_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)P=O")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=P1CCCCC1")


def test_salt_of_partial_ester():
    assert smiles_to_iupac("COP(=O)(O)[O-].[Na+]") == "sodium methyl hydrogen phosphate"
    assert smiles_to_iupac("CCOP(=O)(O)[O-].[Na+]") == "sodium ethyl hydrogen phosphate"
    assert smiles_to_iupac("COP(=O)(O)[O-].[K+]") == "potassium methyl hydrogen phosphate"


def test_salt_of_partial_ester_multivalent_cation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COP(=O)(O)[O-].[Ca+2]")


def test_phosphoric_acid_itself_unaffected():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OP(=O)(O)O")


def test_phosphindole():
    assert smiles_to_iupac("C1C=C2C=CC=CC2=P1") == "2H-phosphindole"


def test_asymmetric_substituents():
    assert smiles_to_iupac("CCP(C)(=O)O") == "ethyl(methyl)phosphinic acid"
    assert smiles_to_iupac("CC(C)P(C)(=O)O") == "methyl(propan-2-yl)phosphinic acid"


def test_rejects_second_phosphinic_acid_group():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OP(C)(=O)CCP(C)(=O)O")


def test_rejects_unrecognized_heteroatom():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCP(C)(=O)O")


def test_phosphinine():
    assert smiles_to_iupac("C1=CC=PC=C1") == "phosphinine"


def test_methylphosphinine():
    assert smiles_to_iupac("CC1=CC=CC=P1") == "2-methylphosphinine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COP(OC)OC", "trimethyl phosphite"),
        ("COP(OC)O", "dimethyl hydrogen phosphite"),
    ],
)
def test_phosphite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phosphate_ester_unaffected():
    # A P=O bond routes to `_phosphate.py` instead, unchanged.
    assert smiles_to_iupac("COP(=O)(OC)OC") == "trimethyl phosphate"


def test_benzene_ring():
    assert smiles_to_iupac("c1ccccc1P(=O)(O)O") == "phenylphosphonic acid"


def test_rejects_second_phosphonic_acid_group():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OP(=O)(O)CCP(=O)(O)O")


def test_rejects_unrecognized_heteroatom__phosphonic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCP(=O)(O)O")


def test_phosphanium():
    assert smiles_to_iupac("[PH4+]") == "phosphanium"


def test_tetramethylphosphanium():
    assert smiles_to_iupac("C[P+](C)(C)C") == "tetramethylphosphanium"


def test_branched_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+](C)(C)C(C)C")


def test_ring_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+]1(C)CCCC1")


def test_halogen_substituted_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+](C)(C)Cl")


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
