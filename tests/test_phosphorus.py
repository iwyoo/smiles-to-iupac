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


def test_multiplied_compound_substituent_sorting_first():
    assert smiles_to_iupac("CCCP(C(C)C)C(C)C") == "di(propan-2-yl)(propyl)phosphane"


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
    "smiles, expected",
    [
        pytest.param("O=Pc1ccccc1C", "(2-methylphenyl)phosphanone", id="substituted_ring_substituent"),
        pytest.param("S=P(c1ccccc1)(c1ccccc1)c1ccccc1", "triphenyl-λ5-phosphanethione", id="thione"),
        pytest.param("C[As](C)(C)=O", "trimethyl-λ5-arsanone", id="arsane_oxide"),
        pytest.param("O=P(CCO)(CCO)CCO", "tris(2-hydroxyethyl)-λ5-phosphanone", id="substituted_alkyl_groups"),
    ],
)
def test_group_15_chalcogenide_with_any_organyl_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


def test_halogen_on_a_quaternary_phosphonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P+](C)(C)Cl")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[P+](C)(C)C(C)C", "trimethyl(propan-2-yl)phosphanium", id="branched_substituent"),
        pytest.param("Cc1ccccc1[P+](C)(C)C", "trimethyl(2-methylphenyl)phosphanium", id="substituted_phenyl"),
        pytest.param("OC(=O)C[PH3+]", "(carboxymethyl)phosphanium", id="acid_group_in_a_substituent"),
    ],
)
def test_phosphonium_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phosphonium_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P@H+](CC)CCC")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[P@@+](CC)(CCC)CCCC")


def test_tetraphenylphosphanium():
    assert smiles_to_iupac("c1ccccc1[P+](c1ccccc1)(c1ccccc1)c1ccccc1") == "tetraphenylphosphanium"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=P(O)(O)CCN", "(2-aminoethyl)phosphonic acid", id="phosphonic_acid_with_an_amino_substituent"),
        pytest.param("O=P(O)(CCN)CCN", "bis(2-aminoethyl)phosphinic acid", id="phosphinic_acid_with_amino_substituents"),
        pytest.param("O=P(O)(O)c1ccc(O)cc1", "(4-hydroxyphenyl)phosphonic acid", id="phosphonic_acid_on_a_substituted_ring"),
    ],
)
def test_phosphorus_acids_with_other_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCOP(=S)(OCC)Oc1ccc([N+](=O)[O-])cc1", "O,O-diethyl O-(4-nitrophenyl) phosphorothioate", id="thiono_ester"),
        pytest.param("CCOP(=O)(OCC)SCCN", "S-(2-aminoethyl) O,O-diethyl phosphorothioate", id="thiolo_ester"),
        pytest.param("CCOP(=S)(OCC)SCSCC", "O,O-diethyl S-[(ethylsulfanyl)methyl] phosphorodithioate", id="dithioate"),
        pytest.param("COP(C)(=S)OC", "O,O-dimethyl methylphosphonothioate", id="phosphonothioate"),
    ],
)
def test_phosphorus_thio_oxoacid_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OS(O)(=O)=O", "sulfuric acid"),
        ("OS(O)=O", "sulfurous acid"),
        ("O[N+]([O-])=O", "nitric acid"),
        ("OP(O)(O)=O", "phosphoric acid"),
        ("OP(O)O", "phosphorous acid"),
        ("OB(O)O", "boric acid"),
        ("OCl(=O)(=O)=O", "perchloric acid"),
        ("O[As](O)(O)=O", "arsoric acid"),
    ],
)
def test_free_mononuclear_oxoacids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=[As](O)(O)c1ccccc1", "phenylarsonic acid"),
        ("CC[Sb](=O)(O)CC", "diethylstibinic acid"),
        ("O=P(Cl)(Cl)Cl", "phosphoryl trichloride"),
        ("O=P(c1ccccc1)(Cl)Cl", "phenylphosphonic dichloride"),
        ("CCP(=S)(CC)Cl", "diethylphosphinothioic chloride"),
        ("CN(C)P(=O)(N(C)C)N(C)C", "hexamethylphosphoric triamide"),
        ("CN(C)P(=O)(C)C", "N,N,P,P-tetramethylphosphinic amide"),
        ("CNP(=O)(NC)c1ccccc1", "N,N'-dimethyl-P-phenylphosphonic diamide"),
        ("CN(C)P(=S)(N(C)C)c1ccccc1", "N,N,N',N'-tetramethyl-P-phenylphosphonothioic diamide"),
    ],
)
def test_group_15_acids_halides_and_amides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CN(C)P(=O)(Cl)Cl", "N,N-dimethylphosphoramidic dichloride"),
        ("CCP(=O)(N(C)C)Cl", "P-ethyl-N,N-dimethylphosphonamidic chloride"),
    ],
)
def test_amidic_halides_of_group_15_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCOP(=O)(OCC)S", "O,O-diethyl hydrogen phosphorothioate"),
        ("[Na+].[O-]P(=O)(O)O", "sodium dihydrogen phosphate"),
        ("[Ca+2].[O-]P(=O)(O)O.[O-]P(=O)(O)O", "calcium bis(dihydrogen phosphate)"),
        ("O=P(N=C=O)(N=C=O)N=C=O", "phosphoryl triisocyanate"),
        ("CP(=O)(C#N)C#N", "methylphosphonic dicyanide"),
        ("N[Pt](N)(Cl)Cl", "diazanidodichloridoplatinum"),
    ],
)
def test_partial_thioesters_oxoanion_salts_pseudohalides_and_amido_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CP(=O)=O", "methyl-λ5-phosphanedione", id="phosphanedione"),
        pytest.param("CC(C)CN[As]=O", "[(2-methylpropyl)amino]arsanone", id="amino_substituted_arsanone"),
        pytest.param("CCC=S=O", "propylidene-λ4-sulfanone", id="thioketone_oxide"),
    ],
)
def test_heterones_of_group_14_15_16_atoms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OP(=O)(O)CC(O)=O", "phosphonoacetic acid"),
        ("COP(=O)(OC)SCCC(O)=O", "3-[(dimethoxyphosphoryl)sulfanyl]propanoic acid"),
        ("CP(=O)(C)CC(O)=O", "(dimethylphosphinoyl)acetic acid"),
        ("OC(=O)c1ccc(cc1)P(=O)(O)c1ccc(cc1)C(O)=O", "4,4'-(hydroxyphosphoryl)dibenzoic acid"),
        ("OC(=O)c1ccc(cc1)P(=O)(C)c1ccc(cc1)C(O)=O", "4,4'-(methylphosphonoyl)dibenzoic acid"),
    ],
)
def test_phosphorus_acid_groups_cited_as_prefixes_under_a_carboxylic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
