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
        pytest.param("P1PPPP1", "pentaphospholane", id="cyclic_phosphane_chain"),
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
    "smiles,expected",
    [
        pytest.param("OP(=O)(O)CCP(=O)(O)O", "(ethane-1,2-diyl)bis(phosphonic acid)", id="carbon_linker"),
        pytest.param("N(CP(O)(O)=O)CP(O)(O)=O", "[azanediylbis(methylene)]bis(phosphonic acid)", id="hetero_linker_outranks_amine"),
        pytest.param(
            "P(CP(O)(O)=O)(CP(O)(O)=O)CP(O)(O)=O", "[phosphanetriyltris(methylene)]tris(phosphonic acid)", id="three_units"
        ),
    ],
)
def test_phosphonic_acids_joined_by_a_linker_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CP(C)O", "dimethylphosphinous acid"),
        ("c1ccccc1P(c1ccccc1)S", "diphenylphosphinothious acid"),
        ("CP(O)(=O)OC#N", "methylphosphonocyanatidic acid"),
        ("P(O)(=O)(OC#N)OC#N", "phosphorodicyanatidic acid"),
        ("c1ccccc1P(O)(=O)Cl", "phenylphosphonochloridic acid"),
        ("CP(O)(=S)Cl", "methylphosphonochloridothioic O-acid"),
        ("CP(O)Cl", "methylphosphonochloridous acid"),
        ("OP(=S)(O)O", "phosphorothioic O,O,O-acid"),
        ("OP(=O)(S)O", "phosphorothioic S-acid"),
        ("SP(=O)(S)O", "phosphorodithioic S,S-acid"),
        ("CP(=[Se])(O)O", "methylphosphonoselenoic O,O-acid"),
        ("CP(=[Se])(C)O", "dimethylphosphinoselenoic O-acid"),
        ("CP(=S)(C)S", "dimethylphosphinodithioic acid"),
        ("c1ccccc1[As](c1ccccc1)S", "diphenylarsinothious acid"),
        ("CB(C)S", "dimethylborinothioic acid"),
        ("CB(O)S", "methylboronothioic acid"),
    ],
)
def test_mononuclear_oxoacids_modified_by_functional_replacement(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[PH4][PH3][PH4]", "1λ5,2λ5,3λ5-triphosphane"),
        ("C[PH3]P", "1-methyl-1λ5-diphosphane"),
        ("[PH4]O[PH4]", "1λ5,3λ5-diphosphoxane"),
        ("POP", "diphosphoxane"),
        ("[SiH3]S[SiH3]", "disilathiane"),
        (
            "C1(=CC=CC=C1)[P+]1(P=P(CC1)(C1=CC=CC=C1)C1=CC=CC=C1)C1=CC=CC=C1",
            "1,1,3,3-tetraphenyl-4,5-dihydro-1H-1,2,3λ5-triphosphol-1-ium",
        ),
    ],
)
def test_heteroatom_chains_cite_nonstandard_bonding_numbers_and_alternate_by_seniority(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CPCl", "methylphosphinous chloride"),
        ("c1ccccc1P(c1ccccc1)Cl", "diphenylphosphinous chloride"),
        ("CP(Cl)Cl", "methylphosphonous dichloride"),
        ("C(C)=[N+](O)[O-]", "ethylideneazinic acid"),
        ("C[N+](C)(O)[O-]", "dimethylazinic acid"),
    ],
)
def test_halides_of_phosphorus_iii_acids_and_azinic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C(C)N=P(C1=CC=CC=C1)(C1=CC=CC=C1)C1=CC=CC=C1", "N-ethyl-P,P,P-triphenyl-λ5-phosphanimine", id="phosphine_imide_from_the_book"),
        pytest.param("C[P](C)(C)=NC", "N,P,P,P-tetramethyl-λ5-phosphanimine", id="all_methyl_phosphine_imide"),
        pytest.param("C[P](C)(C)=N", "trimethyl-λ5-phosphanimine", id="unsubstituted_imino_nitrogen"),
        pytest.param("C[As](C)(C)=N", "trimethyl-λ5-arsanimine", id="arsenic_analogue"),
    ],
)
def test_phosphane_imides_are_lambda5_phosphanimines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[P@](=O)(CCC)C1=CC=CC=C1", "(S)-methyl(phenyl)(propyl)-λ5-phosphanone"),
        ("C[P@@](OC)(=O)C1=CC=CC=C1", "methyl (S)-[methyl(phenyl)phosphinate]"),
    ],
)
def test_stereogenic_phosphorus_oxide_and_phosphinate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CN(C)P(=O)(N=C=S)S", "N,N-dimethylphosphoramid(isothiocyanatido)thioic S-acid"),
        ("ClP(=Nc1ccccc1)(S)c1ccccc1", "N,P-diphenylphosphonochloridimidothioic acid"),
        ("CNP(=S)(O)c1ccccc1", "N-methyl-P-phenylphosphonamidothioic O-acid"),
        ("CN=P(O)(c1ccccc1)c1ccccc1", "N-methyl-P,P-diphenylphosphinimidic acid"),
        ("NNP(=N)(O)O", "phosphorohydrazidimidic acid"),
        ("OP(=O)(O)OS", "phosphoro(thioperoxoic) OS-acid"),
        ("CCP(=[Se])(O)O", "ethylphosphonoselenoic O,O-acid"),
        ("CP(=O)(OC#N)O", "methylphosphonocyanatidic acid"),
        ("BrP(Cl)c1ccccc1", "phenylphosphonous bromide chloride"),
        ("CN(C)P(=O)(N=C=O)Cl", "N,N-dimethylphosphoramidisocyanatidic chloride"),
        ("N=P(N=C=S)(N=C=S)N=C=S", "phosphorimidic triisothiocyanate"),
        ("O=P(N=C=O)(N=C=O)N=C=O", "phosphoryl triisocyanate"),
        ("CNNP(=O)(NNC)NNC", "2,2',2''-trimethylphosphoric trihydrazide"),
        ("CN(C)P(C)(=O)N(C)C", "N,N,N',N',P-pentamethylphosphonic diamide"),
        ("CN(C)P(=O)(O)O", "dimethylphosphoramidic acid"),
        ("CCSP(=O)(OC)c1ccccc1", "S-ethyl O-methyl phenylphosphonothioate"),
        ("CS[As](C)(C)=O", "S-methyl dimethylarsinothioate"),
        ("CCSP(=S)(CC)CC", "ethyl diethylphosphinodithioate"),
        ("COP(Cl)N(C)C", "methyl N,N-dimethylphosphoramidochloridite"),
        ("CCOP(=O)(N=C=S)N(CC)CC", "ethyl N,N-diethylphosphoramid(isothiocyanatidate)"),
        ("CS[Sb](=O)(F)F", "S-methyl stiborodifluoridothioate"),
        ("CN=P(O)(O)OP(=N)(O)O", "N1-methyl-1,3-diimidodiphosphoric acid"),
    ],
)
def test_noncarbon_oxoacids_with_infixes_and_their_halides_amides_and_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[O-]P(=O)([O-])OP(=O)([O-])[O-].[Na+].[Na+].[Na+].[Na+]", "tetrasodium diphosphate"),
        ("[O-]S(=O)(=O)OS(=O)(=O)O.[Na+]", "sodium hydrogen disulfate"),
        ("NNP(=O)OP(=O)N", "3-hydrazidodiphosphonic 1-amide"),
        ("N#CP(I)P(I)I", "cyanohypodiphosphorous triiodide"),
        ("NP(O)OP(N)(=O)O", "{[amino(hydroxy)phosphanyl]oxy}phosphonamidic acid"),
        ("OP(=O)(O)NP(=O)O", "N-(hydroxyphosphonoyl)phosphoramidic acid"),
        ("NP(N)OP(N)(=O)N", "phosphorodiamidic phosphorodiamidous anhydride"),
        ("OS(=O)(=O)O[Se](=O)(=O)O", "selenic sulfuric monoanhydride"),
        ("CC(=O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O", "[({[(acetyloxy)(hydroxy)phosphoryl]oxy}(hydroxy)phosphoryl)oxy]phosphonic acid"),
        ("O[As](O)(=O)CCCCP(O)(O)=O", "(4-arsonobutyl)phosphonic acid"),
        ("O=P(O)(O)CCCP(=O)(O)CC", "{3-[ethyl(hydroxy)phosphoryl]propyl}phosphonic acid"),
        ("COS(=O)(=O)c1ccccc1P(=O)(O)O", "[2-(methoxysulfonyl)phenyl]phosphonic acid"),
        ("CP(C)P(O)P(C)C", "bis(dimethylphosphanyl)phosphinous acid"),
        ("CC(=O)OC(=O)OC(=O)O", "{[(acetyloxy)carbonyl]oxy}formic acid"),
    ],
)
def test_polynuclear_oxoacid_salts_derivatives_and_senior_acid_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_polyacid_substituent_is_flagged_as_not_preferred():
    from smiles_to_iupac import NonPreferredNameWarning

    with pytest.warns(NonPreferredNameWarning, match="P-67.2.6"):
        name = smiles_to_iupac("OC(=O)CCOP(=O)(OP(=O)(O)O)OP(=O)(O)O")
    assert name == "3-{[bis(phosphonooxy)phosphoryl]oxy}propanoic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("BrP(=O)(Br)CCP(=O)(Cl)Cl", "(2-phosphorodibromidoylethyl)phosphonic dichloride"),
        ("BrP(Br)(=O)CCP(=S)(Cl)Cl", "(2-phosphorodichloridothioylethyl)phosphonic dibromide"),
        ("ClP(Cl)(=O)OCCOP(=O)(Cl)N", "2-(phosphoramidochloridoyloxy)ethyl phosphorodichloridate"),
        ("OC(=O)CCP(=O)(N(C)C)N(C)C", "3-(tetramethylphosphorodiamidoyl)propanoic acid"),
        ("OC(=O)CCP(=S)(OC)OC", "3-(dimethoxyphosphorothioyl)propanoic acid"),
        ("OC(=O)CCP(=S)(S)S", "3-trithiophosphonopropanoic acid"),
    ],
)
def test_infix_acyl_prefixes_and_senior_acid_derivative_among_centres(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "O=S(=O)(c1ccccc1)N=P(Nc1ccccc1)(Nc1ccccc1)Nc1ccccc1",
            "N-(trianilino-λ5-phosphanylidene)benzenesulfonamide",
        ),
        ("CC(=O)N=P(C)(C)C", "N-(trimethyl-λ5-phosphanylidene)acetamide"),
        ("O=S(=O)(c1ccccc1)NP(C)C", "N-(dimethylphosphanyl)benzenesulfonamide"),
        ("O=S(=O)(c1ccccc1)N=CC", "N-ethylidenebenzenesulfonamide"),
    ],
)
def test_amides_carrying_phosphanylidene_phosphanyl_and_alkylidene_groups_on_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "OC(=O)CCOP(=O)(O)OP(=O)(O)O",
            "3-[(1,3,3-trihydroxy-1,3-dioxo-1λ5,3λ5-diphosphoxan-1-yl)oxy]propanoic acid",
        ),
        (
            "OC(=O)COS(=O)OS(=O)OS(=O)OC",
            "3,5,7-trioxo-2,4,6,8-tetraoxa-3λ4,5λ4,7λ4-trithiadecan-10-oic acid",
        ),
        (
            "OC(=O)CCSS(=S)SS(=S)SC",
            "3,5-bis(sulfanylidene)-2,3λ4,4,5λ4,6-pentathianonan-9-oic acid",
        ),
        (
            "OC(=O)CS(=O)(=O)OS(=O)(=O)O",
            "(3-hydroxy-1,1,3,3-tetraoxo-1λ6,3λ6-dithioxan-1-yl)acetic acid",
        ),
    ],
)
def test_chains_of_acid_centres_as_skeletal_prefixes_and_replacement_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
