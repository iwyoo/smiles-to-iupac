import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[O-]", "propan-2-olate"),
        ("CCCCC[O-]", "pentan-1-olate"),
    ],
)
def test_alkoxide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("FC[O-]", "fluoromethanolate", id="alkoxide_with_halogen_substituent"),
        pytest.param("C=CC[O-]", "prop-2-en-1-olate", id="alkoxide_with_unsaturation"),
        pytest.param("CC(C)(C)[O-]", "tert-butoxide", id="tert_butoxide"),
        pytest.param("CC(C)C[O-]", "2-methylpropan-1-olate", id="branched_alkoxide"),
    ],
)
def test_alkoxide_with_halogen_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenoxide():
    assert smiles_to_iupac("[O-]c1ccccc1") == "phenoxide"
    assert smiles_to_iupac("[O-]c1ccc(C)cc1") == "4-methylphenoxide"
    assert smiles_to_iupac("[O-]c1ccccc1Cl") == "2-chlorophenoxide"


def test_phenoxide_non_alkyl_ring_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(O)c1ccccc1[O-].[Cu+]")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-]c1ccc(CC=C)cc1", "4-(prop-2-en-1-yl)phenoxide", id="phenoxide_unsaturated_ring_substituent_is_named"),
        pytest.param("[O-]CC[O-]", "ethane-1,2-bis(olate)", id="two_alkoxide_groups_is_named"),
        pytest.param("[O-]CCOC", "2-methoxyethan-1-olate", id="ether_oxygen_alongside_alkoxide_is_named"),
        pytest.param("CC[C@@H](Cl)[O-]", "(1R)-1-chloropropan-1-olate", id="alkoxide_stereocenter_with_coexisting_halogen"),
    ],
)
def test_phenoxide_unsaturated_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_alkoxide():
    assert smiles_to_iupac("c1ccccc1CC[O-]") == "2-phenylethan-1-olate"
    assert smiles_to_iupac("c1ccccc1CCC[O-]") == "3-phenylpropan-1-olate"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Cc1ccccc1CC[O-]", "2-(2-methylphenyl)ethan-1-olate", id="substituted_benzene_ring_alkoxide_is_named"),
        pytest.param("C=Cc1ccccc1CC[O-]", "2-(2-ethenylphenyl)ethan-1-olate", id="chain_alkoxide_unsaturation_is_named"),
    ],
)
def test_phenyl_substituted_benzene_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[NH-]", "propan-2-aminide"),
    ],
)
def test_aminide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[NH-]c1ccccc1", "benzenaminide", id="aromatic_aminide_is_named"),
        pytest.param("[NH-]CC[NH-]", "ethane-1,2-bis(aminide)", id="two_aminide_groups_is_named"),
        pytest.param("C=C[NH-]", "eth-1-en-1-aminide", id="enamine_aminide_is_named"),
        pytest.param("NCC[NH-]", "2-aminoethan-1-aminide", id="second_nitrogen_is_named"),
        pytest.param("Clc1ccc(CCC[NH-])cc1", "3-(4-chlorophenyl)propan-1-aminide", id="phenyl_chain_aminide_ring_halogen"),
        pytest.param("C=Cc1ccccc1CC[NH-]", "2-(2-ethenylphenyl)ethan-1-aminide", id="phenyl_chain_aminide_unsaturation_is_named"),
    ],
)
def test_aromatic_aminide_is_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[NH2+]C1CCCCC1", "N-methylcyclohexanaminium", id="ring_ammonium"),
        pytest.param("C[N+](C)(C)CC(O)CC(=O)O", "3-carboxy-2-hydroxy-N,N,N-trimethylpropan-1-aminium", id="cation_outranks_acid"),
        pytest.param("C[N+](C)(C)CCC(=O)OC", "2-(methoxycarbonyl)-N,N,N-trimethylethan-1-aminium", id="cation_outranks_ester"),
        pytest.param("C[N+](C)(C)c1ccc(C(=O)O)cc1", "4-carboxy-N,N,N-trimethylanilinium", id="aryl_ammonium_with_acid"),
    ],
)
def test_ammonium_parent_with_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aromatic_ammonium():
    assert smiles_to_iupac("c1ccccc1[NH3+]") == "anilinium"


def test_doubly_charged_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH3++]")


K = "κ"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)[O][Cu][O]C(C)=O", "diacetatocopper"),
        ("O=C([O][Zn][O]C(=O)c1ccccc1)c1ccccc1", "dibenzoatozinc"),
        ("O=C[O][Zn][O]C=O", "diformatozinc"),
        ("O=[N+]([O-])[O][Cu][O][N+](=O)[O-]", "dinitratocopper"),
        ("N#C[S][Hg][S]C#N", f"bis(thiocyanato-{K}S)mercury"),
        ("S=C=[N][Hg][N]=C=S", f"bis(thiocyanato-{K}N)mercury"),
        ("CC[S][Hg][S]CC", "diethanethiolatomercury"),
        ("C[N](C)[Ti]([N](C)C)([N](C)C)[N](C)C", "tetrakis(dimethylazanido)titanium"),
        ("C[P](C)[Zr][P](C)C", "bis(dimethylphosphanido)zirconium"),
        ("CC1=CC(C)=[O]->[Cu]<-[O-]1", f"(pentane-2,4-dionato-{K}2O,O')cuprate(1-)"),
        ("O=C1[O-]->[Fe]<-[O-]C1=O", f"(oxalato-{K}2O,O')ferrate(2-)"),
        ("O=C1[O][Cu][O]1", f"(carbonato-{K}2O,O')copper"),
    ],
)
def test_anionic_ligand_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[c-]1ccccc1", "benzenide", id="benzenide_name"),
        pytest.param("[CH-]1CCCCC1", "cyclohexanide", id="carbanide_with_ring_is_named"),
        pytest.param("ClC[CH2-]", "2-chloroethan-1-ide", id="carbanide_with_halogen_is_named"),
        pytest.param("C=CC[CH2-]", "but-3-en-1-ide", id="carbanide_with_unsaturation_is_named"),
        pytest.param("C[C-]=O", "1-oxoethan-1-ide", id="acetyl_anion_with_oxo_on_anion_carbon"),
        pytest.param("CCC(=O)[CH-]C", "3-oxopentan-2-ide", id="oxo_substituent_elsewhere_on_chain"),
        pytest.param("O=CC(=O)[CH-]C", "3,4-dioxobutan-2-ide", id="carbanide_with_second_ketone_is_named"),
        pytest.param("O=C[CH-]C", "1-oxopropan-2-ide", id="carbanide_with_aldehyde_shaped_carbonyl_is_named"),
        pytest.param("[CH3+]", "methylium", id="methylium_name"),
        pytest.param("[CH2+]C", "ethylium", id="ethylium_name"),
        pytest.param("C1C[CH+]C1", "cyclobutylium", id="cyclobutylium_name"),
        pytest.param("C[CH+]CC", "butan-2-ylium", id="branch_point_carbenium_butan_2_ylium_name"),
    ],
)
def test_benzenide_name_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C1CCC2CC[CH+]CC2C1", id="fused_ring_carbenium_raises"),
    ],
)
def test_branch_point_carbenium_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C+]=O", "acetylium", id="acetylium_acylium_cation_name"),
        pytest.param("CC(C)C[C+]=O", "3-methylbutanoylium", id="branched_acylium_name"),
        pytest.param("O=[C+]C1CCCCC1", "cyclohexanecarbonylium", id="cyclohexanecarbonylium_acylium_name"),
        pytest.param("O=[C+]c1ccccc1", "benzoylium", id="benzoylium_acylium_name"),
        pytest.param("[CH-]1C=CC=C1", "cyclopenta-2,4-dien-1-ide", id="cyclopentadienide_name"),
        pytest.param("C1=CC=CC1", "cyclopenta-1,3-diene", id="cyclopentadiene_neutral_parent_unaffected"),
        pytest.param("OP(=O)(O)OP(=O)(O)O", "diphosphoric acid", id="diphosphoric_acid"),
        pytest.param("C=[N+]([H])[O-]", "methanimine N-oxide", id="nitrone_unsubstituted_nitrogen"),
        pytest.param("CCC#[N+][O-]", "propanenitrile oxide", id="nitrile_oxide_propane"),
        pytest.param("C#[N+][O-]", "formonitrile oxide", id="nitrile_oxide_fulminic_acid"),
    ],
)
def test_acetylium_acylium_cation_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC#[N+][N-]C", id="nitrile_imide_still_out_of_scope"),
        pytest.param("NCCN.Cl.Cl", id="two_halide_fragments_raises"),
        pytest.param("CCO.CCO", id="plain_mixture_raises"),
    ],
)
def test_nitrile_imide_still_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_unsupported_base_fragment():
    assert (
        smiles_to_iupac("c1ccccc1CC(=O)Nc1ccccc1C(=O)OCCCC.Cl")
        == "butyl 2-[(phenylacetyl)amino]benzoate;hydrochloride"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[NH4+].[K+].[O-]S(=O)(=O)[O-]", "azanium potassium sulfate"),
        # A +1 and a +2 cation together balancing phosphate's -3 charge.
        ("[NH4+].[Ca+2].[O-]P(=O)([O-])[O-]", "azanium calcium phosphate"),
    ],
)
def test_multi_cation_salt_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_same_cation_type_prefers_multiplying_prefix():
    assert smiles_to_iupac("[Na+].[Na+].[O-]C(=O)[O-]") == "disodium carbonate"


def test_multi_cation_mismatched_charge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Al+3].[K+].[O-]S(=O)(=O)[O-]")


def test_oxidanium():
    assert smiles_to_iupac("[OH3+]") == "oxidanium"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(C)[OH2+]", id="branched_substituent_not_supported"),
        pytest.param("[OH2+]C1CCCCC1", id="ring_substituent_not_supported"),
    ],
)
def test_branched_substituent_not_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_mixed_alkyl_and_phenyl_substituents():
    assert smiles_to_iupac("C[O+](C)c1ccccc1") == "dimethyl(phenyl)oxidanium"


def test_substituted_phenyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[OH2+]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Na+].[O-][N+](=O)[O-]", "sodium nitrate"),
    ],
)
def test_polyatomic_anion_salt_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_methyl_radical_name():
    assert smiles_to_iupac("[CH3]") == "methyl"


def test_branched_chain_radical():
    assert smiles_to_iupac("[CH2]C(C)C") == "2-methylpropyl"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[CH]CC", "butan-2-yl", id="butan_2_yl_name"),
        pytest.param("[C](C)(C)C", "tert-butyl", id="tert_butyl_name"),
        pytest.param("CC[C](C)C", "2-methylbutan-2-yl", id="2_methylbutan_2_yl_name"),
    ],
)
def test_branch_point_radical_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH](C(C)C)C", "3-methylbutan-2-yl", id="branch_point_radical_with_further_branching"),
        pytest.param("CC1CC[CH]C1", "3-methylcyclopentyl", id="substituted_ring_radical"),
    ],
)
def test_branch_point_radical_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_radical_centers_diyl_name():
    assert smiles_to_iupac("[CH2][CH2]") == "ethane-1,2-diyl"


def test_three_radical_centers_name():
    assert smiles_to_iupac("[CH2][CH][CH2]") == "propane-1,2,3-triyl"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH]C", "ethylidene", id="ethylidene_radical_name"),
        pytest.param("[C]1CCC1", "cyclobutylidene", id="cyclobutylidene_radical_name"),
    ],
)
def test_ethylidene_radical_name_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branch_point_divalent_radical_name():
    assert smiles_to_iupac("[C](C)C") == "propan-2-ylidene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[C](=O)C(C)CCC", "2-methylpentanoyl", id="branched_acyl_radical_name"),
        pytest.param("[C](=O)C1CCCCC1", "cyclohexanecarbonyl", id="cyclohexanecarbonyl_acyl_radical_name"),
        pytest.param("[C](=O)c1ccccc1", "benzoyl", id="benzoyl_acyl_radical_name"),
        pytest.param("C[NH]", "methanaminyl", id="aminyl_radical_name"),
        pytest.param("CC=[N]", "ethaniminyl", id="iminyl_radical_name"),
        pytest.param("CC(=O)[NH]", "acetamidyl", id="amidyl_radical_name"),
        pytest.param("C#C[CH]", "prop-2-yn-1-ylidene", id="vinyl_carbyne_name"),
        pytest.param("CCCC[O]", "butoxyl", id="butoxyl_radical_name"),
        pytest.param("C[Se]", "methylselanyl", id="methylselanyl_radical_name"),
    ],
)
def test_branched_acyl_radical_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[CH]1CC[CH]CC1", "cyclohexane-1,4-diyl"),
    ],
)
def test_ring_diradical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diradical_branch_point_name():
    assert smiles_to_iupac("[CH2]C([CH2])C") == "2-methylpropane-1,3-diyl"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC[OH+]", "ethyloxidaniumyl", id="ethyloxidaniumyl"),
        pytest.param("C[N+]", "methanaminyliumyl", id="methanaminyliumyl"),
        pytest.param("CC=[N+]", "ethaniminyliumyl", id="ethaniminyliumyl"),
        pytest.param("CC(=O)[N+]", "acetamidyliumyl", id="ethanamidyliumyl"),
        pytest.param("[Ca+2].CCCCC(C)(C(=O)[O-])c1ccccc1.CCCCC(C)(C(=O)[O-])c1ccccc1", "calcium bis(2-methyl-2-phenylhexanoate)", id="calcium_bis_compound_carboxylate"),
        pytest.param("[Na+].C[O-]", "sodium methoxide", id="sodium_methoxide"),
        pytest.param("[NH4+].CC(=O)[S-]", "azanium ethanethioate", id="ammonium_ethanethioate"),
        pytest.param("[Ca+2].[Cl-].[Cl-]", "calcium dichloride", id="calcium_dichloride"),
    ],
)
def test_ethyloxidaniumyl_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_anions_cited_alphanumerically():
    assert smiles_to_iupac("[Ca+2].[Cl-].[Br-]") == "calcium bromide chloride"


def test_carbanide_salt_names():
    assert smiles_to_iupac("[CH3-].[Li+]") == "lithium methanide"
    assert smiles_to_iupac("C[CH2-].[Na+]") == "sodium ethanide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Valine zwitterion (branched chain).
        ("CC(C)C(C(=O)[O-])[NH3+]", "2-azaniumyl-3-methylbutanoate"),
        ("C[N+](C)(C)CCS(=O)(=O)[O-]", "2-(N,N-dimethylmethanaminiumyl)ethane-1-sulfonate"),
    ],
)
def test_zwitterion_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C1CC([NH3+])C1C(=O)[O-]", id="ionic_center_in_ring_raises"),
        pytest.param("[NH3+]C(S(=O)(=O)[O-])", id="ammonium_bonded_directly_to_sulfonate_carbon_raises"),
    ],
)
def test_zwitterion_ionic_center_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[O-]C(=O)CCC(=O)[O-]", "butanedioate"),
        ("[O-]C(=O)CCCCCC(=O)O", "6-carboxyhexanoate"),
        ("[O-]C(=O)C(=O)O", "hydrogen oxalate"),
        ("[O-]S(=O)(=O)c1ccccc1", "benzenesulfonate"),
        ("OC([O-])=O.[Na+]", "sodium hydrogen carbonate"),
        ("CC(=O)[O-].CC(=O)[O-].[Ca+2]", "calcium diacetate"),
    ],
)
def test_acid_anion_and_salt_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


_P72_ANIONS = [
    ("[CH3-]", "methanide"),
    ("N#C[C-](C#N)C#N", "tricyanomethanide"),
    ("[C-2](c1ccccc1)c1ccccc1", "diphenylmethanediide"),
    ("C[P-]C", "dimethylphosphanide"),
    ("C#[Si-]", "methylidynesilanide"),
    ("[c-]1ccccc1", "benzenide"),
    ("[CH-]1C=CC=C1", "cyclopenta-2,4-dien-1-ide"),
    ("[C-]#[C-]", "ethynediide"),
    ("C[CH-]C", "propan-2-ide"),
    ("[CH2-][CH2-]", "ethane-1,2-diide"),
    ("[SiH3-]", "silanide"),
    ("[NH2-]", "azanide"),
    ("[SH-]", "sulfanide"),
    ("[PH2-]", "phosphanide"),
    ("[NH-2]", "azanediide"),
    ("S[S-]", "disulfanide"),
    ("OO[O-]", "trioxidanide"),
    ("[N-]1C=CC=CC1", "pyridin-1(2H)-ide"),
    ("[CH-]1C=C[CH-]c2ccccc12", "1,4-dihydronaphthalene-1,4-diide"),
    ("CCC(=O)O[O-]", "propaneperoxoate"),
    ("CCC(=O)[S-]", "propanethioate"),
    ("[O-]S(=O)(=O)c1ccccc1", "benzenesulfonate"),
    ("[O-]P(Cc1ccccc1)Cc1ccccc1", "dibenzylphosphinite"),
    ("[O-]C(=O)c1cccc(n1)C(=O)[O-]", "pyridine-2,6-dicarboxylate"),
    ("[O-]C(=N)c1ccc[nH]1", "1H-pyrrole-2-carboximidate"),
    ("CC(=S)O[O-]", "ethaneperoxothioate"),
    ("CC(=O)O[S-]", "ethane(OS-thioperoxoate)"),
    ("OC(=O)CCCCC(=O)[O-]", "5-carboxypentanoate"),
    ("OP(=O)([O-])c1ccccc1", "hydrogen phenylphosphonate"),
    ("CCOC(=O)CCC(=O)[O-]", "4-ethoxy-4-oxobutanoate"),
    ("OP(=O)([O-])Oc1ccccc1", "phenyl hydrogen phosphate"),
    ("C[O-]", "methoxide"),
    ("CC[O-]", "ethoxide"),
    ("CC(C)(C)[O-]", "tert-butoxide"),
    ("[O-]c1ccccc1", "phenoxide"),
    ("CC(C)[O-]", "propan-2-olate"),
    ("[O-]c1ccccc1[O-]", "benzene-1,2-bis(olate)"),
    ("[S-]c1ccccc1[S-]", "benzene-1,2-bis(thiolate)"),
    ("CN(C)[O-]", "dimethylaminoxide"),
    ("CO[O-]", "methaneperoxolate"),
    ("CCS[O-]", "ethane(SO-thioperoxolate)"),
    ("[O-]OCCO[O-]", "ethane-1,2-bis(peroxolate)"),
    ("[S-]Sc1ccc(S[S-])cc1", "benzene-1,4-bis(dithioperoxolate)"),
    ("N[O-]", "aminoxide"),
    ("[OH-]", "hydroxide"),
    ("O[O-]", "hydroperoxide"),
    ("C[NH-]", "methanaminide"),
    ("[NH-]c1ccccc1", "benzenaminide"),
    ("[NH-]CC[NH-]", "ethane-1,2-bis(aminide)"),
    ("CP(C)(C)=[N-]", "trimethyl-λ5-phosphaniminide"),
    ("CC[N-2]", "ethanaminediide"),
    ("[N-2]c1ccccc1", "benzenaminediide"),
    ("CC(C)(C)OO[O-]", "tert-butyltrioxidanide"),
    ("C[C-]=O", "1-oxoethan-1-ide"),
    ("O[NH-]", "hydroxyazanide"),
    ("O[N-2]", "hydroxyazanediide"),
    ("C[SiH4-]", "methylsilanuide"),
    ("C[B-](C)(C)C", "tetramethylboranuide"),
    ("C[P-](C)(C)C", "tetramethylphosphanuide"),
    ("[S-](F)(F)c1ccccc1", "difluoro(phenyl)sulfanuide"),
    ("c1ccccc1[I-]c1ccccc1", "diphenyliodanuide"),
    ("C[B-]1(C)CCCCC1", "1,1-dimethylborinan-1-uide"),
    ("C[B-]1(C)CCC2(C1)CCCCC2", "2,2-dimethyl-2-boraspiro[4.5]decan-2-uide"),
    ("C1C[PH-]2CCC1CC2", "1-phosphabicyclo[2.2.2]octan-1-uide"),
    ("[PH-]c1ccc([PH-])cc1", "(1,4-phenylene)bis(phosphanide)"),
    ("O=C([NH-])CCC(=O)[NH-]", "butanedioylbis(azanide)"),
    ("[CH-]1CCC(CC1)S(=O)(=O)[O-]", "cyclohexan-1-ide-4-sulfonate"),
    ("[O-]C(=O)CCC#[C-]", "pent-1-yn-1-id-5-oate"),
    ("[O-]C(=O)Cc1ccccc1C(=O)[O-]", "2-(carboxylatomethyl)benzoate"),
    ("[O-]c1cc2ccccc2cc1C(=O)[O-]", "3-oxidonaphthalene-2-carboxylate"),
    ("[O-]S(=O)(=O)c1ccc(cc1)C(=O)[O-]", "4-sulfonatobenzoate"),
    ("OC(=O)c1ccc(cc1)C(=O)[O-]", "4-carboxybenzoate"),
    ("[O-]c1ccc(cc1)C(=O)[O-]", "4-oxidobenzoate"),
    ("C#[C-]", "ethynide"),
    ("C=[CH-]", "ethenide"),
    ("ClC[CH2-]", "2-chloroethan-1-ide"),
    ("[c-]1ccc(cc1)C(=O)[O-]", "benzen-1-ide-4-carboxylate"),
    ("[O-]C(=O)CC[CH2-]", "butan-1-id-4-oate"),
    ("CCCC=[N-]", "butan-1-iminide"),
    ("CCC(=[N-])C", "butan-2-iminide"),
    ("[S-]Oc1ccc(O[S-])cc1", "benzene-1,4-bis(OS-thioperoxolate)"),
    ("CCCC(=N)[O-]", "butanimidate"),
    ("[S-]C(=O)c1ccccc1", "benzenecarbothioate"),
    ("CCC(=S)[S-]", "propane(dithioate)"),
    ("N[N-2]", "hydrazine-1,1-diide"),
    ("N[NH-]", "hydrazin-1-ide"),
    ("C[As](C)(C)=[N-]", "trimethyl-λ5-arsaniminide"),
    ("[PH-]CC[AsH3-]", "(2-phosphanidylethyl)arsanuide"),
    ("[SiH2-]CC[PH-]", "(2-silanidylethyl)phosphanide"),
    ("[O-]CC([O-])C1CC[BH2-]CC1", "1-(borinan-1-uid-4-yl)ethane-1,2-bis(olate)"),
    ("[O-]CC([O-])[NH-]", "1-azanidylethane-1,2-bis(olate)"),
    ("[O-]CC([O-])[BH3-]", "1-boranuidylethane-1,2-bis(olate)"),
    ("[O-]P(=O)([O-])c1ccccc1", "phenylphosphonate"),
    ("[O-]P(=O)(O)O", "dihydrogen phosphate"),
    ("COP(=O)(OC)[O-]", "dimethyl phosphate"),
    ("COS(=O)(=O)[O-]", "methyl sulfate"),
    ("N#C[C-](C#N)c1ccc(cc1)[C-](C#N)C#N", "(1,4-phenylene)bis(dicyanomethanide)"),
    ("C[B-]1(C)CC[P-]CC1", "4,4-dimethyl-1,4-phosphaborinan-1-id-4-uide"),
    ("[CH-]1C=CC=CN1", "pyridin-2(1H)-ide"),
    ("CN1C=CC=C[C-]1C", "1,2-dimethylpyridin-2(1H)-ide"),
    ("[Na+].C[BH-](C)C", "sodium trimethylboranuide"),
    ("[Na+].[Na+].[O-]c1ccccc1[O-]", "disodium benzene-1,2-bis(olate)"),
    ("[Li+].CC(C)(C)[AlH-](CC(C)C)CC(C)C", "lithium tert-butylbis(2-methylpropyl)alumanuide"),
    ("[Na+].OP(=O)([O-])c1ccccc1", "sodium hydrogen phenylphosphonate"),
    ("F[I-](F)(F)(F)(F)F", "hexafluoro-λ5-iodanuide"),
]


@pytest.mark.parametrize(("smiles", "expected"), _P72_ANIONS)
def test_p72_anion_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("NC(=O)CC(=O)[O-]", "3-amino-3-oxopropanoate"),
        ("CNC(=O)CCC(=O)[O-]", "4-(methylamino)-4-oxobutanoate"),
        ("[O-]C(=O)c1ccc(cc1)[CH2-]", "4-methanidylbenzoate"),
        ("[O-]C(=O)c1ccc(cc1)[C-](C#N)C#N", "4-(dicyanomethanidyl)benzoate"),
        ("[O-]C(=O)C[C-]1C=CC=C1", "(cyclopenta-2,4-dien-1-id-1-yl)acetate"),
        ("[O-]C(=O)c1ccc(cc1)[BH3-]", "4-boranuidylbenzoate"),
        ("CC(=O)[Te-]", "ethanetelluroate"),
        ("[S-]C(=S)c1ccccc1", "benzenecarbodithioate"),
    ],
)
def test_p72_anion_prefixes_and_acid_variants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


_P72_COVERAGE = [
    ("F[Te-2](F)(F)(F)(F)(F)(F)F", "octafluoro-λ6-tellanediuide"),
    ("C1CCCC[SH3-]C[SH3-]1", "1λ6,3λ6-dithiocane-1,3-diide"),
    ("C1CC[SH3-]CC1", "1λ6-thian-1-ide"),
    ("[O-]P(=O)([O-])c1ccc(cc1)C(=O)[O-]", "4-phosphonatobenzoate"),
    ("[O-][As](=O)([O-])c1ccc(cc1)C(=O)[O-]", "4-arsonatobenzoate"),
    ("[O-]C(=O)c1ccc(cc1)CC[CH2-]", "4-(propan-1-id-3-yl)benzoate"),
    ("[O-]C(=O)c1ccc(cc1)C[CH-]C", "4-(propan-2-id-1-yl)benzoate"),
    ("[O-]C(=O)c1ccc(cc1)C=[CH-]", "4-(eth-1-en-1-id-2-yl)benzoate"),
    ("[O-]C(=O)c1ccc(cc1)S[S-]", "4-disulfanidylbenzoate"),
    ("[O-]C(=O)c1ccc(cc1)SS[S-]", "4-trisulfanidylbenzoate"),
    ("[S-]Sc1ccc(cc1)[S-]", "4-sulfidobenzene-1-dithioperoxolate"),
    ("[O-]C(=O)c1ccc(cc1)C2=C(c3ccc(cc3)C(=O)[O-])[CH-]C=C2", "4,4'-(cyclopenta-2,4-dien-1-id-2,3-diyl)dibenzoate"),
    ("[O-]C(=O)c1ccc(cc1)CCCc1ccc(cc1)C(=O)[O-]", "4,4'-(propane-1,3-diyl)dibenzoate"),
    ("N#C[C-](C#N)C1=C(C(=C(C#N)C#N)1)[C-](C#N)C#N", "[3-(dicyanomethylidene)cycloprop-1-ene-1,2-diyl]bis(dicyanomethanide)"),
    ("OC(=O)c1ccc(cc1)C1=C(c2ccc(cc2)C(=O)O)C1=C(C#N)C#N", "4,4'-[3-(dicyanomethylidene)cycloprop-1-ene-1,2-diyl]dibenzoic acid"),
    ("B1C=Cc2ccccc12", "1H-1-benzoborole"),
    ("C1=Cc2ccccc2[SiH2]1", "1H-1-benzosilole"),
    ("C1=Cc2ccccc2[PH]1", "1H-phosphindole"),
    ("C1C=C2C=CC=CC2=P1", "2H-phosphindole"),
    ("CO[B-]1(C)C=C(C)c2ccccc12", "1-methoxy-1,3-dimethyl-1H-1-benzoborol-1-uide"),
    ("C[B-]1(C)C=C2C=C[CH-]C2=C1", "2,2-dimethyl-2,4-dihydrocyclopenta[c]borol-4-id-2-uide"),
    ("CN1[C-2]C=CC=Cc2ccccc12", "1-methyl-1-benzazocine-2,2(1H)-diide"),
    ("C[N-]c1ccccc1", "N-methylbenzenaminide"),
    ("C[N-]C", "N-methylmethanaminide"),
    ("C[N-]CC[N-]C", "N1,N2-dimethylethane-1,2-bis(aminide)"),
    ("C[N-]CC(=O)[O-]", "(methylazanidyl)acetate"),
    ("CC(=O)[O-]", "acetate"),
    ("ClCC(=O)[O-]", "chloroacetate"),
    ("ClC(O)C(=O)[O-]", "chloro(hydroxy)acetate"),
    ("[O-]C(=O)Cc1ccccc1", "phenylacetate"),
    ("CC(=O)[NH-]", "acetylazanide"),
    ("CC(=O)N[N-2]", "acetylhydrazine-1,1-diide"),
    ("[Ca+2].CC(=O)[O-].CC(=O)[O-]", "calcium diacetate"),
    ("NC(=O)CC(=O)[O-]", "3-amino-3-oxopropanoate"),
]


@pytest.mark.parametrize(("smiles", "expected"), _P72_COVERAGE)
def test_p72_prefixes_chains_fusion_and_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2]CCl", "2-chloroethyl", id="halogen_on_the_radical_chain"),
        pytest.param("[CH2]CC(=O)O", "2-carboxyethyl", id="acid_cited_as_prefix_of_a_radical"),
        pytest.param("C[CH]C(=O)O", "1-carboxyethyl", id="radical_at_the_alpha_carbon"),
        pytest.param("[CH2]CO", "2-hydroxyethyl", id="hydroxy_beside_radical"),
        pytest.param("[CH2]CN", "2-aminoethyl", id="amino_beside_radical"),
        pytest.param("CC(=O)[CH2]", "2-oxopropyl", id="ketone_beside_radical"),
        pytest.param("[CH2]C#N", "cyanomethyl", id="nitrile_beside_radical"),
        pytest.param("N[CH]C(=O)O", "amino(carboxy)methyl", id="two_groups_on_a_one_carbon_radical"),
        pytest.param("C[C](C)C(=O)O", "2-carboxypropan-2-yl", id="tertiary_radical_with_acid"),
        pytest.param("OC1CC[CH]CC1", "4-hydroxycyclohexyl", id="ring_radical_with_hydroxy"),
        pytest.param("[CH2]C(=O)OC", "(methoxycarbonyl)methyl", id="ester_beside_radical"),
        pytest.param("[CH2]c1ccc(cc1)[N+](=O)[O-]", "(4-nitrophenyl)methyl", id="nitro_on_an_aryl_methyl_radical"),
        pytest.param("[O]CC(=O)O", "carboxymethoxyl", id="oxygen_radical_with_acid"),
        pytest.param("[O]CCCl", "2-chloroethoxyl", id="oxygen_radical_with_halogen"),
        pytest.param("CC(C)[O]", "(propan-2-yl)oxyl", id="branched_oxygen_radical"),
    ],
)
def test_radical_beside_characteristic_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[n+]1ccccc1", "1-methylpyridin-1-ium", id="n_alkylpyridinium"),
        pytest.param("c1cc[nH+]cc1", "pyridin-1-ium", id="protonated_pyridine"),
        pytest.param("c1cc[n+](cc1)c1ccccc1", "1-phenylpyridin-1-ium", id="n_arylpyridinium"),
        pytest.param("Cc1cc[n+](C)c(C)c1", "1,2,4-trimethylpyridin-1-ium", id="substituents_numbered_from_the_cationic_nitrogen"),
        pytest.param("C[n+]1ccccc1C(=O)O", "2-carboxy-1-methylpyridin-1-ium", id="cation_outranks_acid"),
        pytest.param("C[n+]1ccccc1[N+](=O)[O-]", "1-methyl-2-nitropyridin-1-ium", id="nitro_beside_the_cation"),
        pytest.param("C[n+]1cccc2ccccc12", "1-methylquinolin-1-ium", id="quinolinium"),
        pytest.param("CC[n+]1ccc2ccccc2c1", "2-ethylisoquinolin-2-ium", id="isoquinolinium"),
        pytest.param("c1ccc2[nH+]cccc2c1", "quinolin-1-ium", id="protonated_quinoline"),
        pytest.param("C[n+]1ccccc1.[Cl-]", "1-methylpyridin-1-ium chloride", id="pyridinium_salt"),
    ],
)
def test_heteroaromatic_ring_cation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-][n+]1ccccc1", "pyridine 1-oxide", id="pyridine_n_oxide"),
        pytest.param("Cc1cc[n+]([O-])cc1", "4-methylpyridine 1-oxide", id="substituted_pyridine_n_oxide"),
        pytest.param("OC(=O)c1cccc[n+]1[O-]", "2-carboxypyridine 1-oxide", id="n_oxide_outranks_acid"),
        pytest.param("[O-][n+]1cccc2ccccc12", "quinoline 1-oxide", id="quinoline_n_oxide"),
        pytest.param("[O-][n+]1ccc2ccccc2c1", "isoquinoline 2-oxide", id="isoquinoline_n_oxide_locant"),
    ],
)
def test_heteroaromatic_n_oxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1CC[NH2+]CC1", "piperidin-1-ium", id="protonated_piperidine"),
        pytest.param("C[NH+]1CCCCC1", "1-methylpiperidin-1-ium", id="n_alkylpiperidinium"),
        pytest.param("C[N+]1(C)CCCCC1", "1,1-dimethylpiperidin-1-ium", id="quaternary_ring_nitrogen"),
        pytest.param("C1COCC[NH2+]1", "morpholin-4-ium", id="cation_locant_from_the_ring_heteroatom_order"),
        pytest.param("C[NH+]1CCOCC1", "4-methylmorpholin-4-ium", id="n_alkylmorpholinium"),
        pytest.param("OC1CC[NH2+]CC1", "4-hydroxypiperidin-1-ium", id="hydroxy_on_a_cationic_ring"),
        pytest.param("C[N+]1(C)CCCC1C(=O)O", "2-carboxy-1,1-dimethylpyrrolidin-1-ium", id="cation_outranks_acid_on_a_saturated_ring"),
        pytest.param("c1ccc2oc[nH+]c2c1", "1,3-benzoxazol-3-ium", id="protonated_fused_azole"),
        pytest.param("[O-][N+]1(C)CCCCC1", "1-methylpiperidine 1-oxide", id="saturated_ring_n_oxide"),
        pytest.param("[O-][n+]1ccoc1", "1,3-oxazole 3-oxide", id="azole_n_oxide_locant"),
    ],
)
def test_ring_cation_and_n_oxide_beyond_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCn1cc[n+](C)c1", "1-ethyl-3-methyl-1H-imidazol-3-ium", id="dialkylimidazolium"),
        pytest.param("C[n+]1ccn(C)c1", "1,3-dimethyl-1H-imidazol-3-ium", id="symmetric_imidazolium"),
        pytest.param("C[n+]1ccsc1", "3-methyl-1,3-thiazol-3-ium", id="n_alkylthiazolium"),
        pytest.param("C[n+]1cccnc1C", "1,2-dimethylpyrimidin-1-ium", id="pyrimidinium_substituents"),
        pytest.param("C[n+]1ccnc2ccccc12", "1-methylquinoxalin-1-ium", id="fused_diazinium"),
    ],
)
def test_heteroaromatic_ring_cation_with_several_nitrogens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OCC[NH]", "2-hydroxyethan-1-aminyl", id="aminyl_with_hydroxy"),
        pytest.param("Clc1ccc([NH])cc1", "4-chlorobenzenaminyl", id="arylaminyl"),
        pytest.param("OCC(=O)[NH]", "2-hydroxyacetamidyl", id="amidyl_with_hydroxy"),
        pytest.param("OC(=O)CC(=O)[NH]", "2-carboxyacetamidyl", id="amidyl_outranks_acid"),
        pytest.param("CCCO[O]", "propylperoxyl", id="alkylperoxyl"),
        pytest.param("CC(=O)O[O]", "acetylperoxyl", id="acylperoxyl"),
    ],
)
def test_nitrogen_and_peroxyl_radicals_beside_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-][n+]1cccnc1", "pyrimidine 1-oxide", id="diazine_n_oxide"),
        pytest.param("Cc1cc[n+]([O-])cn1", "4-methylpyrimidine 1-oxide", id="n_oxide_locant_before_substituents"),
        pytest.param("Cc1ccnc[n+]1[O-]", "6-methylpyrimidine 1-oxide", id="oxide_nitrogen_numbered_lowest"),
        pytest.param("[O-][n+]1ccnc2ccccc12", "quinoxaline 1-oxide", id="fused_diazine_n_oxide"),
        pytest.param("Cc1cnc[nH+]c1", "5-methylpyrimidin-1-ium", id="protonated_diazine"),
        pytest.param("Cc1cnc[nH+]c1Cl", "6-chloro-5-methylpyrimidin-1-ium", id="protonated_diazine_numbered_from_the_cation"),
        pytest.param("c1ccc2nc[nH+]cc2c1", "quinazolin-3-ium", id="protonated_fused_diazine"),
    ],
)
def test_ring_cation_and_n_oxide_with_several_nitrogens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2]CC[CH2]", "butane-1,4-diyl", id="chain_diyl"),
        pytest.param("Cl[CH]C[CH2]", "1-chloropropane-1,3-diyl", id="diyl_with_halogen"),
        pytest.param("[CH2]C(Cl)C[CH2]", "2-chlorobutane-1,4-diyl", id="diyl_with_a_substituent_between"),
        pytest.param("OC(=O)C[CH][CH2]", "3-carboxypropane-1,2-diyl", id="adjacent_centres_with_acid"),
        pytest.param("[CH2]C(C(=O)O)C[CH2]", "2-carboxybutane-1,4-diyl", id="diyl_outranks_acid"),
        pytest.param("C1CC[CH][CH]C1", "cyclohexane-1,2-diyl", id="ring_diyl"),
    ],
)
def test_two_radical_centres_beside_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OCC[CH]", "3-hydroxypropylidene", id="ylidene_with_hydroxy"),
        pytest.param("OC(=O)C[CH]", "2-carboxyethylidene", id="ylidene_outranks_acid"),
        pytest.param("OC1CC[C]CC1", "4-hydroxycyclohexylidene", id="ring_ylidene_with_hydroxy"),
        pytest.param("O=C(O)C[C]", "2-carboxyethylidyne", id="ylidyne_with_acid"),
        pytest.param("C1CC[N]C1", "pyrrolidin-1-yl", id="ring_nitrogen_radical"),
        pytest.param("C1COCC[N]1", "morpholin-4-yl", id="ring_nitrogen_radical_locant_from_the_heteroatom_order"),
    ],
)
def test_ylidene_ylidyne_and_ring_nitrogen_radicals_beside_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2+]CCO", "3-hydroxypropan-1-ylium", id="primary_with_hydroxy"),
        pytest.param("C[CH+]C(=O)C", "3-oxobutan-2-ylium", id="secondary_with_oxo"),
        pytest.param("C[C+](C)CO", "1-hydroxy-2-methylpropan-2-ylium", id="tertiary_with_hydroxy"),
        pytest.param("[CH2+]C=C", "prop-2-en-1-ylium", id="allyl_cation"),
        pytest.param("[CH+]1CCCC1O", "2-hydroxycyclopentan-1-ylium", id="ring_cation_with_hydroxy"),
        pytest.param("OCC[CH+]CCN", "1-amino-5-hydroxypentan-3-ylium", id="two_groups"),
        pytest.param("C[C+](C)C", "2-methylpropan-2-ylium", id="three_branches"),
        pytest.param("[CH2+]C(C)C", "2-methylpropan-1-ylium", id="branched_chain"),
        pytest.param("CC1CC[CH+]C1", "3-methylcyclopentan-1-ylium", id="substituted_ring"),
        pytest.param("[CH2+]C(Cl)", "2-chloroethan-1-ylium", id="halogen"),
        pytest.param("[CH2+]C1CCC2(CC1)CCCC2", "(spiro[4.5]decan-8-yl)methylium", id="spiro_substituent"),
        pytest.param("[CH2+]C1CC2CCC1C2", "(bicyclo[2.2.1]heptan-2-yl)methylium", id="bridged_substituent"),
    ],
)
def test_carbenium_centre_beside_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2-][CH2]", "ethan-2-id-1-yl", id="radical_takes_the_lowest_locant"),
        pytest.param("[CH2+]C[CH2]", "propan-3-ium-1-yl", id="radical_cation"),
        pytest.param("[CH2-][CH][CH2-]", "propane-1,3-diid-2-yl", id="dianion_radical"),
        pytest.param("[CH-]1CCCC[CH]1", "cyclohexan-2-id-1-yl", id="ring_radical_anion"),
        pytest.param("[CH2+]CC[CH][CH2+]", "pentane-1,5-diium-2-yl", id="dication_radical"),
    ],
)
def test_radical_ions_on_a_hydrocarbon_skeleton(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCO[C]([O-])c1ccccc1", "ethoxy(oxido)(phenyl)methyl", id="oxido_prefix_on_the_radical_parent"),
        pytest.param("C[CH][N+](C)(C)C", "1-(trimethylazaniumyl)ethyl", id="cationic_prefix"),
        pytest.param("[CH2]C(=O)[O-]", "carboxylatomethyl", id="carboxylate_prefix"),
        pytest.param("[CH2]C[NH3+]", "2-azaniumylethyl", id="ammonium_prefix"),
        pytest.param("C[CH]C[S-]", "1-sulfidopropan-2-yl", id="sulfido_prefix"),
    ],
)
def test_ionic_groups_as_prefixes_of_a_radical_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
