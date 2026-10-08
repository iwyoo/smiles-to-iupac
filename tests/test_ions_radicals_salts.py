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
        pytest.param("CC#[O+]", "acetylium", id="acylium_oxygen_charged_resonance_form"),
        pytest.param("CCCC[C+]=S", "pentanethioylium", id="chalcogen_acylium_chain"),
        pytest.param("c1ccccc1C#[Se+]", "benzenecarboselenoylium", id="chalcogen_acylium_ring"),
        pytest.param("[CH-]1C=CC=C1", "cyclopenta-2,4-dien-1-ide", id="cyclopentadienide_name"),
        pytest.param("C1=CC=CC1", "cyclopenta-1,3-diene", id="cyclopentadiene_neutral_parent_unaffected"),
        pytest.param(
            "C1(=CCC=C1)[C@H](C1=C[CH-]C=C1)O",
            "3-[(S)-(cyclopenta-1,4-dien-1-yl)(hydroxy)methyl]cyclopenta-2,4-dien-1-ide",
            id="stereocentre_recognised_only_through_the_charged_ring",
        ),
        pytest.param("Cl[C@H](F)c1ccc(C(=O)O)cc1", "4-[(S)-chloro(fluoro)methyl]benzoic acid", id="methyl_substituent_descriptor_without_locant"),
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
        pytest.param("Cl.Cl", id="two_halide_fragments_raises"),
        pytest.param("CCO.CCO", id="plain_mixture_raises"),
    ],
)
def test_nitrile_imide_still_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


# P-74.2.1.1, P-74.2.2: ylides, imides, oxides and azoxy compounds
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[P+](C)(C)[C-](C)C", "2-(trimethylphosphaniumyl)propan-2-ide", id="phosphorus_ylide"),
        pytest.param("C[S+](C)[C-](CC)CC", "3-(dimethylsulfaniumyl)pentan-3-ide", id="sulfur_ylide"),
        pytest.param("CC(C)=[N+](C)[C-](C)C", "2-(N-methylpropan-2-iminiumyl)propan-2-ide", id="azomethine_ylide"),
        pytest.param("CC(C)=[O+][N-]C", "N-[(propan-2-ylidene)oxidaniumyl]methanaminide", id="carbonyl_imide"),
        pytest.param("CC(C)=[O+][C-](C)C", "2-[(propan-2-ylidene)oxidaniumyl]propan-2-ide", id="carbonyl_ylide"),
        pytest.param("CC#[N+][C-](C)C", "2-(acetonitriliumyl)propan-2-ide", id="nitrile_ylide"),
        pytest.param(
            "C[B-](C=C[P+](C1=CC=CC=C1)(C1=CC=CC=C1)C)(C)C",
            "trimethyl{2-[methyldi(phenyl)phosphaniumyl]ethen-1-yl}boranuide",
            id="onium_prefix_on_a_remote_boranuide",
        ),
        pytest.param(
            "C[Se+]1CC2=C(C([CH-]1)=O)C=CC=C2",
            "2-methyl-4-oxo-3,4-dihydro-1H-2-benzoselenopyran-2-ium-3-ide",
            id="ring_chalcogenium_with_ring_carbanide",
        ),
        pytest.param("[O-]c1cc[n+](C)cc1", "1-methylpyridin-1-ium-4-olate", id="ring_cation_with_an_anionic_group_suffix"),
        pytest.param(
            "c1ccccc1[N-]c1nn(-c2ccccc2)c[n+]1-c1ccccc1",
            "N,1,4-triphenyl-1H-1,2,4-triazol-4-ium-3-aminide",
            id="azolium_aminide_with_indicated_hydrogen",
        ),
        pytest.param("CC(C)=[O+][O-]", "2-(propan-2-ylidene)dioxidan-2-ium-1-ide", id="carbonyl_oxide"),
        pytest.param("C[N-][N+](C)=C", "1,2-dimethyl-2-methylidenehydrazin-2-ium-1-ide", id="azomethine_imide"),
        pytest.param("CC#[N+][N-]C", "2-ethylidyne-1-methylhydrazin-2-ium-1-ide", id="nitrile_imide"),
        pytest.param("CN=[N+](C)[N-]C", "1,2,3-trimethyltriaz-2-en-2-ium-1-ide", id="azo_imide"),
        pytest.param("c1ccccc1N=[N+]([O-])c1ccccc1", "diphenyldiazene oxide", id="symmetric_azoxy"),
        pytest.param("Clc1ccccc1N=[N+]([O-])c1ccccc1", "1-(2-chlorophenyl)-2-phenyldiazene 2-oxide", id="azoxy_with_the_oxide_locant"),
        pytest.param("CCC=[S+][O-]", "propylidene-\u03bb4-sulfanone", id="thioaldehyde_s_oxide"),
    ],
)
def test_dipolar_compounds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsupported_base_fragment():
    assert (
        smiles_to_iupac("c1ccccc1CC(=O)Nc1ccccc1C(=O)OCCCC.Cl")
        == "butyl 2-(2-phenylacetamido)benzoate;hydrochloride"
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


def test_pnictogen_cation_is_not_a_radical_centre():
    acid = "[O-]C(=O)CCC(=O)O"
    assert smiles_to_iupac(f"[Sb+3].{acid}.{acid}.{acid}") == "antimony tris(3-carboxypropanoate)"


def test_multi_cation_mismatched_charge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Al+3].[K+].[O-]S(=O)(=O)[O-]")


def test_oxidanium():
    assert smiles_to_iupac("[OH3+]") == "oxidanium"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)[OH2+]", "propan-2-yloxidanium", id="branched_substituent"),
        pytest.param("[OH2+]C1CCCCC1", "cyclohexyloxidanium", id="ring_substituent"),
        pytest.param("Cc1ccccc1[OH2+]", "(2-methylphenyl)oxidanium", id="substituted_phenyl"),
        pytest.param("CC(=O)C[O+](C)C", "dimethyl(2-oxopropyl)oxidanium", id="ketone_group_in_a_substituent"),
    ],
)
def test_oxonium_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_and_phenyl_substituents():
    assert smiles_to_iupac("C[O+](C)c1ccccc1") == "dimethyl(phenyl)oxidanium"


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
        pytest.param("CC=[N+]", "ethaniminylium", id="ethaniminylium"),
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
        pytest.param("C[S+](C)CC(=O)[O-]", "(dimethylsulfaniumyl)acetate", id="sulfonium_prefix_on_a_carboxylate"),
        pytest.param("C[P+](C)(C)CC(=O)[O-]", "(trimethylphosphaniumyl)acetate", id="phosphonium_prefix_on_a_carboxylate"),
        pytest.param("C[S+](C)CCS(=O)(=O)[O-]", "2-(dimethylsulfaniumyl)ethane-1-sulfonate", id="sulfonium_prefix_on_a_sulfonate"),
        pytest.param("C[S+](C)CCCCCC[N+](C)(C)C", "6-(dimethylsulfaniumyl)-N,N,N-trimethylhexan-1-aminium", id="senior_nitrogen_cation_is_the_parent"),
        pytest.param("C[S+](C)CCCCCC[P+](C)(C)C", "[6-(dimethylsulfaniumyl)hexyl]tri(methyl)phosphanium", id="phosphorus_outranks_sulfur"),
        pytest.param("C[O+](C)CCCC[N+](C)(C)C", "4-(dimethyloxidaniumyl)-N,N,N-trimethylbutan-1-aminium", id="oxygen_cation_prefix"),
    ],
)
def test_zwitterion_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("[NH3+]C(=O)[O-]", id="ammonium_bonded_directly_to_carboxylate_carbon_raises"),
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
        ("[Na+].[O-]C(=O)Cc1ccccc1C([O-])=O.[H+]", "sodium hydrogen 2-(carboxylatomethyl)benzoate"),
        ("[K+].[Na+].[O-]C(=O)CC(C(=O)[O-])CC(=O)[O-].[H+]", "potassium sodium hydrogen propane-1,2,3-tricarboxylate"),
        ("[Na+].[O-]P(=O)([O-])[O-].[H+].[H+]", "sodium dihydrogen phosphate"),
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


def test_ammonium_counter_ion_multiplied_with_bis_to_avoid_diazane_reading():
    assert smiles_to_iupac("[NH4+].[NH4+].[O-]S(=O)(=O)[O-]") == "bis(azanium) sulfate"


@pytest.mark.parametrize("smiles", ["CN.[O-]S(=O)(=O)O", "c1ccc(I(c2ccccc2)OS(=O)(=O)C)cc1"])
def test_unaccounted_fragment_or_polyvalent_halogen_is_never_dropped(smiles):
    with pytest.raises(NotImplementedError):
        smiles_to_iupac(smiles)

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
        pytest.param("[CH2]C(=O)OC", "2-methoxy-2-oxoethyl", id="ester_beside_radical"),
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
        pytest.param("[CH3+].[Cl-]", "methylium chloride", id="one_atom_cation_salt"),
        pytest.param(
            "O=C([O-])C[n+]1ccccc1CC(=O)[O-]", "2,2'-(pyridin-1-ium-1,2-diyl)diethanoate", id="ring_cation_linking_two_anions"
        ),
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
        pytest.param("[CH2+]CCC", "butylium", id="terminal_cation_of_unbranched_alkane_takes_the_alkyl_stem"),
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
        pytest.param("[CH2+]C[CH2]", "propan-3-ylium-1-yl", id="radical_cation"),
        pytest.param("[CH2-][CH][CH2-]", "propane-1,3-diid-2-yl", id="dianion_radical"),
        pytest.param("[CH-]1CCCC[CH]1", "cyclohexan-2-id-1-yl", id="ring_radical_anion"),
        pytest.param("[CH2+]CC[CH][CH2+]", "pentane-1,5-diylium-2-yl", id="dication_radical"),
        pytest.param("[CH2-]C(C)C[CH2]", "3-methylbutan-4-id-1-yl", id="branched_skeleton_radical_first"),
        pytest.param("ClC[CH-]C[CH2]", "4-chlorobutan-3-id-1-yl", id="substituted_skeleton"),
        pytest.param("[CH2-]CCC(C)(C)C[CH]CC", "4,4-dimethyloctan-1-id-6-yl", id="centre_set_before_radical_locant"),
        pytest.param("[CH2-]", "methanidyl", id="one_atom_radical_anion"),
        pytest.param("[CH2+]", "methyliumyl", id="one_atom_radical_cation"),
        pytest.param("C[CH+]", "ethan-1-ylium-1-yl", id="both_centres_on_one_atom"),
        pytest.param("[CH2-]C[CH-]", "propane-1,3-diid-1-yl", id="dianion_with_radical_on_an_ionic_atom"),
        pytest.param("c1ccccc1[CH-][CH][CH-]c1ccccc1", "1,3-diphenylpropane-1,3-diid-2-yl", id="aromatic_substituents_do_not_count_as_a_ring_parent"),
        pytest.param("c1ccc2c(c1)[CH-]c1ccccc1[CH]2", "9,10-dihydroanthracen-10-id-9-yl", id="named_ring_system_parent"),
    ],
)
def test_radical_ions_on_a_hydrocarbon_skeleton(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[NH-]", "azanidyl", id="heteroatom_hydride_radical_anion"),
        pytest.param("C[B-](C)C", "trimethylboranuidyl", id="boron_radical_anion"),
        pytest.param("C[N-]", "methanaminidyl", id="amine_derived_radical_anion"),
        pytest.param("CC(=O)[N-]", "acetylazanidyl", id="acyl_radical_anion"),
        pytest.param("c1ccccc1C#[N+]", "benzonitriliumyl", id="nitrilium_radical_cation"),
    ],
)
def test_radical_ions_named_through_the_filled_ion(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCO[C]([O-])c1ccccc1", "ethoxy(oxido)(phenyl)methyl", id="oxido_prefix_on_the_radical_parent"),
        pytest.param("C[CH][N+](C)(C)C", "1-(trimethylazaniumyl)ethyl", id="cationic_prefix"),
        pytest.param("[CH2]C(=O)[O-]", "carboxylatomethyl", id="carboxylate_prefix"),
        pytest.param("[CH2]C[NH3+]", "2-azaniumylethyl", id="ammonium_prefix"),
        pytest.param("C[CH]C[S-]", "1-sulfidopropan-2-yl", id="sulfido_prefix"),
        pytest.param("[CH2]c1cc[n+](C)cc1", "(1-methylpyridin-1-ium-4-yl)methyl", id="cationic_ring_prefix"),
        pytest.param("[CH2]c1ccc[nH+]c1", "(pyridin-1-ium-3-yl)methyl", id="protonated_ring_prefix"),
        pytest.param("[CH2]C[n+]1ccccc1", "2-(pyridin-1-ium-1-yl)ethyl", id="cationic_ring_attached_through_nitrogen"),
        pytest.param("c1ccc2c(c1)ccc[n+]2C[CH2]", "2-(quinolin-1-ium-1-yl)ethyl", id="fused_cationic_ring_through_nitrogen"),
        pytest.param("OO[CH2]", "hydroperoxymethyl", id="hydroperoxy_prefix"),
        pytest.param("CCOO[CH2]", "(ethylperoxy)methyl", id="alkylperoxy_prefix"),
    ],
)
def test_ionic_groups_as_prefixes_of_a_radical_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2][CH]C1CC[CH]CC1", "1-(4-ylocyclohexyl)ethane-1,2-diyl", id="most_centres_in_one_parent"),
        pytest.param("C[C]c1cccc(c1)C[CH2]", "2-[3-(1,1-diyloethyl)phenyl]ethyl", id="yl_outranks_ylidene_parent"),
        pytest.param("[CH2]C(C)(C)[O]", "(2-methyl-1-ylopropan-2-yl)oxyl", id="oxyl_parent_outranks_carbon"),
        pytest.param("[CH2]C1CC[CH]CC1", "4-(ylomethyl)cyclohexyl", id="ring_parent_outranks_chain"),
        pytest.param("[O]CC[N]C", "N-methyl-2-(ylooxidanyl)ethan-1-aminyl", id="nitrogen_outranks_oxygen"),
        pytest.param("[O]c1cccc(c1)C(=O)[O]", "[3-(ylooxidanyl)benzoyl]oxyl", id="acyl_oxyl_outranks_aryl_oxyl"),
        pytest.param("[O]CCC(=O)[O]", "[3-(ylooxidanyl)propanoyl]oxyl", id="acyl_oxyl_outranks_alkoxyl"),
        pytest.param("[O]C(=O)CC(=O)[O]", "[(oxylcarbonyl)acetyl]oxyl", id="radical_acyl_oxygen_prefix"),
        pytest.param("[CH2]CC[N]C(=O)C", "N-(3-ylopropyl)acetamidyl", id="amidyl_parent_outranks_carbon"),
        pytest.param("[O]CC[CH2]", "3-ylopropoxyl", id="oxygen_outranks_carbon"),
    ],
)
def test_choice_of_parent_radical_with_ylo_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[N+](C)=CC", "N,N-dimethylethaniminium", id="iminium_parent"),
        pytest.param("C[NH+]=C(C)CCl", "1-chloro-N-methylpropan-2-iminium", id="n_prefix_alphanumerical_order"),
        pytest.param("CC=[N+](C)C(C)C", "N-methyl-N-(propan-2-yl)ethaniminium", id="n_substituent_longer_than_chain"),
        pytest.param("C[N+](C)=CC(=O)OC", "1-(methoxycarbonyl)-N,N-dimethylmethaniminium", id="cation_outranks_ester"),
        pytest.param("C[N+](C)=CC(=O)[O-]", "(N-methylmethanaminiumylidene)acetate", id="iminium_ylidene_prefix"),
        pytest.param("C=[N+](C)CC(=O)[O-]", "(N-methylmethaniminiumyl)acetate", id="iminium_yl_prefix"),
        pytest.param("C1CC([NH3+])C1C(=O)[O-]", "2-azaniumylcyclobutane-1-carboxylate", id="ring_ammonium_prefix"),
        pytest.param("[CH2-]CC[CH2+]", "butan-4-ylium-1-ide", id="zwitterionic_skeleton"),
        pytest.param("[CH2]C=C[CH2+]", "but-2-en-4-ylium-1-yl", id="unsaturated_radical_cation"),
        pytest.param("[CH2-][CH+]C[CH2]", "butan-2-ylium-1-id-4-yl", id="three_kinds_of_centre"),
    ],
)
def test_iminium_and_mixed_ionic_centres(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1C[OH+]CC[OH+]1", "1,4-dioxane-1,4-diium", id="ring_heteroatoms_keep_the_e_before_diium"),
        pytest.param("C[N+]1(C)CC[N+](C)(C)CC1", "1,1,4,4-tetramethylpiperazine-1,4-diium", id="substituents_on_the_cationic_ring_nitrogens"),
        pytest.param("C[n+]1cc[n+](C)cc1", "1,4-dimethylpyrazine-1,4-diium", id="aromatic_ring_with_two_cationic_nitrogens"),
        pytest.param("C1C[NH+]2CC[NH+]1CC2", "1,4-diazabicyclo[2.2.2]octane-1,4-diium", id="bicyclic_skeleton_centres"),
        pytest.param("C[N+](C)=[N+](C)C", "tetramethyldiazene-1,2-diium", id="fully_substituted_pair_omits_locants"),
        pytest.param("C1CCC[SH3+]C[SH3+]C1", "1λ4,3λ4-dithiocane-1,3-diium", id="lambda4_ring_chalcogen_takes_ium"),
        pytest.param("C[NH2+][NH2+]C", "1,2-dimethylhydrazine-1,2-diium", id="partly_substituted_pair_cites_locants"),
        pytest.param("C[N+](C)(C)CC[N+](C)(C)C", "N1,N1,N1,N2,N2,N2-hexamethylethane-1,2-bis(aminium)", id="quaternary_bis_aminium"),
        pytest.param("C[NH2+]CC[NH2+]C", "N1,N2-dimethylethane-1,2-bis(aminium)", id="secondary_bis_aminium"),
        pytest.param("[NH3+]CC(C(=O)O)C[NH3+]", "2-carboxypropane-1,3-bis(aminium)", id="bis_aminium_with_a_prefix_group"),
        pytest.param("C1C[NH2+]CC[NH2+]1.[Cl-].[Cl-]", "piperazine-1,4-diium dichloride", id="dication_salt"),
    ],
)
def test_several_identical_cationic_centres(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[n+]1ccc(CCc2cc[n+](C)cc2)cc1", "4,4'-(ethane-1,2-diyl)bis(1-methylpyridin-1-ium)", id="ring_cations_joined_through_a_carbon"),
        pytest.param("c1cc[n+](cc1)CCCC[n+]1ccccc1", "1,1'-(butane-1,4-diyl)di(pyridin-1-ium)", id="ring_cations_joined_through_the_cationic_nitrogen"),
        pytest.param("C[n+]1ccc(CC(O)Cc2cc[n+](C)cc2)cc1", "4,4'-(2-hydroxypropane-1,3-diyl)bis(1-methylpyridin-1-ium)", id="substituted_linking_group"),
        pytest.param("[PH3+]c1ccc([PH3+])cc1", "(1,4-phenylene)bis(phosphanium)", id="mononuclear_cations_on_a_ring_linker"),
        pytest.param("C[C+](C)c1cccc(c1)[C+](C)C", "2,2'-(1,3-phenylene)di(propan-2-ylium)", id="carbenium_units_on_a_ring_linker"),
        pytest.param("C[n+]1ccc(CCc2cc[n+](C)cc2)cc1.[Br-].[Br-]", "4,4'-(ethane-1,2-diyl)bis(1-methylpyridin-1-ium) dibromide", id="assembly_as_the_cation_of_a_salt"),
    ],
)
def test_assemblies_of_parent_cations(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C1CC(C(=O)[NH2+]1)C1CC(=O)[NH2+]C1=O", "2,2',5,5'-tetraoxo[3,3'-bipyrrolidine]-1,1'-diium", id="imide_rings_joined_directly"),
        pytest.param("C[n+]1ccc(cc1)-c1cc[n+](C)cc1", "1,1'-dimethyl[4,4'-bipyridine]-1,1'-diium", id="substituted_aromatic_rings_joined_directly"),
        pytest.param("C1C[NH2+]CC1C1CC[NH2+]C1", "[3,3'-bipyrrolidine]-1,1'-diium", id="saturated_rings_joined_directly"),
        pytest.param("C[n+]1ccc(cc1)-c1cc[n+](C)cc1.[Cl-].[Cl-]", "1,1'-dimethyl[4,4'-bipyridine]-1,1'-diium dichloride", id="ring_assembly_as_the_cation_of_a_salt"),
    ],
)
def test_ring_assemblies_of_cations(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[n+]1cccc2[CH+]C=CC=Cc12", "1-methyl-5H-cyclohepta[b]pyridin-1-ium-5-ylium", id="fused_ring_with_ium_and_ylium"),
        pytest.param("C1C[CH+]CC[NH2+]1", "piperidin-1-ium-4-ylium", id="saturated_ring_ium_then_ylium"),
        pytest.param("C1CC[NH2+]C[CH+]1", "piperidin-1-ium-3-ylium", id="all_centres_lowest_then_ylium"),
        pytest.param("C[N+]1(C)CC[CH+]C[CH+]C1", "1,1-dimethylazepan-1-ium-3,5-bis(ylium)", id="several_ylium_centres"),
    ],
)
def test_ium_and_ylium_centres_in_one_ring_system(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[n+]1cc([N+](C)(C)C)cc2ccccc12", "N,N,N,1-tetramethylquinolin-1-ium-3-aminium", id="ring_ium_then_aminium"),
        pytest.param("C[n+]1cc2ccnc([N+](C)(C)C)c2cc1", "N,N,N,2-tetramethyl-2,6-naphthyridin-2-ium-5-aminium", id="skeletal_centre_takes_the_lowest_locant_before_the_suffix"),
        pytest.param("C[n+]1ccc([NH3+])cc1", "1-methylpyridin-1-ium-4-aminium", id="primary_aminium_on_the_ring"),
        pytest.param("C[n+]1ccc(cc1)C(=O)[NH3+]", "1-methylpyridin-1-ium-4-carboxamidium", id="ring_ium_then_carboxamidium"),
        pytest.param("C[n+]1ccc(cc1)C#[NH+]", "1-methylpyridin-1-ium-4-carbonitrilium", id="ring_ium_then_carbonitrilium"),
        pytest.param("C[n+]1ccc(N(C)C)cc1", "4-(dimethylamino)-1-methylpyridin-1-ium", id="neutral_substituted_amine_is_a_prefix_of_the_ring_cation"),
    ],
)
def test_skeletal_cationic_centre_with_a_cationic_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[N+]1(C)CCOCC1", "4,4-dimethylmorpholin-4-ium", id="quaternary_ring_nitrogen_beside_a_ring_oxygen"),
        pytest.param("C[N+]1(C)CCN(C)CC1", "1,1,4-trimethylpiperazin-1-ium", id="quaternary_ring_nitrogen_beside_a_neutral_ring_nitrogen"),
        pytest.param("C[N+]1(C)CCCOCC1.[Cl-]", "4,4-dimethyl-1,4-oxazepan-4-ium chloride", id="quaternary_ring_nitrogen_in_a_salt"),
        pytest.param("C[N+](C)(C)C1CCOCC1", "N,N,N-trimethyloxan-4-aminium", id="quaternary_ammonium_on_a_heterocycle"),
    ],
)
def test_quaternary_ammonium_on_or_in_a_heterocycle(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-][O-]", "dioxidanediide", id="peroxide_dianion"),
        pytest.param("[Na+].[O-][O-].[Na+]", "disodium dioxidanediide", id="peroxide_dianion_in_a_salt"),
    ],
)
def test_adjacent_anionic_oxygens_are_not_a_nitro_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CC[S+](C)C1", "1-methylthiolan-1-ium", id="ring_sulfonium"),
        pytest.param("C1CC[OH+]CC1", "oxan-1-ium", id="ring_oxonium"),
        pytest.param("C1CC[P+](C)(C)C1", "1,1-dimethylphospholan-1-ium", id="ring_phosphonium"),
        pytest.param("C[S+]1CCOCC1", "4-methyl-1,4-oxathian-4-ium", id="cationic_atom_numbered_among_two_ring_heteroatoms"),
        pytest.param("c1ccc2c(c1)C=C[S+]2C", "1-methyl-1-benzothiophen-1-ium", id="fused_ring_sulfonium"),
        pytest.param("C1CC[S+](C)C1.[I-]", "1-methylthiolan-1-ium iodide", id="ring_sulfonium_salt"),
    ],
)
def test_ring_cations_on_heteroatoms_other_than_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("c1cc[o+]cc1", "pyrylium", id="pyrylium"),
        pytest.param("Cc1cc(C)[o+]c(C)c1", "2,4,6-trimethylpyrylium", id="substituted_pyrylium"),
        pytest.param("c1cc[te+]cc1", "telluropyrylium", id="chalcogen_analogue_of_pyrylium"),
        pytest.param("c1ccc2[s+]c3ccccc3cc2c1", "thioxanthylium", id="xanthylium_family"),
        pytest.param("c1ccc2[o+]cccc2c1", "1λ4-benzopyran-1-ylium", id="lambda_name_of_a_benzopyrylium"),
        pytest.param("c1ccc2c[o+]ccc2c1", "2λ4-benzopyran-2-ylium", id="isobenzopyrylium"),
        pytest.param("c1ccc(cc1)-c1cc2ccccc2[o+]c1-c1ccccc1", "2,3-diphenyl-1λ4-benzopyran-1-ylium", id="substituted_lambda_name"),
        pytest.param("C1C=CC=[O+]1", "2H-1λ4-furan-1-ylium", id="five_membered_ring_with_indicated_hydrogen"),
        pytest.param("[S+]1=CCC=C1", "3H-1λ4-thiophen-1-ylium", id="indicated_hydrogen_on_the_third_atom"),
        pytest.param("c1cc[n+]2ccccc2c1", "5λ5-quinolizin-5-ylium", id="bridgehead_nitrogen_cation"),
        pytest.param("c1cc[s+]cc1.[Cl-]", "thiopyrylium chloride", id="ylium_ring_in_a_salt"),
        pytest.param("C1C=C[N+]2=CC=CC=C12", "1H-4λ5-indolizin-4-ylium", id="general_fused_parent_with_indicated_hydrogen"),
        pytest.param("C1=CC2=CC=CC=C2[I+]1", "1λ3-benziodol-1-ylium", id="halogen_ring_centre"),
        pytest.param("c1ccc2c(c1)[S+]=CCSC=C2", "3H-1λ4,4-benzodithiocin-1-ylium", id="lambda_joins_the_cited_heteroatom_locant"),
        pytest.param("c1cc[n+]2cc[n+]3ccccc3c2c1", "5λ5,8λ5-dipyrido[1,2-a:2',1'-c]pyrazine-5,8-diylium", id="two_cationic_centres"),
        pytest.param("C1=C2C(=CC=C1)[N-]C1=CC=3C=CC=C[N+]3C=C12", "5H-11λ5-indolo[2,3-b]quinolizin-11-ylium-5-ide", id="ylium_with_ide_centre"),
        pytest.param(
            "C1(=CC=CC=C1)[B-]1(OC2=[N+](C=CC=C2)O1)C1=CC=CC=C1",
            "2,2-diphenyl-4λ5-[1,3,4,2]dioxazaborolo[4,5-a]pyridin-4-ylium-2-uide",
            id="ylium_with_uide_centre_drops_its_indicated_hydrogen",
        ),
        pytest.param(
            "O[B-]1([N+]2=C(C3=C1C=CC=C3)NC3=C2C=CC=C3)O",
            "6,6-dihydroxy-6,11-dihydro-5λ5-[1,3]benzimidazolo[1,2-b][2,1]benzazaborol-5-ylium-6-uide",
            id="uide_centre_pairs_with_the_remaining_indicated_hydrogen_as_dihydro",
        ),
    ],
)
def test_ring_heteroatom_ylium_cations(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C(c1ccc(cc1)C(=O)[NH+])[NH+]", "benzene-1,4-bis(carboxamidylium)", id="carboxamidylium_centres"),
        pytest.param("O=[S+](=O)c1ccc(cc1)[S+](=O)=O", "(1,4-phenylene)bis(dioxo-λ6-sulfanylium)", id="sulfonylium_centres"),
        pytest.param("O=C(c1ccccc1C(=O)S[S+])S[S+]", "(benzene-1,2-dicarbonyl)bis(disulfanylium)", id="disulfanylium_centres_on_acyl_groups"),
    ],
)
def test_polycations_of_ylium_groups_on_one_skeleton(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC[S+](C)CCOCCOCC[S+](C)CC", "3,12-dimethyl-6,9-dioxa-3,12-dithiatetradecane-3,12-diium", id="two_cationic_centres_in_a_replacement_chain"),
        pytest.param("CC[S+](C)CCOCCOCCOCC", "12-methyl-3,6,9-trioxa-12-thiatetradecan-12-ium", id="one_cationic_centre_in_a_replacement_chain"),
        pytest.param("CC[SH+]CCOCCOCCOCC", "3,6,9-trioxa-12-thiatetradecan-12-ium", id="cationic_centre_without_substituent"),
        pytest.param("CCOP(=O)(C1CCCCC1)OCC[S+](C)CCCCCC", "4-cyclohexyl-8-methyl-4-oxo-3,5-dioxa-8-thia-4λ5-phosphatetradecan-8-ium", id="lambda_heterounit_with_oxo_and_ring_substituent"),
        pytest.param("C1CC[As+]2(C1)CCCC2", "5λ5-arsaspiro[4.4]nonan-5-ylium", id="spiro_atom_cation_by_the_lambda_convention"),
        pytest.param("C1CC[N+]2(C1)CCCC2", "5λ5-azaspiro[4.4]nonan-5-ylium", id="spiro_ammonium_by_the_lambda_convention"),
        pytest.param("C1CC2CC[S+]1CC2", "1λ4-thiabicyclo[2.2.2]octan-1-ylium", id="bridgehead_cation_by_the_lambda_convention"),
    ],
)
def test_cationic_centres_named_by_skeletal_replacement_and_lambda_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("N[n+]1ccccc1", "1-aminopyridin-1-ium", id="amino_on_the_cationic_nitrogen"),
        pytest.param("CN(C)[n+]1ccccc1", "1-(dimethylamino)pyridin-1-ium", id="substituted_amino_on_the_cationic_nitrogen"),
        pytest.param("NN[n+]1ccccc1", "1-hydrazinylpyridin-1-ium", id="hydrazinyl_on_the_cationic_nitrogen"),
        pytest.param("N[n+]1ccn(C)c1", "3-amino-1-methyl-1H-imidazol-3-ium", id="amino_beside_a_second_ring_nitrogen"),
        pytest.param("C[n+]1ccccc1NN", "2-hydrazinyl-1-methylpyridin-1-ium", id="hydrazinyl_on_a_ring_carbon"),
        pytest.param("C[n+]1ccc(NNC)cc1", "1-methyl-4-(2-methylhydrazin-1-yl)pyridin-1-ium", id="substituted_hydrazinyl_on_a_ring_carbon"),
        pytest.param("C[n+]1ccc(N=Nc2ccccc2)cc1", "1-methyl-4-(phenyldiazenyl)pyridin-1-ium", id="diazenyl_on_a_ring_carbon"),
    ],
)
def test_ring_cations_with_nitrogen_bearing_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-]Cl=O", "chlorite", id="chlorite"),
        pytest.param("[O-]Cl(=O)=O", "chlorate", id="chlorate"),
        pytest.param("[O-]Cl(=O)(=O)=O", "perchlorate", id="perchlorate"),
        pytest.param("[O-]Br(=O)=O", "bromate", id="bromate"),
        pytest.param("[O-]I(=O)(=O)=O", "periodate", id="periodate"),
        pytest.param("[Na+].[O-]Cl(=O)(=O)=O", "sodium perchlorate", id="halogen_oxoanion_in_a_metal_salt"),
        pytest.param("C[n+]1ccccc1.[O-]Cl(=O)(=O)=O", "1-methylpyridin-1-ium perchlorate", id="halogen_oxoanion_in_an_organic_salt"),
        pytest.param("[O-]C#N", "cyanate", id="cyanate"),
        pytest.param("[S-]C#N", "thiocyanate", id="thiocyanate"),
        pytest.param("[Se-]C#N", "selenocyanate", id="selenocyanate"),
        pytest.param("[K+].[S-]C#N", "potassium thiocyanate", id="cyanate_family_in_a_salt"),
    ],
)
def test_anions_of_halogen_oxoacids_and_cyanic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[CH2+]C[CH2+]", "propane-1,3-bis(ylium)", id="two_hydride_losses_on_a_chain"),
        pytest.param("C[C+2]C", "propane-2,2-bis(ylium)", id="two_hydride_losses_on_one_carbon"),
        pytest.param("C1=C[CH+][CH+]1", "cyclobut-3-ene-1,2-bis(ylium)", id="two_hydride_losses_on_a_ring"),
        pytest.param("O=[C+]CC[C+]=O", "1,4-dioxobutane-1,4-bis(ylium)", id="diacylium_cations_as_oxo_substituted_ylium_centres"),
        pytest.param("CC(=[OH+])CC(C)=[OH+]", "(pentane-2,4-diylidene)bis(oxidanium)", id="two_protonated_carbonyl_groups"),
        pytest.param("C[N+](C)(C)C(=O)c1ccccc1", "N,N,N-trimethylbenzamidium", id="quaternary_acylammonium"),
        pytest.param("CC(=O)[NH2+]C", "N-methylacetamidium", id="protonated_secondary_acylammonium"),
        pytest.param("C[N+](C)(C)C(=O)CC(=O)[N+](C)(C)C", "N1,N1,N1,N3,N3,N3-hexamethylpropanebis(amidium)", id="two_acylammonium_groups"),
        pytest.param("c1ccccc1C#[NH+]", "benzonitrilium", id="nitrilium"),
        pytest.param("C(CC#[NH+])C#[NH+]", "butanebis(nitrilium)", id="two_nitrilium_groups"),
    ],
)
def test_cationic_centres_on_characteristic_groups_and_hydride_losses(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C[n+]1ccccc1CC[n+]1ccccc1", id="units_differ"),
        pytest.param("C[n+]1ccc(OCc2cc[n+](C)cc2)cc1", id="unsymmetrical_linking_group"),
    ],
)
def test_cations_that_are_not_a_multiplicative_assembly(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "C[n+]1c(-c2ccccc2)cc(-c2ccccc2)cc1C(=O)[O-]",
            "1-methyl-4,6-diphenylpyridin-1-ium-2-carboxylate",
            id="ring_cation_with_anionic_suffix",
        ),
        pytest.param("C[n+]1ccn(C)c1C(=O)[O-]", "1,3-dimethyl-1H-imidazol-3-ium-2-carboxylate", id="indicated_hydrogen"),
        pytest.param("C[n+]1ccccc1CC(=O)[O-]", "(1-methylpyridin-1-ium-2-yl)acetate", id="ring_cation_as_prefix"),
        pytest.param(
            "O=C([O-])CCn1cc[n+](C)c1",
            "3-(3-methyl-1H-imidazol-3-ium-1-yl)propanoate",
            id="ring_cation_as_prefix_with_second_nitrogen",
        ),
        pytest.param("[CH2+]CC(=O)[O-]", "propan-3-ylium-1-oate", id="carbocation_with_chain_anion"),
        pytest.param("[O-]C(=O)CC[CH+]CC(=O)[O-]", "hexan-4-ylium-1,6-dioate", id="carbocation_with_two_anion_groups"),
        pytest.param("[CH2+]CS(=O)(=O)[O-]", "ethan-2-ylium-1-sulfonate", id="carbocation_with_sulfonate"),
    ],
)
def test_cation_on_the_parent_hydride_of_an_anionic_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC[Si]C[N+](C)(C)C.[I-]",
        "CC[Ge]C[N+](C)(C)C.[I-]",
    ],
)
def test_valence_deficient_heteroatom_is_not_named_as_its_hydride(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


# P-71.2.1.1, P-71.2.2.1, P-71.2.3 radicals on parent hydrides of other elements
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[NH2]", "azanyl", id="mononuclear_pnictogen"),
        pytest.param("[SiH3]", "silyl", id="group_14_ane_becomes_yl"),
        pytest.param("[BH2]", "boranyl", id="boron"),
        pytest.param("[OH]", "hydroxyl", id="retained_hydroxyl"),
        pytest.param("O[O]", "hydroperoxyl", id="retained_hydroperoxyl"),
        pytest.param("[SiH2]", "silylidene", id="divalent_centre"),
        pytest.param("[SiH]", "silylidyne", id="trivalent_centre"),
        pytest.param("C[Sn]C", "dimethylstannylidene", id="substituted_divalent_centre"),
        pytest.param("CC[P]C", "ethyl(methyl)phosphanyl", id="substituted_phosphorus"),
        pytest.param("[SiH3][SiH][SiH3]", "trisilan-2-yl", id="chain_locant"),
        pytest.param("CN(C)[N]", "dimethylhydrazinylidene", id="isodiazene_as_the_divalent_nitrogen"),
        pytest.param("C[N+](C)=[N-]", "dimethylhydrazinylidene", id="isodiazene_as_the_zwitterion"),
        pytest.param("CC[N+](C)=[N-]", "ethyl(methyl)hydrazinylidene", id="isodiazene_prefixes_alphabetized"),
        pytest.param("[NH2][NH]", "hydrazinyl", id="retained_hydrazine_root"),
        pytest.param("[NH][NH]", "hydrazine-1,2-diyl", id="two_radical_centres_on_a_chain"),
    ],
)
def test_radicals_on_heteroatom_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "CC[NH+](CC)Cc1ccc(-c2nc(C(=O)[O-])cs2)cc1",
            "2-{4-[(diethylazaniumyl)methyl]phenyl}-1,3-thiazole-4-carboxylate",
        ),
        ("C[NH+](C)Cc1nc(C(=O)[O-])cs1", "2-[(dimethylazaniumyl)methyl]-1,3-thiazole-4-carboxylate"),
        ("C[NH+](C)Cc1ccc(-c2ccc(C(=O)[O-])o2)cc1", "5-{4-[(dimethylazaniumyl)methyl]phenyl}furan-2-carboxylate"),
    ],
)
def test_ammonium_on_a_heterocyclic_anionic_parent_stays_an_azaniumyl_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[N+]1(CC(=O)c2ccccc2)CCCCC1", "1-methyl-1-(2-oxo-2-phenylethyl)piperidin-1-ium", id="ketone_side_chain_on_a_quaternary_ring_nitrogen"),
        pytest.param("C[N+](C)(C)C(=O)O", "carboxytri(methyl)azanium", id="acid_group_on_the_cationic_nitrogen"),
        pytest.param("C[N+](C)(C)C(=O)Cl", "carbonochloridoyltri(methyl)azanium", id="acid_halide_on_the_cationic_nitrogen"),
    ],
)
def test_onium_centres_that_carry_characteristic_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("COC(=O)c1cc[n+](Cc2ccccc2)cc1.[Br-]", "1-benzyl-4-(methoxycarbonyl)pyridin-1-ium bromide", id="ester_is_a_prefix_of_the_ring_cation"),
        pytest.param("COC(=O)c1cc[n+](Cc2ccccc2)cc1F.[Br-]", "1-benzyl-3-fluoro-4-(methoxycarbonyl)pyridin-1-ium bromide", id="halogen_beside_the_ester_prefix"),
    ],
)
def test_ring_cation_with_an_ester_group_names_the_cation_as_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=[S-]C1=CC=CC=C1", "oxo(phenyl)-λ4-sulfanide"),
        ("C(C)(=O)[NH+]", "acetamidylium"),
        ("C(C)[NH+]", "ethanaminylium"),
        ("[NH3+]CC(CC[NH3+])CC[NH3+]", "3-(azaniumylmethyl)pentane-1,5-bis(aminium)"),
        ("NCCC[NH3+]", "3-aminopropan-1-aminium"),
    ],
)
def test_ion_endings_hydrogen_free_cations_and_cationic_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[NH2+]", "azanylium", id="nitrenium"),
        pytest.param("[SiH3+]", "silylium", id="silylium"),
        pytest.param("C[Si+](C)C", "trimethylsilylium", id="trimethylsilylium"),
        pytest.param("C1(=CC=CC=C1)[Si+](C1=CC=CC=C1)C1=CC=CC=C1", "triphenylsilylium", id="triphenylsilylium"),
        pytest.param("[PH2+]", "phosphanylium", id="phosphanylium"),
        pytest.param("C[PH+]", "methylphosphanylium", id="substituted_phosphanylium"),
        pytest.param("C1(=CC=CC=C1)[S+]", "phenylsulfanylium", id="phenylsulfanylium"),
        pytest.param("C[O+]", "methoxylium", id="methoxylium"),
        pytest.param("c1ccccc1[O+]", "phenoxylium", id="phenoxylium"),
        pytest.param("ClCC(=O)[O+]", "(chloroacetyl)oxylium", id="acyloxylium"),
        pytest.param("[O+]CC[O+]", "(ethane-1,2-diyl)bis(oxylium)", id="bis_oxylium"),
        pytest.param("[S+]c1cccc([S+])n1", "(pyridine-2,6-diyl)bis(sulfanylium)", id="bis_sulfanylium"),
        pytest.param("C[N+]1(C)CC[N+]CC1", "4,4-dimethylpiperazin-4-ium-1-ylium", id="ring_ium_and_nitrenium_ylium"),
        pytest.param("[BH2+]", "boranylium", id="boranylium"),
        pytest.param("[Cl+]", "chloranylium", id="chloranylium"),
        pytest.param("C(C)=[OH+]", "ethylideneoxidanium", id="ylidene_oxonium"),
        pytest.param("C(C)[O+]=C(C)C", "ethyl(propan-2-ylidene)oxidanium", id="mixed_ylidene_and_alkyl_oxonium"),
        pytest.param("OC(C)=[OH+]", "(1-hydroxyethylidene)oxidanium", id="hydroxyethylidene_oxonium"),
        pytest.param("OC(C)=[N+](C)C", "(1-hydroxyethylidene)di(methyl)azanium", id="ylidene_azanium"),
        pytest.param("C(C)(=O)[Cl+]C", "acetyl(methyl)chloranium", id="chloronium_with_an_acyl_group"),
    ],
)
def test_ylium_and_onium_cations_of_mononuclear_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC(=O)C1=CC=C(C#[N+][O-])C=C1", "4-(methoxycarbonyl)benzonitrile oxide"),
        ("c1ccccc1C#[N+][S-]", "benzonitrile sulfide"),
        ("[O-][N+]#CC1CCC(CC1)C#[N+][O-]", "cyclohexane-1,4-dicarbonitrile dioxide"),
        ("O=N#CC1=CC=C(C(=O)[O-])C=C1", "4-[(oxo-λ5-azanylidyne)methyl]benzoate"),
        ("[O-]C(=O)c1ccc(cc1)C#[N+][O-].[Na+]", "sodium 4-[(oxo-λ5-azanylidyne)methyl]benzoate"),
    ],
)
def test_nitrile_oxides_and_their_prefix_beside_an_anion(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CN(C)C(N)=[NH2+]", "N,N-dimethylguanidinium"),
        ("CNC(=[O+]c1ccccc1)NC", "N,N'-dimethyl-O-phenyluronium"),
        ("CNC(=[S+]C)Nc1ccccc1", "N,S-dimethyl-N'-phenylthiouronium"),
    ],
)
def test_guanidinium_and_uronium_cations(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[c+]1ccccc1", "benzenylium"),
        ("Cc1cc[c+]cc1", "4-methylbenzen-1-ylium"),
        ("c1ccccc1C(=O)O[OH2+]", "2-benzoyldioxidan-1-ium"),
        ("[NH3+]N", "hydrazin-1-ium"),
        ("C[NH2+]NC", "1,2-dimethylhydrazin-1-ium"),
        ("[NH3+][NH3+]", "hydrazine-1,2-diium"),
        ("[OH2+]C(C)=[OH+]", "(1-oxidaniumylethylidene)oxidanium"),
    ],
)
def test_aryl_cations_chain_onium_cations_and_two_centre_oxonium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_added_hydrogen_of_a_mancude_ring_cation():
    assert smiles_to_iupac("C=1CC=C[C+]2C=C3C=CC=CC3=CC12") == "anthracen-4a(2H)-ylium"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[N+]1(C)CC[N-]CC1", "1,1-dimethylpiperazin-1-ium-4-ide"),
        ("C[Se+]1CC[N-]CC1", "1-methyl-1,4-selenomorpholin-1-ium-4-ide"),
    ],
)
def test_zwitterionic_ring_with_a_ring_cation_and_a_ring_anion(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
