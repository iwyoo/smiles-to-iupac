import pytest
from smiles_to_iupac import NonPreferredNameWarning, smiles_to_alternative_names, smiles_to_iupac
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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-]c1ccc(CC=C)cc1", "4-(prop-2-en-1-yl)phenoxide", id="phenoxide_unsaturated_ring_substituent_is_named"),
        pytest.param("[O-]CC[O-]", "ethane-1,2-bis(olate)", id="two_alkoxide_groups_is_named"),
        pytest.param("[O-]CCOC", "2-methoxyethan-1-olate", id="ether_oxygen_alongside_alkoxide_is_named"),
        pytest.param("C(S(=O)(=O)Cl)CC[O-]", "3-(chlorosulfonyl)propan-1-olate", id="acid_halide_below_the_anion_is_a_prefix"),
        pytest.param("C([NH3+])CCOC#N", "3-cyanatopropan-1-aminium", id="cyanate_below_the_cation_is_a_prefix"),
        pytest.param("C([N+](C)(C)C)CCSC#N", "N,N,N-trimethyl-3-(thiocyanato)propan-1-aminium", id="thiocyanate_below_a_quaternary_cation_is_a_prefix"),
        pytest.param("C(C(=S)N)CC[O-]", "3-carbamothioylpropan-1-olate", id="thioamide_below_the_anion_is_a_prefix"),
        pytest.param("C([O-])CCSC#N", "3-(thiocyanato)propan-1-olate", id="thiocyanate_below_the_anion_is_a_prefix"),
        pytest.param("CC[C@@H](Cl)[O-]", "(1R)-1-chloropropan-1-olate", id="alkoxide_stereocenter_with_coexisting_halogen"),
        pytest.param("[O-]C(O)c1ccccc1", "hydroxy(phenyl)methanolate", id="neutral_hydroxy_on_the_alkoxide_carbon_is_a_prefix"),
        pytest.param("CC[C@H](C)CC[NH-]", "(3S)-3-methylpentan-1-aminide", id="primary_aminide_stereocenter"),
        pytest.param("CC[C@H](C)C[CH2-]", "(3R)-3-methylpentan-1-ide", id="carbanide_stereocenter_ranks_the_anionic_carbon_below_its_twin"),
    ],
)
def test_phenoxide_unsaturated_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COC(=O)CC[O-].[Na+]", "sodium 2-(methoxycarbonyl)ethan-1-olate", id="alkoxycarbonyl_prefix_on_olate"),
        pytest.param("CC(=O)OCC[O-].[Na+]", "sodium 2-(acetyloxy)ethan-1-olate", id="acyloxy_prefix_on_olate"),
        pytest.param("COC(=O)C(C)[S-].[K+]", "potassium 1-(methoxycarbonyl)ethane-1-thiolate", id="ester_prefix_on_thiolate"),
        pytest.param("C[C@H]([O-])C(=O)OC.[Na+]", "sodium (1S)-1-(methoxycarbonyl)ethan-1-olate", id="ester_prefix_with_stereodescriptor"),
        pytest.param("COC(=O)CC[N-]C.[Na+]", "sodium 2-(methoxycarbonyl)-N-methylethan-1-aminide", id="ester_prefix_on_secondary_aminide"),
    ],
)
def test_anion_outranks_ester(smiles, expected):
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
        pytest.param("C=C[NH-]", "ethenaminide", id="enamine_aminide_is_named"),
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
        pytest.param("C[N+](C)(C)CCC(=O)OC", "3-methoxy-N,N,N-trimethyl-3-oxopropan-1-aminium", id="cation_outranks_ester"),
        pytest.param("CCOC(=O)CC[NH3+]", "3-ethoxy-3-oxopropan-1-aminium", id="ester_carbon_in_primary_aminium_chain"),
        pytest.param("CC(=O)OCC[N+](C)(C)C", "2-(acetyloxy)-N,N,N-trimethylethan-1-aminium", id="acyloxy_prefix_on_aminium"),
        pytest.param("CCOS(=O)(=O)CC[N+](C)(C)C", "2-(ethoxysulfonyl)-N,N,N-trimethylethan-1-aminium", id="sulfonate_ester_prefix_on_aminium"),
        pytest.param("CC(=O)OCC[n+]1ccccc1", "1-[2-(acetyloxy)ethyl]pyridin-1-ium", id="acyloxy_prefix_on_ring_cation"),
        pytest.param("COCC[n+]1cccc(CC)c1", "3-ethyl-1-(2-methoxyethyl)pyridin-1-ium", id="ether_in_the_substituent_of_a_ring_cation"),
        pytest.param("CCOS(=O)(=O)CC[NH+]1CCCCC1", "1-[2-(ethoxysulfonyl)ethyl]piperidin-1-ium", id="sulfonate_ester_prefix_on_saturated_ring_cation"),
        pytest.param("CCOC(=O)CC[NH+]1CCCC1", "1-(3-ethoxy-3-oxopropyl)pyrrolidin-1-ium", id="ester_carbon_in_chain_on_saturated_ring_cation"),
        pytest.param("C[N+](C)(C)c1ccc(C(=O)O)cc1", "4-carboxy-N,N,N-trimethylanilinium", id="aryl_ammonium_with_acid"),
        pytest.param("Nc1ccccc1[NH3+]", "2-aminoanilinium", id="neutral_amino_prefix_on_aniline"),
        pytest.param("NC1=CC=CC=C1[NH2+]C", "2-amino-N-methylanilinium", id="neutral_amino_prefix_n_substituted"),
        pytest.param("Nc1cccnc1[NH3+]", "3-aminopyridin-2-aminium", id="neutral_amino_prefix_on_heteroring"),
        pytest.param(
            "C[N+](C)(CCNC=O)CCNC=O", "2-formamido-N-(2-formamidoethyl)-N,N-dimethylethan-1-aminium", id="two_formamido_arms"
        ),
        pytest.param(
            "C[N+](C)(CCNC(C)=O)CCNC(C)=O", "2-acetamido-N-(2-acetamidoethyl)-N,N-dimethylethan-1-aminium", id="two_acetamido_arms"
        ),
        pytest.param(
            "C[N+](C)(CCNC(=O)c1ccccc1)CCNC(=O)c1ccccc1",
            "2-benzamido-N-(2-benzamidoethyl)-N,N-dimethylethan-1-aminium",
            id="two_benzamido_arms",
        ),
        pytest.param("C[NH+](CCNC=O)CCNC=O", "2-formamido-N-(2-formamidoethyl)-N-methylethan-1-aminium", id="two_arms_on_a_tertiary_ammonium"),
        pytest.param("CC(=O)NCC[NH2+]CCNC(C)=O", "2-acetamido-N-(2-acetamidoethyl)ethan-1-aminium", id="two_arms_on_a_secondary_ammonium"),
        pytest.param("C[N+](C)(C)c1ccc(C=NO)cc1", "4-(N-hydroxymethanimidoyl)-N,N,N-trimethylanilinium", id="oxime_on_an_aryl_ammonium"),
        pytest.param("C[N+](C)(C)c1ccccc1C=NO", "2-(N-hydroxymethanimidoyl)-N,N,N-trimethylanilinium", id="oxime_ortho_to_the_ammonium"),
        pytest.param("C[N+](C)(C)c1ccc(/C=N/O)cc1", "4-[(E)-N-hydroxymethanimidoyl]-N,N,N-trimethylanilinium", id="e_oxime_on_an_aryl_ammonium"),
        pytest.param("C[N+](C)(C)c1ccc(/C=N\\O)cc1", "4-[(Z)-N-hydroxymethanimidoyl]-N,N,N-trimethylanilinium", id="z_oxime_on_an_aryl_ammonium"),
        pytest.param("C[N+](C)(C)c1ccc(C=N)cc1", "4-(methanimidoyl)-N,N,N-trimethylanilinium", id="imine_on_an_aryl_ammonium"),
        pytest.param("C[N+](C)(C)c1ccc(C=NC)cc1", "N,N,N-trimethyl-4-(N-methylmethanimidoyl)anilinium", id="alphanumerical_order_with_an_imino_prefix"),
        pytest.param("C[N+](C)(C)C1CCCCC1=NO", "2-(hydroxyimino)-N,N,N-trimethylcyclohexan-1-aminium", id="oxime_on_a_ring_ammonium"),
        pytest.param("C[N+](C)(C)C/C=N/O", "(2E)-2-(hydroxyimino)-N,N,N-trimethylethan-1-aminium", id="e_oxime_on_a_chain_ammonium"),
        pytest.param("C[N+](C)(C)C/C=N\\O", "(2Z)-2-(hydroxyimino)-N,N,N-trimethylethan-1-aminium", id="z_oxime_on_a_chain_ammonium"),
        # P-41, P-73.1.1.2: the cation is the parent, so thiocarbonyl groups are prefixes as their oxygen analogues are
        pytest.param("C[N+](C)(C)CCC(N)=S", "3-amino-N,N,N-trimethyl-3-sulfanylidenepropan-1-aminium", id="thioamide_carbon_on_a_chain_ammonium"),
        pytest.param("S=C(C)NCC[NH3+]", "2-(ethanethioamido)ethan-1-aminium", id="thioamide_nitrogen_on_a_chain_ammonium"),
        pytest.param("S=C(C)N1CC[NH2+]CC1", "4-(ethanethioyl)piperazin-1-ium", id="thioacyl_on_a_ring_cation_nitrogen"),
        pytest.param("O=C(c1ccoc1)N1CC[NH2+]CC1", "4-(furan-3-carbonyl)piperazin-1-ium", id="heteroaroyl_on_a_ring_cation_nitrogen"),
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
        pytest.param("C[CH+]/C=C/C", "(3E)-pent-3-en-2-ylium", id="double_bond_descriptor_of_a_carbenium_parent"),
        pytest.param("C[C@@H](CC)[CH+]C", "(3R)-3-methylpentan-2-ylium", id="stereocentre_descriptor_of_a_carbenium_parent"),
    ],
)
def test_benzenide_name_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
            "Cl[Ti+](Cl)Cl.c1cc[cH-]c1.O",
            "trichloridotitanium(1+) cyclopenta-2,4-dien-1-ide—water (1/1)",
            id="metal_cation_cyclopentadienide_salt_with_neutral_component",
        ),
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
        pytest.param("c1cc[n+](cc1)[N-]C", "N-(pyridin-1-ium-1-yl)methanaminide", id="ring_nitrogen_aminide_ylide"),
        pytest.param("c1cc[n+](cc1)[S-]", "(pyridin-1-ium-1-yl)sulfanide", id="ring_nitrogen_sulfanide_ylide"),
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
            "C[BH2-][NH2+][BH2-][NH2+]C",
            "{[(methanaminiumyl)boranuidyl]azaniumyl}(methyl)boranuide",
            id="alternating_ammonium_boranuide_chain",
        ),
        pytest.param(
            "C[BH2-][NH2+][BH3-]", "[(boranuidyl)azaniumyl](methyl)boranuide", id="alternating_chain_with_a_spare_boranuide"
        ),
        pytest.param(
            "C[Se+]1CC2=C(C([CH-]1)=O)C=CC=C2",
            "2-methyl-4-oxo-3,4-dihydro-1H-2-benzoselenopyran-2-ium-3-ide",
            id="ring_chalcogenium_with_ring_carbanide",
        ),
        pytest.param("[O-]c1cc[n+](C)cc1", "1-methylpyridin-1-ium-4-olate", id="ring_cation_with_an_anionic_group_suffix"),
        pytest.param("C[N+]1(C)C=CC([O-])=N1", "1,1-dimethyl-1H-pyrazol-1-ium-3-olate", id="olate_on_a_ring_carbon_double_bonded_to_a_ring_nitrogen"),
        pytest.param("C[N+]1(C)C=CC([S-])=N1", "1,1-dimethyl-1H-pyrazol-1-ium-3-thiolate", id="thiolate_beside_a_disubstituted_ring_nitrogen"),
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
        pytest.param("c1ccccc1/N=[N+](\\[O-])c1ccccc1", "(Z)-diphenyldiazene oxide", id="azoxy_descriptor_ranks_the_oxide_oxygen"),
        pytest.param("CC/N=[N+](/[O-])C", "(E)-1-ethyl-2-methyldiazene 2-oxide", id="azoxy_descriptor_with_the_oxide_locant"),
        pytest.param("Clc1ccccc1N=[N+]([O-])c1ccccc1", "1-(2-chlorophenyl)-2-phenyldiazene 2-oxide", id="azoxy_with_the_oxide_locant"),
        pytest.param("[O-]C(=O)c1ccc(cc1)N=[N+]([O-])c1ccccc1", "4-(2-oxo-2-phenyl-2λ5-diazenyl)benzoate", id="azoxy_prefix_beside_a_carboxylate"),
        pytest.param("[O-]C(=O)CCN=[N+]([O-])c1ccccc1", "3-(2-oxo-2-phenyl-2λ5-diazenyl)propanoate", id="azoxy_prefix_on_a_chain_anion"),
        pytest.param("[O-]C(=O)c1ccc(cc1)[N+]([O-])=Nc1ccccc1", "4-(1-oxo-2-phenyl-1λ5-diazenyl)benzoate", id="azoxy_prefix_attached_through_the_oxide_nitrogen"),
        pytest.param("C[N+](C)(C)c1ccc(cc1)N=[N+]([O-])c1ccccc1", "N,N,N-trimethyl-4-(2-oxo-2-phenyl-2λ5-diazenyl)anilinium", id="azoxy_prefix_beside_an_ammonium"),
        pytest.param("CCC=[S+][O-]", "propylidene-\u03bb4-sulfanone", id="thioaldehyde_s_oxide"),
        pytest.param("[O-][S+]1C=CC=C1", "1H-1\u03bb4-thiophen-1-one", id="heterone_s_oxide_of_a_ring"),
        pytest.param("[O-][S+]1N=CC=C1", "1H-1\u03bb4,2-thiazol-1-one", id="heterone_s_oxide_with_two_heteroatoms"),
        pytest.param("[O-][Se+]1C=CC=C1", "1H-1\u03bb4-selenophen-1-one", id="selenium_analogue_of_the_ring_s_oxide"),
        pytest.param("C[S+]([O-])C", "(methanesulfinyl)methane", id="sulfoxide_written_with_charges"),
        pytest.param("CC(C)=[N+](C)[O-]", "N-methylpropan-2-imine N-oxide", id="nitrone"),
        pytest.param("C[N+](C)([O-])CC(=O)O", "carboxy-N,N-dimethylmethanamine N-oxide", id="amine_oxide_outranks_a_carboxylic_acid"),
        pytest.param("C[N+](C)([O-])c1ccc(O)cc1", "4-hydroxy-N,N-dimethylaniline N-oxide", id="aniline_oxide_outranks_a_hydroxy_group"),
        pytest.param("C[N+]([O-])=CCO", "2-hydroxy-N-methylethan-1-imine N-oxide", id="nitrone_outranks_a_hydroxy_group"),
        pytest.param("C[N+]([O-])=CCC#N", "2-cyano-N-methylethan-1-imine N-oxide", id="nitrone_outranks_a_nitrile"),
        pytest.param("O=C(O)c1ccc(C=[N+](C)[O-])cc1", "1-(4-carboxyphenyl)-N-methylmethanimine N-oxide", id="nitrone_outranks_a_carboxylic_acid"),
        pytest.param("C/C=[N+](\\C)[O-]", "(1E)-N-methylethanimine N-oxide", id="nitrone_descriptor_ranks_the_oxide_ligand"),
        pytest.param("C/C=[N+](/[O-])C", "(1E)-N-methylethanimine N-oxide", id="nitrone_oxide_written_before_the_substituent"),
        pytest.param("C/C=[N+](/c1ccccc1)[O-]", "(1Z)-N-phenylethanimine N-oxide", id="nitrone_with_an_aryl_nitrogen_substituent"),
        pytest.param("[O-]/[NH+]=C\\c1ccccc1", "(Z)-phenylmethanimine N-oxide", id="nitrone_hydrogen_on_nitrogen_ranks_the_oxide"),
        pytest.param("[O-]/[NH+]=C/CC", "(1E)-propan-1-imine N-oxide", id="nitrone_hydrogen_on_nitrogen_with_a_locant"),
        pytest.param("CC(C)=[O+][N-]C", "N-[(propan-2-ylidene)oxidaniumyl]methanaminide", id="carbonyl_imide"),
        pytest.param("CC(C)=[O+][C-](C)C", "2-[(propan-2-ylidene)oxidaniumyl]propan-2-ide", id="carbonyl_ylide"),
        pytest.param("CC(C)=[N+](C)[C-](C)C", "2-(N-methylpropan-2-iminiumyl)propan-2-ide", id="azomethine_ylide"),
        pytest.param("CC#[N+][O-]", "acetonitrile oxide", id="nitrile_oxide"),
        pytest.param("CC#[N+][C-](C)C", "2-(acetonitriliumyl)propan-2-ide", id="nitrile_ylide"),
        pytest.param("[CH2]C([CH2])c1ccc(cc1)[P+]([O-])(c1ccccc1)c1ccccc1", "2-{4-[oxidodi(phenyl)phosphaniumyl]phenyl}propane-1,3-diyl", id="dipolar_phosphorus_substituent_group"),
        pytest.param("[CH2]C([CH2])c1ccc(cc1)CC#[N+][O-]", "2-{4-[2-(oxidoazaniumylidyne)ethyl]phenyl}propane-1,3-diyl", id="dipolar_nitrile_oxide_substituent_group"),
        pytest.param("[CH2]C[N+]([O-])(C)C", "2-[dimethyl(oxo)-\u03bb5-azanyl]ethyl", id="amine_oxide_substituent_group_as_oxo_azanyl"),
        pytest.param("C[B-](C)(C)c1ccc([CH2])cc1", "[4-(trimethylboranuidyl)phenyl]methyl", id="ate_centre_cited_as_prefix_of_a_carbon_radical"),
        pytest.param("C[N]CC[B-](C)(C)C", "N-methyl-2-(trimethylboranuidyl)ethan-1-aminyl", id="ate_centre_cited_as_prefix_of_an_aminyl_radical"),
    ],
)
def test_dipolar_compounds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [("[NH3+]C#N", "cyanoazanium"), ("N#[N+]C#N", "cyanodiazyn-1-ium")],
)
def test_cation_on_a_nitrile_carbon(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereogenic_azanium_with_a_heteroatom_ligand():
    assert smiles_to_iupac("C[N@+](CC)(CCC)Cl") == "(S)-chloro(ethyl)(methyl)(propyl)azanium"


def test_unsupported_base_fragment():
    assert (
        smiles_to_iupac("c1ccccc1CC(=O)Nc1ccccc1C(=O)OCCCC.Cl")
        == "butyl 2-(2-phenylacetamido)benzoate—hydrogen chloride (1/1)"
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].[Fe+3]", "iron(3+) triacetate", id="charge_number_after_metal"),
        pytest.param("[CH3-].[CH3-].[Cu+2]", "copper(2+) dimethanide", id="carbanide_anion"),
        pytest.param("[Fe+3].[Fe+3].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-]", "diiron(3+) tris(sulfate)", id="multiplied_cation"),
        pytest.param("[Na+].[Cu+2].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-]", "copper(2+) sodium triacetate", id="cations_alphabetical"),
        pytest.param("CC(=O)[O-].[Ag+]", "silver acetate", id="single_common_charge_without_number"),
    ],
)
def test_transition_metal_cation_salt_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_same_cation_type_prefers_multiplying_prefix():
    assert smiles_to_iupac("[Na+].[Na+].[O-]C(=O)[O-]") == "disodium carbonate"


def test_pnictogen_cation_is_not_a_radical_centre():
    acid = "[O-]C(=O)CCC(=O)O"
    assert smiles_to_iupac(f"[Sb+3].{acid}.{acid}.{acid}") == "antimony tris(3-carboxypropanoate)"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[Na+].[Cl-].[Na+].[Cl-].O", "sodium chloride—water (2/1)", id="repeated_formula_unit_of_the_salt"),
        pytest.param("[CH3-].[HH].[Na+]", "sodium methanide—dihydrogen (1/1)", id="organic_salt_with_dihydrogen"),
        pytest.param("C.[Al+3]", "methane—aluminium(3+) (1/1)", id="cation_with_a_neutral_molecule"),
        pytest.param("[CH3-].[Ta]", "methanide—tantalum (1/1)", id="anion_with_a_neutral_atom"),
        pytest.param("[Al+3].[K+].[O-]S(=O)(=O)[O-]", "aluminium(3+)—potassium(1+)—sulfate (1/1/1)", id="ions_that_do_not_balance"),
    ],
)
def test_adducts_of_ions_with_other_components(smiles, expected):
    with pytest.warns(NonPreferredNameWarning, match="P-14.8.2"):
        assert smiles_to_iupac(smiles) == expected


def test_radical_beside_an_ion_is_not_named_as_a_closed_shell_adduct():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH3].[Na+]")


def test_oxidanium():
    assert smiles_to_iupac("[OH3+]") == "oxidanium"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("[AsH4+]", "arsanium"),
        ("[SeH3+]", "selanium"),
        ("[TeH3+]", "tellanium"),
        ("C[Te+](C)C.[Cl-]", "trimethyltellanium chloride"),
        ("C[Te+](C)CC(=O)[O-]", "(dimethyltellaniumyl)acetate"),
        ("[ClH2+]", "chloranium"),
    ],
)
def test_onium_cations_of_the_other_mononuclear_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)[OH2+]", "(propan-2-yl)oxidanium", id="branched_substituent"),
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
    "smiles, expected",
    [
        pytest.param("C/C=C\\C[P+](C)(C)C", "[(2Z)-but-2-en-1-yl]tri(methyl)phosphanium", id="phosphonium"),
        pytest.param("C/C=C/C[S+](C)C", "[(2E)-but-2-en-1-yl]di(methyl)sulfanium", id="sulfonium"),
        pytest.param("C/C=C/C[As+](C)(C)C", "[(2E)-but-2-en-1-yl]tri(methyl)arsanium", id="arsonium"),
        pytest.param("C/C=C/[SH2+]", "[(1E)-prop-1-en-1-yl]sulfanium", id="unsubstituted_sulfonium"),
        pytest.param("C/C=C/C[S+](C)C.[Br-]", "[(2E)-but-2-en-1-yl]di(methyl)sulfanium bromide", id="salt"),
    ],
)
def test_alkenyl_group_on_an_onium_centre_cites_its_double_bond_descriptor(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Na+].[O-][N+](=O)[O-]", "sodium nitrate"),
        ("N[Co+3](N)(N)(N)(N)N.[O-]C(C)=O.[O-]C(C)=O.[O-]C(C)=O", "hexaazanidocobalt(3+) triacetate"),
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
        pytest.param("CS(=O)(=O)[NH]", "methanesulfonamidyl", id="sulfonamide_radical"),
        pytest.param("CC(=O)N[NH]", "2-acetylhydrazin-1-yl", id="hydrazine_radical_with_acyl"),
        pytest.param("[NH]S(=O)(=O)CCS(=O)(=O)[NH]", "ethane-1,2-bis(sulfonamidyl)", id="two_sulfonamide_radicals"),
    ],
)
def test_nitrogen_radical_beside_a_heteroatom_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[CH2][C@H](Cl)CC", "(2S)-2-chlorobutyl", id="stereocentre_of_a_radical_group"),
        pytest.param("C/C=C/C(=O)[CH2]", "(3E)-2-oxopent-3-en-1-yl", id="double_bond_of_a_radical_group"),
        pytest.param("[CH2][C@H](Cl)C[CH2]", "(2S)-2-chlorobutane-1,4-diyl", id="stereocentre_of_a_polyradical"),
        pytest.param("[CH2]C(C)/C=C/C(C)[CH2]", "(3E)-2,5-dimethylhex-3-ene-1,6-diyl", id="double_bond_between_two_branched_carbons"),
        pytest.param("[CH2][C@H](C)CC", "(2R)-2-methylbutyl", id="centre_stereogenic_only_through_the_radical_centre"),
        pytest.param("[CH2][C@@H](C)CC", "(2S)-2-methylbutyl", id="opposite_configuration_of_the_same_group"),
    ],
)
def test_stereo_elements_of_a_radical_are_cited(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polyradical_stereocentre_is_never_dropped():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2][C@H](C)C[CH2]")


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1-c1[c]cccc1", "[1,1'-biphenyl]-2-yl", id="radical_centre_of_a_bare_ring_assembly"),
        pytest.param("Brc1ccccc1-c1[c]cccc1", "2'-bromo[1,1'-biphenyl]-2-yl", id="halogen_beside_the_radical_centre_of_a_ring_assembly"),
        pytest.param("c1ccccc1-c1ccccc1-c1[c]cccc1", "[1\u00b9,2\u00b9:2\u00b2,3\u00b9-terphenyl]-1\u00b2-yl", id="radical_centre_of_a_ring_assembly_chain"),
    ],
)
def test_radical_centre_of_a_ring_assembly_is_cited(smiles, expected):
    # P-71.2.1.2, P-28.2
    assert smiles_to_iupac(smiles) == expected


def test_two_radical_centers_diyl_name():
    assert smiles_to_iupac("[CH2][CH2]") == "ethane-1,2-diyl"


def test_three_radical_centers_name():
    assert smiles_to_iupac("[CH2][CH][CH2]") == "propane-1,2,3-triyl"


def test_five_ring_valences_beside_a_cationic_substituent():
    assert (
        smiles_to_iupac("[CH2+][C]1[C](C)[C](C)[C](C)[C]1C")
        == "1,2,3,4-tetramethyl-5-methyliumylcyclopentane-1,2,3,4,5-pentayl"
    )


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccc2c(c1)[CH]c1ccccc1C2", "anthracene-9-hydryl", id="fused_hydrocarbon"),
        pytest.param("c1ccc2c(c1)[CH]c1ccccc1C2C", "9-methylanthracene-9-hydryl", id="hydrogen_added_to_a_substituted_atom"),
        pytest.param("Clc1ccc2c(c1)[CH]c1ccccc1C2", "3-chloroanthracene-9-hydryl", id="lowest_locant_to_the_hydryl_atom"),
        pytest.param("[CH]1C=CC2=CC=CC=C2C1", "naphthalene-1-hydryl", id="radical_beside_the_added_hydrogen"),
        pytest.param("[CH]1C=CC=C2C=CC=CC12", "naphthalene-4a-hydryl", id="hydrogen_added_to_a_fusion_atom"),
        pytest.param("c1ccc2c(c1)[CH]c1ccccc1N2", "acridine-10-hydryl", id="hydrogen_added_to_a_ring_nitrogen"),
        pytest.param("[N]1C=CC=CC1", "pyridine-2-hydryl", id="radical_on_a_ring_nitrogen"),
        pytest.param("C1=CSC[CH]1", "thiophene-2-hydryl", id="monocycle_with_a_lone_pair_atom"),
        pytest.param("CC1C=CC=C[CH]1", "1-methylbenzene-1-hydryl", id="substituted_benzene"),
        pytest.param("[CH]1C=CC=CC1", "cyclohexa-2,4-dien-1-yl", id="every_position_equivalent"),
        pytest.param("C1CC2=CC=CC=C2[CH]1", "2,3-dihydro-1H-inden-1-yl", id="indicated_hydrogen_parent"),
        pytest.param("c1ccc2c(c1)[CH]c1ccccc1C2(C)C", "10,10-dimethyl-9,10-dihydroanthracen-9-yl", id="no_hydrogen_to_receive"),
    ],
)
def test_hydryl_radicals_of_mancude_ring_systems(smiles, expected):
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
    assert smiles_to_iupac("[CH2-][C@@H](C)CCC.[Li+]") == "lithium (2S)-2-methylpentan-1-ide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Valine zwitterion (branched chain).
        ("CC(C)C(C(=O)[O-])[NH3+]", "2-azaniumyl-3-methylbutanoate"),
        ("C[N+](C)(C)CCS(=O)(=O)[O-]", "2-(N,N-dimethylmethanaminiumyl)ethane-1-sulfonate"),
        pytest.param("C[S+](C)CC(=O)[O-]", "(dimethylsulfaniumyl)acetate", id="sulfonium_prefix_on_a_carboxylate"),
        pytest.param(
            "C[S+](C)CC([NH3+])C(=O)[O-]", "2-azaniumyl-3-(dimethylsulfaniumyl)propanoate", id="cationic_prefixes_of_two_classes_on_a_carboxylate"
        ),
        pytest.param(
            "NC(=[NH2+])NCCC[C@H](N)C(=O)[O-]", "(2S)-2-amino-5-{[amino(azaniumylidene)methyl]amino}pentanoate", id="guanidinium_prefix_on_a_nitrogen_of_a_carboxylate"
        ),
        pytest.param("NC(=[NH2+])c1ccc(cc1)C(=O)[O-]", "4-[amino(azaniumylidene)methyl]benzoate", id="amidinium_prefix_on_a_ring_of_a_carboxylate"),
        pytest.param(
            "[NH3+][C@@H](Cc1ccc[nH+]c1)C(=O)[O-]", "(2S)-2-azaniumyl-3-(pyridin-1-ium-3-yl)propanoate", id="ium_ring_beside_an_azaniumyl_prefix_on_a_carboxylate"
        ),
        pytest.param("C[P+](C)(C)CC(=O)[O-]", "(trimethylphosphaniumyl)acetate", id="phosphonium_prefix_on_a_carboxylate"),
        pytest.param("C[S+](C)CCS(=O)(=O)[O-]", "2-(dimethylsulfaniumyl)ethane-1-sulfonate", id="sulfonium_prefix_on_a_sulfonate"),
        pytest.param("C[S+](C)CCCCCC[N+](C)(C)C", "6-(dimethylsulfaniumyl)-N,N,N-trimethylhexan-1-aminium", id="senior_nitrogen_cation_is_the_parent"),
        pytest.param("C[S+](C)CCCCCC[P+](C)(C)C", "[6-(dimethylsulfaniumyl)hexyl]tri(methyl)phosphanium", id="phosphorus_outranks_sulfur"),
        pytest.param("C[O+](C)CCCC[N+](C)(C)C", "4-(dimethyloxidaniumyl)-N,N,N-trimethylbutan-1-aminium", id="oxygen_cation_prefix"),
        # P-74.1.3, P-72.2.2.2.1.2: the ammonium of an acid-ester anion is a prefix of the alkyl group, the anion name ends the compound
        pytest.param("C[N+](C)(C)CCOS(=O)(=O)[O-]", "2-(trimethylazaniumyl)ethyl sulfate", id="sulfate_ester"),
        pytest.param("[NH3+]CCOS(=O)(=O)[O-]", "2-azaniumylethyl sulfate", id="primary_ammonium_on_a_sulfate"),
        pytest.param("C[N+](C)(C)CCOS(=O)[O-]", "2-(trimethylazaniumyl)ethyl sulfite", id="sulfite_ester"),
        pytest.param("C[N+](C)(C)CCOC(=O)[O-]", "2-(trimethylazaniumyl)ethyl carbonate", id="carbonate_ester"),
        pytest.param("C[N+](C)(C)CCOP([O-])OC", "methyl 2-(trimethylazaniumyl)ethyl phosphite", id="phosphite_ester"),
        pytest.param("CP(=O)([O-])OCC[N+](C)(C)C", "2-(trimethylazaniumyl)ethyl methylphosphonate", id="methylphosphonate_ester"),
        pytest.param("O=S(=O)([O-])OCC[n+]1cccc(C)c1", "2-(3-methylpyridin-1-ium-1-yl)ethyl sulfate", id="ring_cation_prefix_on_a_sulfate"),
        pytest.param("O=P([O-])(OC)OCC[n+]1ccccc1", "methyl 2-(pyridin-1-ium-1-yl)ethyl phosphate", id="ring_cation_prefix_on_a_phosphate"),
        pytest.param("C[n+]1ccccc1C[P-]C", "methyl[(1-methylpyridin-1-ium-2-yl)methyl]phosphanide", id="ring_cation_prefix_on_a_phosphanide"),
        pytest.param("O=S(=O)([O-])OCC[P+](C)(C)C", "2-(trimethylphosphaniumyl)ethyl sulfate", id="phosphonium_prefix_on_a_sulfate"),
        pytest.param("CCP(=O)([O-])OCC[N+](C)(C)C", "2-(trimethylazaniumyl)ethyl ethylphosphonate", id="ethylphosphonate_ester"),
        pytest.param("CC(C)P(=O)([O-])OCC[N+](C)(C)C", "2-(trimethylazaniumyl)ethyl (propan-2-yl)phosphonate", id="branched_phosphonate_group"),
        pytest.param("C[N+](C)(C)CCOP(=O)([O-])c1ccccc1", "2-(trimethylazaniumyl)ethyl phenylphosphonate", id="phenylphosphonate_ester"),
        pytest.param("CP(=O)([O-])OCC[NH3+]", "2-azaniumylethyl methylphosphonate", id="primary_ammonium_on_a_phosphonate"),
        pytest.param("CP(=O)([O-])OCC[NH+](C)C", "2-(dimethylazaniumyl)ethyl methylphosphonate", id="tertiary_ammonium_on_a_phosphonate"),
        pytest.param("ClCCOS(=O)(=O)[O-]", "2-chloroethyl sulfate", id="substituted_group_of_a_plain_sulfate_anion"),
        # P-103.2.4, P-72.7: the carboxylate ends the compound and the other residue is an acyl prefix with its azaniumyl group
        pytest.param("[NH3+]CC(=O)N1CCCC1C(=O)[O-]", "1-(azaniumylacetyl)pyrrolidine-2-carboxylate", id="proline_c_terminus"),
        pytest.param("[NH3+]CC(=O)N1CCC[C@H]1C(=O)[O-]", "(2S)-1-(azaniumylacetyl)pyrrolidine-2-carboxylate", id="proline_with_a_descriptor"),
        pytest.param(
            "C[C@@H]([NH3+])C(=O)N1CCC[C@H]1C(=O)[O-]",
            "(2S)-1-[(2R)-2-azaniumylpropanoyl]pyrrolidine-2-carboxylate",
            id="descriptors_of_the_acyl_group_and_the_ring",
        ),
        pytest.param("CC([NH3+])C(=O)N1CCCC1C(=O)[O-]", "1-(2-azaniumylpropanoyl)pyrrolidine-2-carboxylate", id="alanyl_proline_unspecified"),
        pytest.param(
            "CC(C)[C@H]([NH3+])C(=O)NCC(=O)[O-]", "{[(2S)-2-azaniumyl-3-methylbutanoyl]amino}acetate", id="valine_n_terminus"
        ),
        pytest.param("CC(C)C([NH3+])C(=O)NCC(=O)[O-]", "[(2-azaniumyl-3-methylbutanoyl)amino]acetate", id="valine_n_terminus_unspecified"),
        pytest.param(
            "CC(C)C[C@H]([NH3+])C(=O)NCC(=O)[O-]", "{[(2S)-2-azaniumyl-4-methylpentanoyl]amino}acetate", id="leucine_n_terminus"
        ),
        pytest.param(
            "CC(C)[C@H]([NH3+])C(=O)N[C@@H](C)C(=O)[O-]",
            "(2S)-2-{[(2S)-2-azaniumyl-3-methylbutanoyl]amino}propanoate",
            id="descriptors_of_both_residues",
        ),
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
    ("C=C[S-]", "ethenethiolate"),
    ("P=C[S-]", "phosphanylidenemethanethiolate"),
    ("CN(C)[O-]", "dimethylaminoxide"),
    ("C/C=C\\C(C)N(C)[O-]", "methyl[(3Z)-pent-3-en-2-yl]aminoxide"),
    ("CO[O-]", "methaneperoxolate"),
    ("CCS[O-]", "ethane(SO-thioperoxolate)"),
    ("[O-]OCCO[O-]", "ethane-1,2-bis(peroxolate)"),
    ("[O-]OC1CCCC1", "cyclopentaneperoxolate"),
    ("[O-]Oc1ccc(O[O-])cc1", "benzene-1,4-bis(peroxolate)"),
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
    ("C/C(=N/[O-])c1ccccc1", "{[(1Z)-1-phenylethylidene]amino}oxidanide"),
    ("N/C(=N/[O-])c1ccc(NC(=O)O)cc1", "({(E)-amino[4-(carboxyamino)phenyl]methylidene}amino)oxidanide"),
    ("O[N-2]", "hydroxyazanediide"),
    ("C[SiH4-]", "methylsilanuide"),
    ("C[B-](C)(C)C", "tetramethylboranuide"),
    ("[C-]#[N+]C[B-](F)(F)F", "trifluoro(isocyanomethyl)boranuide"),
    ("O=[N+]([O-])O[B-](F)(F)F", "trifluoro(nitrooxy)boranuide"),
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
    ("[N-]=N", "diazenide"),
    ("CN=[N-]", "methyldiazenide"),
    ("N=N[NH-]", "triaz-2-en-1-ide"),
    ("CN=N[NH-]", "3-methyltriaz-2-en-1-ide"),
    ("C[N-]NC", "1,2-dimethylhydrazin-1-ide"),
    ("CN[NH-]", "2-methylhydrazin-1-ide"),
    ("C[SiH2][SiH2-]", "2-methyldisilan-1-ide"),
    ("C[Si](C)(C)[Si-](C)C", "pentamethyldisilan-1-ide"),
    ("C[As](C)(C)=[N-]", "trimethyl-λ5-arsaniminide"),
    ("[PH-]CC[AsH3-]", "(2-phosphanidylethyl)arsanuide"),
    ("[SiH2-]CC[PH-]", "(2-silanidylethyl)phosphanide"),
    ("[O-]CC([O-])C1CC[BH2-]CC1", "1-(borinan-1-uid-4-yl)ethane-1,2-bis(olate)"),
    ("[O-]CC([O-])[NH-]", "1-azanidylethane-1,2-bis(olate)"),
    ("[O-]CC([O-])[BH3-]", "1-boranuidylethane-1,2-bis(olate)"),
    ("[O-]P(=O)([O-])c1ccccc1", "phenylphosphonate"),
    ("OP(=O)(O)CP(=O)([O-])[O-]", "(phosphonomethyl)phosphonate"),
    ("OP(=O)(O)CCP(=O)([O-])[O-]", "(2-phosphonoethyl)phosphonate"),
    ("OP(=O)(O)C(O)P(=O)([O-])[O-]", "[hydroxy(phosphono)methyl]phosphonate"),
    ("OP(=O)(O)c1ccc(P(=O)([O-])[O-])cc1", "(4-phosphonophenyl)phosphonate"),
    ("OP(=O)(O)CP(=O)([O-])O", "hydrogen (phosphonomethyl)phosphonate"),
    ("OP(=O)([O-])CP(=O)(O)[O-]", "methylenebis(hydrogen phosphonate)"),
    ("OP(=O)([O-])CCP(=O)(O)[O-]", "(ethane-1,2-diyl)bis(hydrogen phosphonate)"),
    ("OP(=O)([O-])c1ccc(P(=O)(O)[O-])cc1", "(1,4-phenylene)bis(hydrogen phosphonate)"),
    ("[O-]P(=O)([O-])CP(=O)([O-])[O-]", "methylenebis(phosphonate)"),
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
        ("[O-]C(=O)c1ccc(cc1)[CH2-]", "(4-carboxylatophenyl)methanide"),
        ("[O-]C(=O)c1ccc(cc1)[C-](C#N)C#N", "(4-carboxylatophenyl)di(cyano)methanide"),
        ("[O-]C(=O)C[C-]1C=CC=C1", "1-(carboxylatomethyl)cyclopenta-2,4-dien-1-ide"),
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
    ("C1CC[SH3-]CC1", "1λ4-thian-1-uide"),
    ("[PH2-]12(CCCCC1)CCCCC2", "6λ5-phosphaspiro[5.5]undecan-6-uide"),
    ("[O-]P(=O)([O-])c1ccc(cc1)C(=O)[O-]", "4-phosphonatobenzoate"),
    ("[O-][As](=O)([O-])c1ccc(cc1)C(=O)[O-]", "4-arsonatobenzoate"),
    ("[O-]C(=O)CCP(=O)(O)[O-]", "3-[hydroxy(oxidanidyl)phosphoryl]propanoate"),
    ("[O-]CCP(=O)(O)[O-]", "hydrogen (2-oxidoethyl)phosphonate"),
    ("[O-]C(=O)c1ccc(cc1)CC[CH2-]", "3-(4-carboxylatophenyl)propan-1-ide"),
    ("[O-]C(=O)c1ccc(cc1)C[CH-]C", "1-(4-carboxylatophenyl)propan-2-ide"),
    ("[O-]C(=O)c1ccc(cc1)C=[CH-]", "2-(4-carboxylatophenyl)eth-1-en-1-ide"),
    ("[O-]C(=O)c1ccc(cc1)S[S-]", "4-disulfanidylbenzoate"),
    ("[O-]C(=O)c1ccc(cc1)SS[S-]", "4-trisulfanidylbenzoate"),
    ("[S-]Sc1ccc(cc1)[S-]", "4-disulfanidylbenzene-1-thiolate"),
    ("[O-]OCC[S-]", "2-sulfidoethane-1-peroxolate"),
    ("[O-]C(=O)c1ccc(cc1)C2=C(c3ccc(cc3)C(=O)[O-])[CH-]C=C2", "4,4'-(cyclopenta-1,4-dien-3-ide-1,2-diyl)dibenzoate"),
    ("[O-]C(=O)c1ccc(cc1)CCCc1ccc(cc1)C(=O)[O-]", "4,4'-(propane-1,3-diyl)dibenzoate"),
    ("OC(=O)CC1CC([CH-]C1)CC(=O)O", "2,4-bis(carboxymethyl)cyclopentan-1-ide"),
    ("OCC1CC([CH+]C1)CO", "2,4-bis(hydroxymethyl)cyclopentan-1-ylium"),
    ("OC(=O)CC[N-]CCC(=O)O", "2-carboxy-N-(2-carboxyethyl)ethan-1-aminide"),
    ("[O-]C(=O)CC=[B-]1Cc2ccccc2C1", "2-(2-carboxylatoethylidene)-1,3-dihydro-2H-2-benzoborol-2-uide"),
    ("[O-]C(=O)Cc1c[n-]cn1", "4-(carboxylatomethyl)-1H-imidazol-1-ide"),
    ("[O-]C(=O)Cc1c[n-]c(C(=O)[O-])n1", "4-(carboxylatomethyl)-1H-imidazol-1-ide-2-carboxylate"),
    ("N#C[C-](C#N)C1=C(C(=C(C#N)C#N)1)[C-](C#N)C#N", "[3-(dicyanomethylidene)cycloprop-1-ene-1,2-diyl]bis(dicyanomethanide)"),
    ("OC(=O)c1ccc(cc1)C1=C(c2ccc(cc2)C(=O)O)C1=C(C#N)C#N", "4,4'-[3-(dicyanomethylidene)cycloprop-1-ene-1,2-diyl]dibenzoic acid"),
    ("B1C=Cc2ccccc12", "1H-1-benzoborole"),
    ("C1=Cc2ccccc2[SiH2]1", "1H-1-benzosilole"),
    ("C1=Cc2ccccc2[PH]1", "1H-phosphindole"),
    ("C1C=C2C=CC=CC2=P1", "2H-phosphindole"),
    ("CO[B-]1(C)C=C(C)c2ccccc12", "1-methoxy-1,3-dimethyl-1H-1-benzoborol-1-uide"),
    ("C=[B-]1Cc2ccccc2C1", "2-methylidene-1,3-dihydro-2H-2-benzoborol-2-uide"),
    ("C=[B-]1C=CC=C1", "1-methylidene-1H-borol-1-uide"),
    ("C[B-]1(C)C=C2C=C[CH-]C2=C1", "2,2-dimethyl-2,4-dihydrocyclopenta[c]borol-4-id-2-uide"),
    ("CN1[C-2]C=CC=Cc2ccccc12", "1-methyl-1-benzazocine-2,2(1H)-diide"),
    ("C[N-]c1ccccc1", "N-methylbenzenaminide"),
    ("C[N-]c1ccc[cH-]1", "2-(methylazanidyl)cyclopenta-2,4-dien-1-ide"),
    ("c1cc[cH-]c1[CH2-]", "2-methanidylcyclopenta-2,4-dien-1-ide"),
    ("[CH2-]CC1=CC=C[CH-]1", "2-(ethan-1-id-2-yl)cyclopenta-2,4-dien-1-ide"),
    ("C[CH-]c1ccc[cH-]1", "2-(ethan-1-id-1-yl)cyclopenta-2,4-dien-1-ide"),
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


def test_charge_of_a_linking_aminide_is_not_dropped_beside_two_carboxylates():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]C(=O)CC[N-]CCC(=O)[O-]")


def test_ammonium_counter_ion_multiplied_with_bis_to_avoid_diazane_reading():
    assert smiles_to_iupac("[NH4+].[NH4+].[O-]S(=O)(=O)[O-]") == "bis(azanium) sulfate"


def test_unaccounted_fragment_is_never_dropped():
    with pytest.warns(NonPreferredNameWarning):
        assert smiles_to_iupac("CN.[O-]S(=O)(=O)O") == "methanamine—hydrogen sulfate (1/1)"

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
        pytest.param("OCc1ccc[n+](CO)c1CO", "1,2,3-tris(hydroxymethyl)pyridin-1-ium", id="cation_outranks_identical_alcohol_groups"),
        pytest.param("C[n+]1ccccc1[N+](=O)[O-]", "1-methyl-2-nitropyridin-1-ium", id="nitro_beside_the_cation"),
        pytest.param("C[n+]1cccc2ccccc12", "1-methylquinolin-1-ium", id="quinolinium"),
        pytest.param("CC[n+]1ccc2ccccc2c1", "2-ethylisoquinolin-2-ium", id="isoquinolinium"),
        pytest.param("c1ccc2[nH+]cccc2c1", "quinolin-1-ium", id="protonated_quinoline"),
        pytest.param("c1c[nH+]c[nH]1", "1H-imidazol-3-ium", id="indicated_hydrogen_before_the_cationic_centre"),
        pytest.param("C1=NCC[NH2+]1", "4,5-dihydro-1H-imidazol-1-ium", id="protonated_nitrogen_carries_the_indicated_hydrogen"),
        pytest.param("C1=[NH+]CNC1", "2,5-dihydro-1H-imidazol-3-ium", id="centre_numbered_before_the_hydro_prefixes"),
        pytest.param("c1ccc2[nH]c[nH+]c2c1", "1H-1,3-benzimidazol-3-ium", id="fused_azolium_keeps_the_indicated_hydrogen_locant"),
        pytest.param("C[n+]1ccc(-c2ncccn2)cc1", "1-methyl-4-(pyrimidin-2-yl)pyridin-1-ium", id="pyrimidine_on_the_cation"),
        pytest.param("C[n+]1ccc(-c2ccncn2)cc1", "1-methyl-4-(pyrimidin-4-yl)pyridin-1-ium", id="pyrimidin_4_yl_on_the_cation"),
        pytest.param("C[n+]1ccc(-c2cnccn2)cc1", "1-methyl-4-(pyrazin-2-yl)pyridin-1-ium", id="pyrazine_on_the_cation"),
        pytest.param("C[n+]1ccc(-c2ccnnc2)cc1", "1-methyl-4-(pyridazin-4-yl)pyridin-1-ium", id="pyridazine_on_the_cation"),
        pytest.param("C[n+]1ccccc1-c1ncccn1", "1-methyl-2-(pyrimidin-2-yl)pyridin-1-ium", id="diazine_beside_the_cationic_nitrogen"),
        pytest.param("C[n+]1cccc(-c2ncccn2)c1", "1-methyl-3-(pyrimidin-2-yl)pyridin-1-ium", id="diazine_at_position_3"),
        pytest.param("C[n+]1ccc(-c2nccs2)cc1", "1-methyl-4-(1,3-thiazol-2-yl)pyridin-1-ium", id="azole_on_the_cation"),
        pytest.param("C[n+]1ccc(-c2ccccn2)cc1", "1'-methyl-2,4'-bipyridin-1'-ium", id="identical_rings_are_an_assembly"),
        pytest.param("c1ccc(-c2ncccn2)nc1", "2-(pyridin-2-yl)pyrimidine", id="neutral_pair_keeps_the_senior_ring"),
        pytest.param("C[n+]1ccccc1.[Cl-]", "1-methylpyridin-1-ium chloride", id="pyridinium_salt"),
        pytest.param("[CH3+].[Cl-]", "methylium chloride", id="one_atom_cation_salt"),
        pytest.param(
            "O=C([O-])C[n+]1ccccc1CC(=O)[O-]", "2,2'-(pyridin-1-ium-1,2-diyl)diethanoate", id="ring_cation_linking_two_anions"
        ),
    ],
)
def test_heteroaromatic_ring_cation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-72.2.2.2.3, P-72.7, P-15.3.2, P-28.2.1: azolide anions joined by a linker or a bond are a multiplied parent, each unit named as the anion
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(C)(c1ccc[n-]1)c1ccc[n-]1", "2,2'-(propane-2,2-diyl)bis(1H-pyrrol-1-ide)", id="branched_linker"),
        pytest.param("C(c1ccc[n-]1)c1ccc[n-]1", "2,2'-methylenebis(1H-pyrrol-1-ide)", id="methylene_linker"),
        pytest.param("CCC(CC)(c1ccc[n-]1)c1ccc[n-]1", "2,2'-(pentane-3,3-diyl)bis(1H-pyrrol-1-ide)", id="longer_linker"),
        pytest.param("c1ccc(C(c2ccc[n-]2)c2ccc[n-]2)cc1", "2,2'-(phenylmethylene)bis(1H-pyrrol-1-ide)", id="ring_substituent_on_the_linker"),
        pytest.param("c1cc[n-]c1-c1ccc[n-]1", "2,2'-bi(1H-pyrrol-1-ide)", id="direct_bond_between_the_rings"),
        pytest.param("c1cn[n-]c1-c1cc[n-]n1", "3,5'-bi(1H-pyrazol-1-ide)", id="pyrazolides_with_different_junction_locants"),
        pytest.param("c1nc[n-]c1-c1ncc[n-]1", "2,5'-bi(1H-imidazol-1-ide)", id="imidazolides_with_different_junction_locants"),
        pytest.param("c1cc[n-]c1C(C)(C)c1ccc[nH]1", "2-[2-(1H-pyrrol-2-yl)propan-2-yl]-1H-pyrrol-1-ide", id="one_anionic_ring_beside_a_neutral_one"),
        pytest.param("c1cc[nH]c1C(C)(C)c1ccc[nH]1", "2,2'-(propane-2,2-diyl)bis(1H-pyrrole)", id="neutral_pair_for_comparison"),
    ],
)
def test_two_azolide_anions_joined_by_a_linker_or_a_bond(smiles, expected):
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
        pytest.param(
            "C(=[N+]1CCCCC1)[NH+]1CCCCC1",
            "1-[(piperidin-1-ium-1-yl)methylidene]piperidin-1-ium",
            id="ring_cation_ylidene_of_a_second_ring_cation",
        ),
        pytest.param("C[N+]1(C)CCCCC1", "1,1-dimethylpiperidin-1-ium", id="quaternary_ring_nitrogen"),
        pytest.param("C1COCC[NH2+]1", "morpholin-4-ium", id="cation_locant_from_the_ring_heteroatom_order"),
        pytest.param("C[NH+]1CCOCC1", "4-methylmorpholin-4-ium", id="n_alkylmorpholinium"),
        pytest.param("OC1CC[NH2+]CC1", "4-hydroxypiperidin-1-ium", id="hydroxy_on_a_cationic_ring"),
        pytest.param("C[N+]1(C)CCCC1C(=O)O", "2-carboxy-1,1-dimethylpyrrolidin-1-ium", id="cation_outranks_acid_on_a_saturated_ring"),
        pytest.param("c1ccc2oc[nH+]c2c1", "1,3-benzoxazol-3-ium", id="protonated_fused_azole"),
        pytest.param("[O-][N+]1(C)CCCCC1", "1-methylpiperidine 1-oxide", id="saturated_ring_n_oxide"),
        pytest.param("[O-][NH+]1CCOCC1", "morpholine 4-oxide", id="saturated_ring_nh_oxide"),
        pytest.param("CC(N)C[N+]1([O-])CCCC1", "1-(2-aminopropyl)pyrrolidine 1-oxide", id="ring_oxide_is_the_parent_beside_a_chain_amine"),
        pytest.param("[O-][n+]1ccoc1", "1,3-oxazole 3-oxide", id="azole_n_oxide_locant"),
        pytest.param("COc1ccc2no[n+]([O-])c2c1", "6-methoxy-2,1,3-benzoxadiazole 1-oxide", id="n_oxide_nitrogen_bonded_to_a_ring_oxygen_is_not_nitro"),
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
        pytest.param("Cc1cc(C)[s+]s1", "3,5-dimethyl-1,2-dithiol-1-ium", id="substituted_dithiolium"),
        pytest.param("C1=CC=[O+]S1", "1,2-oxathiol-1-ium", id="oxygen_cation_of_an_azole_type_ring"),
        pytest.param("C1=CC=[S+]N1", "1,2-thiazol-1-ium", id="cation_beside_a_ring_nitrogen_hydrogen"),
    ],
)
def test_chalcogen_cation_of_a_five_membered_mancude_ring(smiles, expected):
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
        pytest.param("[CH2]c1ccc([CH2-])cc1", "(4-methanidylphenyl)methyl", id="carbanide_on_a_ring_beside_a_radical"),
        pytest.param("[CH2]c1cc([CH2])cc([CH2-])c1", "[3-methanidyl-5-(ylomethyl)phenyl]methyl", id="carbanide_on_a_ring_beside_two_radicals"),
        pytest.param("[CH2]C1CC[CH+]CC1", "(cyclohexan-4-ylium-1-yl)methyl", id="ring_carbocation_beside_a_radical"),
        pytest.param("[CH]1[CH-]CCC1C1C[CH+]CC1", "[1,1'-bi(cyclopentan)]-3'-ylium-3-id-2-yl", id="ring_assembly_with_ionic_centres_beside_a_radical"),
        pytest.param("[CH2]C[C-]1CC[CH+]C1", "2-(cyclopentan-3-ylium-1-id-1-yl)ethyl", id="ring_ions_entered_through_an_ionic_atom_beside_a_radical"),
        pytest.param("[CH2]c1ccc([CH2-])c([CH2+])c1", "(4-methanidyl-3-methyliumylphenyl)methyl", id="carbanide_and_carbocation_beside_a_radical"),
        pytest.param("[CH2]C(C[CH-]C)C[CH+]C", "4-(propan-2-id-1-yl)pentan-2-ylium-5-yl", id="carbanide_arm_and_parent_carbocation_beside_a_radical"),
        pytest.param("[CH2]c1ccccc1[NH-]", "(2-azanidylphenyl)methyl", id="anionic_nitrogen_beside_a_radical"),
    ],
)
def test_ionic_centres_in_a_different_parent_than_the_radical(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-][n+]1cccnc1", "pyrimidine 1-oxide", id="diazine_n_oxide"),
        pytest.param("[O-][N+]1=NCCCc2ccccc21", "4,5-dihydro-3H-1,2-benzodiazepine 1-oxide", id="azoxy_ring_with_hydro_prefixes"),
        pytest.param("[O-][N+]1=NCCCC1", "3,4,5,6-tetrahydropyridazine 1-oxide", id="saturated_azoxy_ring"),
        pytest.param("Cc1cc[n+]([O-])cn1", "4-methylpyrimidine 1-oxide", id="n_oxide_locant_before_substituents"),
        pytest.param("Cc1ccnc[n+]1[O-]", "6-methylpyrimidine 1-oxide", id="oxide_nitrogen_numbered_lowest"),
        pytest.param("[O-][n+]1ccnc2ccccc12", "quinoxaline 1-oxide", id="fused_diazine_n_oxide"),
        pytest.param("Cc1c[n+]([O-])c2ccccc2[n+]1[O-]", "2-methylquinoxaline 1,4-dioxide", id="two_n_oxides_on_a_fused_diazine"),
        pytest.param("Nc1n[n+]([O-])c2ccccc2[n+]1[O-]", "3-amino-1,2,4-benzotriazine 1,4-dioxide", id="two_n_oxides_beside_an_unoxidized_nitrogen"),
        pytest.param("Cc1cc[n+]([O-])cc1-c1cccc[n+]1[O-]", "4'-methyl-2,3'-bipyridine 1,1'-dioxide", id="two_n_oxides_in_a_ring_assembly"),
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
        pytest.param("[c+]1cccc2ccccc12", "naphthalen-1-ylium", id="mancude_fused_ring_position_takes_no_added_hydrogen"),
        pytest.param("[CH+]1CC=Cc2ccccc12", "1,2-dihydronaphthalen-1-ylium", id="hydro_pair_on_fused_ring_stays_a_prefix"),
        pytest.param("C1CCC2CC[CH+]CC2C1", "decahydronaphthalen-2-ylium", id="saturated_fused_ring"),
        pytest.param("[CH+]1C=CC=CC=C1", "cyclohepta-2,4,6-trien-1-ylium", id="monocycle_drawn_aromatic_by_the_toolkit"),
        pytest.param("C[C+]1C=CC=CC=C1", "1-methylcyclohepta-2,4,6-trien-1-ylium", id="substituted_monocycle_drawn_aromatic"),
        pytest.param("[CH+]1C=CC=C2C1=CC=CC=C2", "benzo[7]annulenylium", id="cation_delocalized_over_a_fused_ring_system"),
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
        pytest.param("c1ccccc1[NH2+]", "benzenaminiumyl", id="aminiumyl_on_benzene_is_not_anilinium"),
        pytest.param("Cc1ccc(cc1)[NH2+]", "4-methylbenzenaminiumyl", id="substituted_benzenaminiumyl"),
        pytest.param("c1ccccc1[NH+]", "benzenaminylium", id="hydride_loss_cation_of_an_arylamine"),
        pytest.param("CC[N+](CC)(CC)[BH2-]", "(N,N-diethylethanaminiumyl)boranuidyl", id="zwitterionic_radical_on_boron"),
        pytest.param("C[Si](C)(C)[N-][N+]=NC", "3-methyl-1-(trimethylsilyl)triaz-2-en-2-ium-1-id-2-yl", id="zwitterionic_radical_on_the_cationic_atom"),
        pytest.param("[S-]1CC1", "1\u03bb4-thiiran-1-id-1-yl", id="ring_chalcogen_with_both_centres"),
        pytest.param("COC(=O)[C-]", "(methoxycarbonyl)methanidylidene", id="carbanidylidene_with_an_ester_prefix"),
        pytest.param("C[N-][NH]", "2-methylhydrazin-2-id-1-yl", id="radical_nitrogen_takes_the_lowest_locant_before_the_anion"),
        pytest.param("c1ccccc1[N-][N]", "2-phenylhydrazin-2-id-1-ylidene", id="two_radical_electrons_of_a_nitrogen_beside_an_anion"),
        pytest.param("CC[N+]1C(=O)C=CC=C1", "1-ethyl-2-oxopyridin-1-ium-1(2H)-yl", id="added_hydrogen_beside_the_ring_oxo"),
        pytest.param("[CH2]C#[O+]", "2-(oxidaniumylidyne)ethyl", id="cationic_prefix_on_a_radical_parent"),
        pytest.param("CC[n+]1c(sc2ccccc12)C=C[CH]C=Cc1sc2ccccc2[n+]1CC", "1,5-bis(3-ethyl-1,3-benzothiazol-3-ium-2-yl)penta-1,4-dien-3-yl", id="two_cationic_rings_on_a_radical_parent"),
        pytest.param("FC(F)(F)[C]1SSS[C+]1C(F)(F)F", "4,5-bis(trifluoromethyl)-1,2,3-trithiolan-5-ylium-4-yl", id="radical_cation_of_a_heterocycle"),
        pytest.param("C1=CC=C2C(=C1)[CH]c1ccccc1[CH-]2", "9,10-dihydroanthracen-10-id-9-yl", id="added_hydrogen_radical_anion"),
        pytest.param("c1ccc2c(c1)[CH+][CH]c1ccccc12", "9,10-dihydrophenanthren-10-ylium-9-yl", id="added_hydrogen_radical_cation"),
    ],
)
def test_radical_ions_named_through_the_filled_ion(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-72.6, P-72.7, P-73.6, P-73.7: anionic and cationic centres in parents and in prefixes
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param(
            "N#C[C-](C#N)C(CCC(=O)[O-])CC[C-](C#N)C#N",
            "2-(2-carboxylatoethyl)-1,1,5,5-tetracyanopentane-1,5-diide",
            id="more_ide_centres_outrank_an_idoate_parent",
        ),
        pytest.param("[O-]C(C1CC[BH2-]CC1)C[O-]", "1-(borinan-1-uid-4-yl)ethane-1,2-bis(olate)", id="more_anionic_centres_in_the_parent"),
        pytest.param("[PH-]CC[AsH3-]", "(2-phosphanidylethyl)arsanuide", id="uide_outranks_ide"),
        pytest.param("[NH3+]C(=O)CCC1CCC([NH3+])C1", "3-(3-azaniumylcyclopentyl)propanamidium", id="amidium_outranks_aminium"),
        pytest.param("NC(=[NH2+])CC(N)C(N)=[NH2+]", "2-aminobutanebis(imidamidium)", id="amidinium_units_with_an_amino_prefix"),
        pytest.param("[NH3+]CC(=[NH2+])CCN", "4-amino-2-azaniumylidenebutan-1-aminium", id="aminium_outranks_iminium_and_amino"),
        pytest.param("NC(=O)C(N)C(=[NH2+])N", "2,3-diamino-3-oxopropanimidamidium", id="amidinium_outranks_the_neutral_amide"),
        pytest.param("C[S+](C)CCCCCC[NH3+]", "6-(dimethylsulfaniumyl)hexan-1-aminium", id="sulfaniumyl_prefix_on_a_primary_aminium"),
        pytest.param("C[S+](C)CC[NH2+]C", "2-(dimethylsulfaniumyl)-N-methylethan-1-aminium", id="sulfaniumyl_prefix_on_a_secondary_aminium"),
        pytest.param("C[S+](C)CCCCCC[N+](C)(C)C", "6-(dimethylsulfaniumyl)-N,N,N-trimethylhexan-1-aminium", id="sulfaniumyl_prefix_on_a_quaternary_aminium"),
        pytest.param("C[S+](C)CCCCCC[P+](C)(C)C", "[6-(dimethylsulfaniumyl)hexyl]tri(methyl)phosphanium", id="phosphorus_outranks_sulfur"),
        pytest.param(
            "C[C+]1CC(C)=CC(=[O+]C)C1",
            "1,3-dimethyl-5-(methyloxidaniumylidene)cyclohex-3-en-1-ylium",
            id="ylium_outranks_ium",
        ),
        pytest.param("[N-2]", "azanediidyl", id="dianionic_radical_centre"),
        pytest.param("[C]1=[C][CH-]C=C1", "cyclopenta-1,4-dien-3-ide-1,2-diyl", id="ide_kept_before_a_diyl_suffix"),
        pytest.param("[CH-]([CH2])[CH2]", "propan-2-ide-1,3-diyl", id="ide_between_two_free_valences"),
        pytest.param("[O-][O]", "dioxidanidyl", id="radical_and_anion_on_a_peroxide"),
        pytest.param("[S-][S]", "disulfanidyl", id="radical_and_anion_on_a_disulfide"),
        pytest.param("[B-]1C=C2C=CC=CC2=C1", "2H-2-benzoborol-2-uid-2-ylidene", id="uide_and_ylidene_on_one_ring_atom"),
        pytest.param("[O-]C(=O)Cc1ccccc1C(=O)[O-]", "2-(carboxylatomethyl)benzoate", id="carboxylato_prefix"),
        pytest.param("[O-]S(=O)(=O)Cc1ccccc1C(=O)[O-]", "2-(sulfonatomethyl)benzoate", id="sulfonato_prefix"),
        pytest.param("[O-]C(=O)CC[Se](=O)[O-]", "3-seleninatopropanoate", id="seleninato_prefix"),
        pytest.param("[O-]c1cc2ccccc2cc1C(=O)[O-]", "3-oxidonaphthalene-2-carboxylate", id="oxido_prefix"),
        pytest.param("[O-]CC[NH3+]", "2-azaniumylethan-1-olate", id="ammonium_prefix_on_an_olate"),
        pytest.param("C[N+](C)(C)CC[O-]", "2-(N,N-dimethylmethanaminiumyl)ethan-1-olate", id="quaternary_cation_prefix_on_an_olate"),
        pytest.param("C[N+](C)(N)CC(=O)[O-]", "(1,1-dimethylhydrazin-1-ium-1-yl)acetate", id="hydrazinium_prefix_on_a_carboxylate"),
        pytest.param("C[NH2+]CC[O-]", "2-(methanaminiumyl)ethan-1-olate", id="secondary_cation_prefix_on_an_olate"),
        pytest.param("C[N+](C)(C)c1ccc([O-])cc1", "4-(N,N-dimethylmethanaminiumyl)phenoxide", id="cation_prefix_on_a_phenoxide"),
        pytest.param("[S-]CC[N+](C)(C)C", "2-(N,N-dimethylmethanaminiumyl)ethane-1-thiolate", id="ammonium_prefix_on_a_thiolate"),
        pytest.param("C[N+](C)(C)C1CC[CH-]CC1", "4-(N,N-dimethylmethanaminiumyl)cyclohexan-1-ide", id="ammonium_prefix_on_a_ring_carbanide"),
        pytest.param("C[N+](C)(C)C[PH3-]", "[(N,N-dimethylmethanaminiumyl)methyl]phosphanuide", id="ammonium_prefix_on_a_phosphanuide"),
        pytest.param("[O-]C(=O)C[NH2+]CC(=O)[O-]", "[(carboxylatomethyl)azaniumyl]acetate", id="cationic_centre_between_identical_carboxylates"),
        pytest.param("[O-]C(=O)C[NH+](CC(=O)[O-])CC(=O)[O-]", "[bis(carboxylatomethyl)azaniumyl]acetate", id="cationic_centre_between_three_identical_carboxylates"),
        pytest.param("[O-]CCC1=CC=C[CH-]1", "2-(2-oxidoethyl)cyclopenta-2,4-dien-1-ide", id="ide_ring_outranks_an_olate_chain"),
        pytest.param("[O-]CCc1cccc[c-]1", "2-(2-oxidoethyl)benzen-1-ide", id="benzenide_outranks_an_olate_chain"),
        pytest.param("[CH-2]c1ccc[cH-]1", "(cyclopenta-3,5-dien-2-id-1-yl)methanediide", id="ring_ide_prefix_on_a_diide_parent"),
        pytest.param("[CH-2]c1cccc[c-]1", "(benzen-2-id-1-yl)methanediide", id="benzene_ide_prefix_free_valence_first"),
        pytest.param("[O-]CC(C1CCC[CH-]1)C[O-]", "2-(cyclopentan-2-id-1-yl)propane-1,3-bis(olate)", id="cycloalkane_ide_prefix_on_bis_olate"),
        pytest.param("c1cc[cH-]c1[N-2]", "cyclopenta-2,4-dien-1-id-2-aminediide", id="ring_ide_beside_a_suffix_on_a_diene"),
        pytest.param("[O-]C(=O)C1=CC=C[CH-]1", "cyclopenta-2,4-dien-1-ide-2-carboxylate", id="ring_ide_before_a_consonant_suffix"),
    ],
)
def test_ionic_centres_in_parents_and_prefixes(smiles, expected):
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
        pytest.param("C[N+](C)=CC(=O)OC", "2-methoxy-N,N-dimethyl-2-oxoethan-1-iminium", id="cation_outranks_ester"),
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
        pytest.param("C[N+](C)(C)C1CCC([NH3+])CC1", "N1,N1,N1-trimethylcyclohexane-1,4-bis(aminium)", id="quaternary_and_primary_aminium_on_a_ring"),
        pytest.param("C[N+](C)(C)c1ccc([NH3+])cc1C", "N1,N1,N1,2-tetramethylbenzene-1,4-bis(aminium)", id="quaternary_and_primary_aminium_on_a_substituted_benzene"),
        pytest.param("[NH3+]CC(C(=O)O)C[NH3+]", "2-carboxypropane-1,3-bis(aminium)", id="bis_aminium_with_a_prefix_group"),
        pytest.param("C[N+](C)(C)c1ccc(cc1)[NH+](C)C", "N1,N1,N1,N4,N4-pentamethylbenzene-1,4-bis(aminium)", id="substituted_aminium_pair_on_a_benzene"),
        pytest.param("C[N+](C)(C)C1CCC(CC1)[N+](C)(C)C", "N1,N1,N1,N4,N4,N4-hexamethylcyclohexane-1,4-bis(aminium)", id="substituted_aminium_pair_on_a_cycloalkane"),
        pytest.param("C[N+](C)(C)CC(O)C[NH+](C)C", "2-hydroxy-N1,N1,N1,N3,N3-pentamethylpropane-1,3-bis(aminium)", id="substituted_aminium_pair_on_a_chain_with_a_hydroxy_group"),
        pytest.param("CC[N+](C)(C)c1cc(Cl)cc(c1)[N+](C)(C)C", "5-chloro-N1-ethyl-N1,N1,N3,N3,N3-pentamethylbenzene-1,3-bis(aminium)", id="prefixes_in_alphanumerical_order"),
        pytest.param("C[Si](C)(C)[NH3+]", "1,1,1-trimethylsilanaminium", id="aminium_on_a_silane_parent"),
        pytest.param(
            "C[Si](C)(C)[NH2+][SiH2][NH2+][Si](C)(C)C",
            "N,N'-bis(trimethylsilyl)silanebis(aminium)",
            id="bis_aminium_on_a_silane_parent_with_n_substituents",
        ),
        pytest.param("C1C[NH2+]CC[NH2+]1.[Cl-].[Cl-]", "piperazine-1,4-diium dichloride", id="dication_salt"),
    ],
)
def test_several_identical_cationic_centres(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-15.3.2, P-73.1.1.2, P-41: identical aminium or anilinium units joined by a divalent group are multiplied with 'bis', the cation
# outranking the ester or ether that joins them
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[NH3+]CCOCC[NH3+]", "2,2'-oxybis(ethan-1-aminium)", id="chain_units_joined_by_oxygen"),
        pytest.param("C[NH2+]CCOCC[NH2+]C", "2,2'-oxybis(N-methylethan-1-aminium)", id="secondary_chain_units"),
        pytest.param("C[N+](C)(C)CCSSCC[N+](C)(C)C", "2,2'-disulfanediylbis(N,N,N-trimethylethan-1-aminium)", id="quaternary_units_joined_by_a_dichalcogen"),
        pytest.param("[NH3+]CCOc1ccc(cc1)OCC[NH3+]", "2,2'-[1,4-phenylenebis(oxy)]bis(ethan-1-aminium)", id="concatenated_linker"),
        pytest.param("[NH3+]c1ccc(cc1)C(C)(C)c1ccc(cc1)[NH3+]", "4,4'-(propane-2,2-diyl)bis(anilinium)", id="anilinium_units_on_rings"),
        pytest.param("[NH3+]C1CCC(CC1)OC1CCC(CC1)[NH3+]", "4,4'-oxybis(cyclohexan-1-aminium)", id="cycloalkane_units"),
        pytest.param(
            "C[N+](C)(C)c1ccc(cc1)Cc1ccc(cc1)[N+](C)(C)C",
            "4,4'-methylenebis(N,N,N-trimethylanilinium)",
            id="quaternary_anilinium_units",
        ),
        pytest.param("C[NH2+]c1ccc(cc1)Oc1ccc(cc1)[NH2+]C", "4,4'-oxybis(N-methylanilinium)", id="secondary_anilinium_units"),
        pytest.param(
            "C[N+](C)(C)CCOC(=O)CCC(=O)OCC[N+](C)(C)C",
            "2,2'-[butanedioylbis(oxy)]bis(N,N,N-trimethylethan-1-aminium)",
            id="cation_outranks_the_ester_linker",
        ),
        pytest.param(
            "C[N+](C)(C)CCN(C)CC[N+](C)(C)C",
            "2,2'-(methylazanediyl)bis(N,N,N-trimethylethan-1-aminium)",
            id="neutral_amine_linker",
        ),
        pytest.param(
            "C[N+](C)(C)CC[N+](C)(C)CC[N+](C)(C)C",
            "N1,N1,N1,N2,N2-pentamethyl-N2-[2-(trimethylazaniumyl)ethyl]ethane-1,2-bis(aminium)",
            id="quaternary_ammonium_linker",
        ),
    ],
)
def test_aminium_units_joined_multiplicatively(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-73.1.1.2, P-31.2: a protonated ring nitrogen of a hydrogenated fused or bicyclic ring is the 'ium' of its parent hydride with the hydro prefixes kept
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[NH+]1CCc2ccccc2C1", "2-methyl-1,2,3,4-tetrahydroisoquinolin-2-ium", id="tertiary_nitrogen"),
        pytest.param("C1C[NH2+]Cc2ccccc21", "1,2,3,4-tetrahydroisoquinolin-2-ium", id="secondary_nitrogen_of_isoquinoline"),
        pytest.param("C1C[NH2+]c2ccccc21", "2,3-dihydro-1H-indol-1-ium", id="indoline"),
        pytest.param("C1CC[NH2+]c2ccccc21", "1,2,3,4-tetrahydroquinolin-1-ium", id="tetrahydroquinoline"),
        pytest.param("C1[NH2+]Cc2ccccc21", "2,3-dihydro-1H-isoindol-2-ium", id="isoindoline"),
        pytest.param("C1CC2CCCCC2C[NH2+]1", "decahydroisoquinolin-2-ium", id="fully_saturated_isoquinoline"),
        pytest.param("C1C[NH2+]CC2CCCCC21", "decahydroisoquinolin-2-ium", id="saturated_isoquinoline_other_atom_order"),
        pytest.param("C1C[NH+]2CCCCC2C1", "octahydroindolizin-4-ium", id="bridgehead_nitrogen_of_indolizine"),
        pytest.param("C1C[NH+]2CCCC2C1", "hexahydro-1H-pyrrolizin-4-ium", id="bridgehead_nitrogen_of_pyrrolizine"),
        pytest.param("C1C[NH+]2CCC1CC2.[Cl-]", "1-azabicyclo[2.2.2]octan-1-ium chloride", id="bridgehead_nitrogen_of_a_von_baeyer_system_with_its_salt"),
        pytest.param("[NH3+]CC[NH+]1CCc2ccccc2C1", "2-(1,2,3,4-tetrahydroisoquinolin-2-ium-2-yl)ethan-1-aminium", id="aminium_beside_a_hydro_fused_ring_cation"),
        pytest.param("[NH3+]CC[NH+]1CCC2CCCCC2C1", "2-(decahydroisoquinolin-2-ium-2-yl)ethan-1-aminium", id="aminium_beside_a_saturated_fused_ring_cation"),
        pytest.param("[NH3+]CC[NH+]1CCCCC1", "2-(piperidin-1-ium-1-yl)ethan-1-aminium", id="aminium_beside_a_monocyclic_ring_cation"),
        pytest.param("[NH3+]CC[NH+]1CC[NH2+]CC1", "1-(2-azaniumylethyl)piperazine-1,4-diium", id="ring_with_more_centres_is_the_parent"),
        pytest.param("C1C[NH2+]CC[NH+]1CC[NH3+]", "1-(2-azaniumylethyl)piperazine-1,4-diium", id="ring_with_more_centres_other_atom_order"),
    ],
)
def test_protonated_nitrogen_of_a_hydrogenated_fused_ring(smiles, expected):
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
        pytest.param("NC(=[NH2+])c1ccc[nH+]c1", "pyridin-1-ium-3-carboximidamidium", id="ring_ium_then_carboximidamidium"),
        pytest.param("[NH3+]c1ccc[n+]c1", "pyridin-1-ylium-3-aminium", id="ring_nitrenium_then_aminium"),
        pytest.param("C1=C(C=C[N+]2=CC=CC=C12)C(=[NH2+])N", "5λ5-quinolizin-5-ylium-2-carboximidamidium", id="lambda_ylium_ring_then_carboximidamidium"),
        pytest.param("c1cc([NH3+])c[n+]2ccccc12", "5λ5-quinolizin-5-ylium-3-aminium", id="lambda_ylium_ring_then_aminium"),
        pytest.param("C[n+]1ccc(N(C)C)cc1", "4-(dimethylamino)-1-methylpyridin-1-ium", id="neutral_substituted_amine_is_a_prefix_of_the_ring_cation"),
        pytest.param("CCNc1cc[n+](CC[NH3+])cc1", "2-[4-(ethylamino)pyridin-1-ium-1-yl]ethan-1-aminium", id="neutral_substituted_amine_beside_a_cationic_amine_is_a_prefix"),
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
        pytest.param("C[C@H]1CCC[N+]2(C1)CCCC2", "(7S)-7-methyl-5λ5-azaspiro[4.5]decan-5-ylium", id="spiro_ammonium_with_a_stereocentre"),
        pytest.param("[P+]12(c3ccccc3-c3ccccc31)c1ccccc1-c1ccccc12", "5H-5λ5,5′-spirobi[benzo[b]phosphindol]-5-ylium", id="spiro_atom_cation_of_mancude_components_keeps_its_indicated_hydrogen"),
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
        pytest.param("CC[NH3+].[C-]#N", "ethanaminium cyanide", id="cyanide_with_an_organic_cation"),
        pytest.param("[Ca+2].[C-]#N.[C-]#N", "calcium dicyanide", id="cyanide_with_a_multiplied_anion"),
        pytest.param("C[N+](C)(C)C.[N-]=[N+]=[N-]", "N,N,N-trimethylmethanaminium azide", id="azide_with_an_organic_cation"),
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
    "smiles,expected",
    [
        pytest.param("C[n+]1ccccc1CC[n+]1ccccc1", "1-[2-(1-methylpyridin-1-ium-2-yl)ethyl]pyridin-1-ium", id="units_differ"),
        pytest.param(
            "C[n+]1ccc(OCc2cc[n+](C)cc2)cc1",
            "1-methyl-4-[(1-methylpyridin-1-ium-4-yl)methoxy]pyridin-1-ium",
            id="unsymmetrical_linking_group",
        ),
        pytest.param(
            "C[n+]1ccc(CC[P+](C)(C)C)cc1", "1-methyl-4-[2-(trimethylphosphaniumyl)ethyl]pyridin-1-ium", id="ring_cation_senior_to_phosphonium"
        ),
        pytest.param(
            "C[n+]1ccc(/C=C/c2cc[n+](C)cc2)cc1", "4,4'-[(1E)-ethene-1,2-diyl]bis(1-methylpyridin-1-ium)", id="stereo_double_bond_of_the_linking_group"
        ),
    ],
)
def test_cations_that_are_not_a_multiplicative_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
            "C[n+]1[nH]n(C)c(=O)c1C(=O)[O-]",
            "1,3-dimethyl-5-oxo-2,5-dihydro-1H-1,2,3-triazol-3-ium-4-carboxylate",
            id="ring_cation_with_oxo_group_beside_the_anionic_suffix",
        ),
        pytest.param(
            "C[n+]1ccn(C)c(=O)c1CC(=O)[O-]",
            "(1,4-dimethyl-3-oxo-3,4-dihydropyrazin-1-ium-2-yl)acetate",
            id="ring_cation_with_oxo_group_beside_an_anionic_prefix",
        ),
        pytest.param(
            "O=C([O-])CCn1cc[n+](C)c1",
            "3-(3-methyl-1H-imidazol-3-ium-1-yl)propanoate",
            id="ring_cation_as_prefix_with_second_nitrogen",
        ),
        pytest.param(
            "O=C([O-])CCCC[NH+]1CCS(=O)CC1",
            "5-(1-oxo-1\u03bb4,4-thiazinan-4-ium-4-yl)pentanoate",
            id="ring_cation_prefix_with_a_sulfoxide_sulfur",
        ),
        pytest.param(
            "O=C([O-])CCC[NH+]1CCS(=O)(=O)CC1",
            "4-(1,1-dioxo-1\u03bb6,4-thiazinan-4-ium-4-yl)butanoate",
            id="ring_cation_prefix_with_a_sulfone_sulfur",
        ),
        pytest.param(
            "[O-]S(=O)(=O)CC[NH+]1CCS(=O)CC1",
            "2-(1-oxo-1\u03bb4,4-thiazinan-4-ium-4-yl)ethane-1-sulfonate",
            id="ring_cation_prefix_with_a_sulfoxide_beside_a_sulfonate",
        ),
        pytest.param(
            "O=C([O-])CCC[NH+]1CCCS(=O)C1",
            "4-(1-oxo-1\u03bb4,3-thiazinan-3-ium-3-yl)butanoate",
            id="sulfoxide_sulfur_and_cationic_nitrogen_in_a_1_3_relation",
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
        pytest.param("CN(C)[N]", "2,2-dimethylhydrazin-1-ylidene", id="isodiazene_as_the_divalent_nitrogen"),
        pytest.param("C[N+](C)=[N-]", "dimethylhydrazinylidene", id="isodiazene_as_the_zwitterion"),
        pytest.param("CC[N+](C)=[N-]", "ethyl(methyl)hydrazinylidene", id="isodiazene_prefixes_alphabetized"),
        pytest.param("[NH2][NH]", "hydrazinyl", id="retained_hydrazine_root"),
        pytest.param("[NH][NH]", "hydrazine-1,2-diyl", id="two_radical_centres_on_a_chain"),
        pytest.param("[BH][BH]", "diborane(4)-1,2-diyl", id="two_radical_centres_on_a_polyborane"),
        pytest.param("[B][BH]", "diboran(4)-1-yl-2-ylidene", id="polyborane_mixed_valences_elide_before_the_first_suffix"),
        pytest.param("CC[B][B]CC", "1,2-diethyldiborane(4)-1,2-diyl", id="substituted_polyborane_radical_centres"),
        pytest.param("[SiH2][SiH3]", "disilanyl", id="diatomic_chain_cites_no_locant"),
        pytest.param("P[P]", "diphosphanylidene", id="diatomic_chain_ylidene"),
        pytest.param("CC(C)(C)O[P](c1ccccc1)(c1ccccc1)c1ccccc1", "tert-butoxytri(phenyl)-λ5-phosphanyl", id="lambda5_phosphanyl"),
        pytest.param("C[P](C)=O", "dimethylphosphinoyl", id="phosphinoyl_acyl_radical"),
        pytest.param("c1ccccc1[S](=O)=O", "benzenesulfonyl", id="sulfonyl_radical"),
        pytest.param("[C](=O)O", "hydroxycarbonyl", id="carbonic_acid_acyl_radical"),
        pytest.param("C[P](C)(C)=[N]", "P,P,P-trimethyl-λ5-phosphaniminyl", id="phosphaniminyl"),
        pytest.param("c1ccccc1[N]", "benzenaminylidene", id="nitrene_on_a_ring"),
        pytest.param("CC(=O)[N]", "acetamidylidene", id="acylnitrene"),
        pytest.param("CC(C)(C)S[S]", "tert-butyldisulfanyl", id="disulfanyl_radical"),
        pytest.param("N[O]", "aminoxyl", id="aminoxyl"),
        pytest.param("ClCN(CCl)[O]", "bis(chloromethyl)aminoxyl", id="substituted_aminoxyl"),
        pytest.param("CS(=O)(=O)[O]", "(methanesulfonyl)oxyl", id="oxyl_on_a_sulfonyl_group"),
        pytest.param("CCC[B]C#N", "cyano(propyl)boranyl", id="radical_outranks_nitrile"),
        pytest.param("CCC[B]O", "hydroxy(propyl)boranyl", id="radical_outranks_borinic_acid"),
        pytest.param("O[B]O", "dihydroxyboranyl", id="radical_centre_is_not_a_borono_group"),
        pytest.param("C[SiH]NO", "(hydroxyamino)(methyl)silyl", id="nitrogen_oxygen_chain_on_the_centre"),
        pytest.param("C[SiH]OS", "methyl(sulfanyloxy)silyl", id="oxygen_sulfur_chain_on_the_centre"),
        pytest.param("C[SiH]SO", "(hydroxysulfanyl)(methyl)silyl", id="sulfur_oxygen_chain_on_the_centre"),
        pytest.param("O[Si]C(O)O", "(dihydroxymethyl)(hydroxy)silylidene", id="heteroatom_groups_on_divalent_centre"),
        pytest.param("C[SiH2][SiH]", "2-methyldisilanylidene", id="divalent_centre_on_a_substituted_chain"),
        pytest.param("C1CC[SiH]CC1", "silinan-1-yl", id="ring_heteroatom_centre"),
        pytest.param("C[Si](C)(C)O[O]", "(trimethylsilyl)peroxyl", id="peroxyl_on_a_heteroatom_group"),
        pytest.param("C[Sn](C)[Sn](C)C", "1,1,2,2-tetramethyldistannane-1,2-diyl", id="two_centres_on_a_substituted_chain"),
        pytest.param("C[Si](C)[Si]C", "1,1,2-trimethyldisilan-1-yl-2-ylidene", id="mixed_centres_on_a_substituted_chain"),
        pytest.param("[SiH2]C[SiH2]", "methylenebis(silyl)", id="centres_on_separate_parents"),
        pytest.param("[SiH]1CC[SiH]C1", "1,3-disilolane-1,3-diyl", id="centres_on_skeletal_atoms_of_a_replacement_ring"),
        pytest.param("[CH3][Al][O][Al][CH3]", "1,3-dimethyldialuminoxane-1,3-diyl", id="two_centres_on_an_alternating_chain"),
        pytest.param("C[Al](C)S[Al]C", "1,3,3-trimethyldialuminathianyl", id="one_centre_on_a_chain_bridged_by_sulfur"),
        pytest.param("[CH3][Tl][CH3]", "dimethylthallanyl", id="thallium_with_two_bonds_has_a_radical_centre"),
        pytest.param("Cl[Tl]", "chlorothallanylidene", id="thallium_with_one_bond_has_two_radical_centres"),
        pytest.param("C[SiH]CC[SiH]C", "(ethane-1,2-diyl)bis(methylsilyl)", id="identical_substituted_hydride_radicals"),
    ],
)
def test_radicals_on_heteroatom_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-71.2.3, P-71.2.4, P-71.3.1, P-71.3.3, P-71.4, P-71.6: several radical centres
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C]C[C]C", "pentane-2,4-diylidene", id="two_ylidene_centres"),
        pytest.param("[CH2][CH]", "ethan-1-yl-2-ylidene", id="yl_before_ylidene"),
        pytest.param("[CH2][C]", "ethan-1-yl-2-ylidyne", id="yl_before_ylidyne"),
        pytest.param("[CH2]C[CH]", "propan-1-yl-3-ylidene", id="mixed_centres_far_apart"),
        pytest.param("[CH2]C1CC1[CH2]", "(cyclopropane-1,2-diyl)dimethyl", id="assembly_of_methyl_groups"),
        pytest.param("[CH2]c1ccc(cc1)[CH2]", "(benzene-1,4-diyl)dimethyl", id="assembly_on_benzene"),
        pytest.param("[CH2]C1CCC2(CC1)CCCC2[CH2]", "(spiro[4.5]decane-1,8-diyl)dimethyl", id="assembly_on_spiro_skeleton"),
        pytest.param("[O]C(C)(C)CC(C)(C)[O]", "(2,4-dimethylpentane-2,4-diyl)bis(oxyl)", id="assembly_of_oxyl_groups"),
        pytest.param("[NH]CC[NH]", "(ethane-1,2-diyl)bis(aminyl)", id="assembly_of_aminyl_groups"),
        pytest.param("C[N]CC[N]C", "N1,N2-dimethylethane-1,2-bis(aminyl)", id="substituted_aminyl_groups_on_one_parent"),
        pytest.param("C[N]C[N]C", "N,N'-dimethylmethanebis(aminyl)", id="substituted_aminyl_groups_without_a_multiplicative_amine_name"),
        pytest.param("[O][O]C1CC(C1)O[O]", "(cyclobutane-1,3-diyl)bis(peroxyl)", id="assembly_of_peroxyl_groups"),
        pytest.param("[S]Sc1ccc2cc(S[S])ccc2c1", "(naphthalene-2,6-diyl)bis(disulfanyl)", id="assembly_of_disulfanyl_groups_on_a_fused_ring"),
        pytest.param("O=[S]c1ccc(cc1)[S]=O", "benzene-1,4-disulfinyl", id="two_sulfinyl_radicals"),
        pytest.param("O=[C]c1ccc(cc1)[C]=O", "benzene-1,4-dicarbonyl", id="two_carbonyl_radicals"),
        pytest.param("O=[C]CC[C]=O", "butanedioyl", id="two_acyl_radicals_on_a_chain"),
        pytest.param("[N]=C=[N]", "methanebis(iminyl)", id="two_iminyl_radicals"),
        pytest.param("O=C([NH])c1ccccc1C(=O)[NH]", "benzene-1,2-bis(carboxamidyl)", id="two_carboxamidyl_radicals"),
        pytest.param("O=C([NH])CCC(=O)[NH]", "butanebis(amidyl)", id="two_amidyl_radicals"),
        pytest.param("O=C([N])CCCCC(=O)[N]", "hexanebis(amidylidene)", id="two_amidylidene_radicals"),
        pytest.param("[CH2][CH]C1CC[CH]CC1", "1-(4-ylocyclohexyl)ethane-1,2-diyl", id="parent_with_most_centres"),
        pytest.param("[CH2]Cc1cccc(c1)[C]C", "2-[3-(1,1-diyloethyl)phenyl]ethyl", id="parent_with_most_yl_centres"),
        pytest.param("C[CH]C1CC[CH]C1", "3-(1-yloethyl)cyclopentyl", id="ring_outranks_chain"),
        pytest.param("C[C]C(=NC)C", "3-(methylimino)butan-2-ylidene", id="imidoyl_carbene"),
        pytest.param("CC(=NC)[N]", "(N-methylethanimidoyl)azanylidene", id="imidoyl_nitrene"),
        pytest.param("CCC(=O)[C]C", "3-oxopentan-2-ylidene", id="acyl_carbene"),
        pytest.param("C=C[CH]", "prop-2-en-1-ylidene", id="vinyl_carbene"),
        pytest.param("[C]1C[C]=Cc2ccccc12", "naphthalen-3-yl-1(2H)-ylidene", id="yl_and_ylidene_with_added_hydrogen"),
        pytest.param("[CH]1CCC(C)C[CH]1", "4-methylcyclohexane-1,2-diyl", id="ring_diyl_with_an_alkyl_prefix"),
        pytest.param("[CH]1CC(Cl)C(O)C[CH]1", "4-chloro-5-hydroxycyclohexane-1,2-diyl", id="ring_diyl_with_a_hydroxy_prefix"),
    ],
)
def test_radicals_with_several_centres(smiles, expected):
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


# P-62.5, P-73.1.1.2, P-73.1.2.1, P-73.5, P-74.1.3: ionic centres on heteroatoms of carbon skeletons
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-][NH+]1CCCO1", "1,2-oxazolidine 2-oxide", id="n_oxide_of_a_secondary_ring_nitrogen"),
        pytest.param("C[N+](C)=C=N", "(iminomethylidene)di(methyl)azanium", id="imino_group_on_an_ylidene"),
        pytest.param("C#[N+]C=C=C", "methylidyne(propa-1,2-dien-1-yl)azanium", id="formonitrile_is_not_substituted"),
        pytest.param("CC#[N+]C", "N-methylacetonitrilium", id="n_substituted_nitrilium"),
        pytest.param("OCCC#[N+]C", "3-hydroxy-N-methylpropanenitrilium", id="n_prefix_cited_among_the_chain_prefixes"),
        pytest.param("Cc1ccc(C#[N+]C)cc1", "N,4-dimethylbenzonitrilium", id="n_prefix_merged_with_a_ring_prefix"),
        pytest.param("N#CCC#[NH+]", "cyanoacetonitrilium", id="one_cationic_group_of_two_nitriles"),
        pytest.param("C=[P+](O)OC", "hydroxy(methoxy)(methylidene)phosphanium", id="phosphanium_with_ylidene_and_heteroatom_groups"),
        pytest.param("CCSC#[O+]", "[(ethylsulfanyl)methylidyne]oxidanium", id="oxidanium_with_ylidyne"),
        pytest.param("C[N+](C)=CC=[N+](C)C", "N1,N1,N2,N2-tetramethylethane-1,2-bis(iminium)", id="two_iminium_centres"),
        pytest.param("C[N+](C)=CC=N", "(2-iminoethylidene)di(methyl)azanium", id="one_iminium_centre_of_two_imines"),
        pytest.param("C[S+](C)CC[O-]", "2-(dimethylsulfaniumyl)ethan-1-olate", id="sulfonium_beside_an_anionic_group"),
        pytest.param("C[P+](C)(C)CC[O-]", "2-(trimethylphosphaniumyl)ethan-1-olate", id="phosphonium_beside_an_anionic_group"),
    ],
)
def test_ionic_centres_on_skeletal_heteroatoms(smiles, expected):
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
        ("[NH3+]CC1CCC([NH3+])CO1", "6-(azaniumylmethyl)oxan-3-aminium"),
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
        pytest.param("C/C=C\\[Si+](C)C", "dimethyl[(1Z)-prop-1-en-1-yl]silylium", id="silylium-alkenyl-descriptor"),
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
        pytest.param("[AlH2+]", "alumanylium", id="alumanylium"),
        pytest.param("P[PH+]", "diphosphanylium", id="diphosphanylium"),
        pytest.param("[GeH3][GeH2+]", "digermanylium", id="digermanylium"),
        pytest.param("CP[PH+]", "2-methyldiphosphan-1-ylium", id="substituted_dinuclear_chain_cites_locants"),
        pytest.param("C[P+]P", "1-methyldiphosphan-1-ylium", id="substituent_on_the_cationic_atom"),
        pytest.param("C[Al+]Cl", "chloro(methyl)alumanylium", id="substituted_alumanylium"),
        pytest.param("C[Si+2]C", "dimethylsilanebis(ylium)", id="silane_bis_ylium"),
        pytest.param("C[Si+3]", "methylsilanetris(ylium)", id="silane_tris_ylium"),
        pytest.param("[AlH+2]", "alumanebis(ylium)", id="alumane_bis_ylium"),
        pytest.param("[Cl+]", "chloranylium", id="chloranylium"),
        pytest.param("C(C)=[OH+]", "ethylideneoxidanium", id="ylidene_oxonium"),
        pytest.param("C(C)[O+]=C(C)C", "ethyl(propan-2-ylidene)oxidanium", id="mixed_ylidene_and_alkyl_oxonium"),
        pytest.param("OC(C)=[OH+]", "(1-hydroxyethylidene)oxidanium", id="hydroxyethylidene_oxonium"),
        pytest.param("OC(C)=[N+](C)C", "(1-hydroxyethylidene)di(methyl)azanium", id="ylidene_azanium"),
        pytest.param("C(C)(=O)[Cl+]C", "acetyl(methyl)chloranium", id="chloronium_with_an_acyl_group"),
        pytest.param("CCC[N+](C)(C)O", "N-hydroxy-N,N-dimethylpropan-1-aminium", id="aminium_with_a_hydroxy_group_on_the_cationic_nitrogen"),
        pytest.param("CC[N+](C)(C)OC", "N-methoxy-N,N-dimethylethanaminium", id="aminium_with_an_alkoxy_group_on_the_cationic_nitrogen"),
        pytest.param("[NH3+]OCC", "ethoxyazanium", id="azanium_with_only_oxygen_substituents"),
        pytest.param("C[S+](C)OCC", "ethoxydi(methyl)sulfanium", id="sulfanium_with_an_alkoxy_group"),
        pytest.param("CNN=[N+]", "3-methyltriaz-1-en-1-ylium", id="triazenylium"),
        pytest.param("C[Si](C)(C)[Si+](C)[Si](C)(C)C", "heptamethyltrisilan-2-ylium", id="trisilanylium"),
        pytest.param("CN(C)[N+2]", "2,2-dimethylhydrazine-1,1-bis(ylium)", id="hydrazine_bis_ylium"),
        pytest.param("O=C(c1ccccc1)N[NH+]", "2-benzoylhydrazin-1-ylium", id="acylhydrazinylium"),
        pytest.param("CC(C)(C)O[O+]", "tert-butylperoxylium", id="peroxylium"),
        pytest.param("c1ccccc1S[S+]", "phenyldisulfanylium", id="disulfanylium"),
        pytest.param("CCS[O+]", "(ethylsulfanyl)oxylium", id="sulfanyl_oxylium"),
        pytest.param("CN(C)[O+]", "N-methylmethanaminoxylium", id="aminoxylium"),
        pytest.param("CCO[Te+]", "ethoxytellanylium", id="alkoxy_tellanylium"),
        pytest.param("CC[N+2]", "ethylazanebis(ylium)", id="alkylazane_bis_ylium"),
        pytest.param("O=C(c1ccccc1)[N+2]", "benzoylazanebis(ylium)", id="acylazane_bis_ylium"),
        pytest.param("C=C[S+]=O", "ethenesulfinylium", id="sulfinylium"),
        pytest.param("C[P+](C)=O", "dimethylphosphinoylium", id="phosphinoylium"),
        pytest.param("C[P+](=O)CP", "methyl(phosphanylmethyl)phosphinoylium", id="phosphinoylium_beside_a_phosphanyl_group"),
        pytest.param("C[P+](=O)PC", "methyl(methylphosphanyl)phosphinoylium", id="phosphinoylium_bonded_to_phosphorus"),
        pytest.param("C[P+2]=O", "methylphosphonobis(ylium)", id="phosphono_bis_ylium"),
        pytest.param("c1cc[c+]o1", "furan-2-ylium", id="furan_ylium"),
        pytest.param("O=C(c1ccc[nH]1)[NH+]", "1H-pyrrole-2-carboxamidylium", id="heteroaryl_carboxamidylium"),
        pytest.param("CC(=O)C(C(C)=O)[N+]#N", "2,4-dioxopentane-3-diazonium", id="diazonium_beside_ketones"),
        pytest.param("N#[N+]c1ccc(cc1)[N+]#N", "benzene-1,4-bis(diazonium)", id="bis_diazonium"),
        pytest.param("C[C@H](CC)[N+]#N.[Cl-]", "(2R)-butane-2-diazonium chloride", id="diazonium_salt_cites_the_configuration"),
        pytest.param("C/C=C/C[N+]#N", "(2E)-but-2-ene-1-diazonium", id="diazonium_cites_the_double_bond"),
        pytest.param("N#[N+]c1cccc(c1)C[N+]#N", "3-[(diazyn-1-ium-1-yl)methyl]benzenediazonium", id="diazonium_on_a_side_chain_of_a_ring_diazonium"),
        pytest.param("N#[N+]CS", "1-sulfanylmethanediazonium", id="cation_outranks_a_sulfanyl_on_the_same_carbon"),
        pytest.param("N#[N+]c1ccc(cc1)Oc1ccc(cc1)[N+]#N", "4,4'-oxybis(benzenediazonium)", id="diazonium_rings_joined_by_a_linker"),
        pytest.param("N#[N+]c1ccc(cc1)C(=O)c1ccc(cc1)[N+]#N", "4,4'-carbonylbis(benzenediazonium)", id="cations_outrank_a_ketone_linker"),
        pytest.param("[NH3+]c1ccc(cc1)C(=O)c1ccc(cc1)[NH3+]", "4,4'-carbonylbis(anilinium)", id="aminium_rings_joined_by_a_carbonyl"),
        pytest.param("N#[N+]c1ccc(cc1)C(=O)OCC[N+]#N", "4-{[2-(diazyn-1-ium-1-yl)ethoxy]carbonyl}benzenediazonium", id="diazonium_in_the_alkoxy_of_an_ester"),
        pytest.param("N#[N+]c1ccc(cc1)C(=O)Oc1ccc(cc1)[N+]#N", "4-{[4-(diazyn-1-ium-1-yl)benzoyl]oxy}benzenediazonium", id="diazonium_in_the_acyl_of_an_ester"),
        pytest.param("N#[N+]c1ccc(cc1)C(=O)c1ccc(cc1)C(=O)O", "4-(4-carboxybenzoyl)benzenediazonium", id="carboxy_beside_the_acyl_group_of_a_cation"),
        pytest.param("N#[N+]CC(=O)OC(=O)CC[N+]#N", "3-{[(diazyn-1-ium-1-yl)acetyl]oxy}-3-oxopropane-1-diazonium", id="diazonium_in_the_acyl_of_an_anhydride"),
        pytest.param("N#[N+]c1ccc(cc1)OP(=O)(O)Oc1ccc(cc1)[N+]#N", "4,4'-[(hydroxyphosphoryl)bis(oxy)]bis(benzenediazonium)", id="cations_outrank_a_phosphate_linker"),
        pytest.param("N#[N+]CCC(=O)NCCNC(=O)CC[N+]#N", "3,3'-[ethane-1,2-diylbis(azanediyl)]bis(3-oxopropane-1-diazonium)", id="cations_outrank_a_diamide_linker"),
        pytest.param("N#[N+]CCOC(=O)P(=O)(O)O", "2-[(phosphonocarbonyl)oxy]ethane-1-diazonium", id="carbonyl_on_a_phosphorus_acid_group_of_a_cation"),
        pytest.param("N#[N+]c1ccc(cc1)C(=O)c1ccc(cc1)C(=O)OC", "4-[4-(methoxycarbonyl)benzoyl]benzenediazonium", id="ester_beside_the_acyl_group_of_a_cation"),
        pytest.param("N#[N+]c1ccc(cc1)Nc1ccc(cc1)[N+]#N", "4,4'-azanediylbis(benzenediazonium)", id="neutral_amine_linking_two_diazonium_rings"),
        pytest.param("[NH3+]CCNCC[NH3+]", "2,2'-azanediylbis(ethan-1-aminium)", id="neutral_amine_linking_two_aminium_chains"),
        pytest.param("N#[N+]c1ccc(cc1)C(C(=O)S)c1ccc(cc1)[N+]#N", "4,4'-(2-oxo-2-sulfanylethane-1,1-diyl)bis(benzenediazonium)", id="thioic_acid_group_on_the_methylene_between_cations"),
        pytest.param("N#[N+]CCOP(=O)(O)CC[N+]#N", "2-{[2-(diazyn-1-ium-1-yl)ethoxy](hydroxy)phosphoryl}ethane-1-diazonium", id="unsymmetrical_phosphinate_linker_between_cations"),
        pytest.param("O=C(c1ccccc1)N[N+](C)(C)C", "N',N',N'-trimethylbenzohydrazid-N'-ium", id="hydrazide_ium"),
    ],
)
def test_ylium_and_onium_cations_of_mononuclear_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diazonium_prefix_does_not_leak_into_an_ammonium_name():
    smiles_to_iupac("N#[N+]CCS(=O)(=O)NCCC[N+]#N")
    assert smiles_to_iupac("[NH3+]CCS(=O)(=O)NCCC[NH3+]") == "3-(2-azaniumylethane-1-sulfonamido)propan-1-aminium"


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
        ("CN([N+](C)(C)C)C", "pentamethylhydrazinium"),
        ("C[N+](C)(C)N", "1,1,1-trimethylhydrazin-1-ium"),
        ("CC/C=N/[N+](C)(C)C", "1,1,1-trimethyl-2-[(1E)-propylidene]hydrazin-1-ium"),
        ("[OH2+]C(C)=[OH+]", "(1-oxidaniumylethylidene)oxidanium"),
        ("C[N+](C)=N", "1,1-dimethyldiazen-1-ium"),
        ("CC(=O)N=[NH2+]", "2-acetyldiazen-1-ium"),
        ("C[N+](C)(C)P", "N,N-dimethyl-N-phosphanylmethanaminium"),
        ("C[S+](C)N", "aminodi(methyl)sulfanium"),
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[K+].C[As](C)(=O)[O-]", "potassium dimethylarsinate"),
        ("[Na+].OB(OC#N)[O-]", "sodium hydrogen borocyanatidate"),
    ],
)
def test_salts_of_arsenic_and_cyanato_boron_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)OCCC(=O)[O-].[Na+]", "sodium 3-(acetyloxy)propanoate", id="carboxylate_outranks_ester"),
        pytest.param(
            "CC(=O)OCCC(=O)[O-].CC(=O)OCCC(=O)[O-].[Ca+2]",
            "calcium bis[3-(acetyloxy)propanoate]",
            id="multiplied_anion_with_inner_parentheses_takes_square_brackets",
        ),
    ],
)
def test_carboxylate_salts_of_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NC(=[NH2+])c1ccccc1", "benzenecarboximidamidium", id="aryl_amidinium"),
        pytest.param("CC(N)=[NH2+].[Cl-]", "ethanimidamidium chloride", id="chain_amidinium_salt"),
        pytest.param("CNCc1cncc(C(N)=[NH2+])c1", "5-[(methylamino)methyl]pyridine-3-carboximidamidium", id="amidinium_with_a_neutral_amine"),
        pytest.param("CN(C)C(=[NH2+])c1ccccc1", "N,N-dimethylbenzenecarboximidamidium", id="n_substituted_amidinium"),
    ],
)
def test_amidinium_cations_take_the_imidamidium_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-]C1=NCCCC1", "3,4,5,6-tetrahydropyridin-2-olate", id="hydro_ring_olate"),
        pytest.param(
            "CCC1(CC)C(=O)N=C([O-])NC1=O.[Na+]",
            "sodium 5,5-diethyl-4,6-dioxo-1,4,5,6-tetrahydropyrimidin-2-olate",
            id="lactim_anion_of_a_barbiturate",
        ),
    ],
)
def test_lactim_and_hydro_ring_anions_take_the_olate_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_a_phenolate_ester_never_loses_its_charge():
    assert smiles_to_iupac("CC(=O)Oc1ccc([O-])cc1") == "4-(acetyloxy)phenoxide"


# P-73.6, P-74.1.3: the anion is the parent of a zwitterion and cites the diazonium group as a prefix
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-]C(=O)c1ccc(cc1)[N+]#N", "4-(diazyn-1-ium-1-yl)benzoate", id="carboxylate_of_a_ring"),
        pytest.param("[O-]C(=O)CC[N+]#N", "3-(diazyn-1-ium-1-yl)propanoate", id="carboxylate_of_a_chain"),
        pytest.param("[O-]C(=O)c1ccc(C#[N+]C)cc1", "4-[(methylazaniumylidyne)methyl]benzoate", id="nitrilium_beside_a_carboxylate"),
        pytest.param("[O-]C(=O)CCSC(=[OH+])CC", "3-{[1-(oxidaniumylidene)propyl]sulfanyl}propanoate", id="protonated_acyl_beside_a_carboxylate"),
        pytest.param("[O-]c1ccc(cc1)[N+]#N", "4-(diazyn-1-ium-1-yl)phenoxide", id="phenoxide"),
        pytest.param("[O-]S(=O)(=O)c1ccc(cc1)[N+]#N", "4-(diazyn-1-ium-1-yl)benzene-1-sulfonate", id="sulfonate_outranks_the_diazonium_cation"),
        pytest.param("[O-][N+](=O)c1ccc(cc1)[N+]#N", "4-nitrobenzenediazonium", id="nitro_charges_are_not_an_anion"),
    ],
)
def test_diazonium_prefix_on_an_anionic_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-73.1.1.2, P-73.6: the cationic ring nitrogen takes the lowest locant when another ring atom is substituted
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1C[NH2+]CCN1CC", "4-ethylpiperazin-1-ium", id="protonated_ring_with_a_substituent_ahead_of_the_marker"),
        pytest.param("[CH2-][NH+]1CCN(CC)CC1", "(4-ethylpiperazin-1-ium-1-yl)methanide", id="ylide_on_a_ring_nitrogen_of_a_diazine"),
        pytest.param("CC(=O)N1CC[NH2+]CC1", "4-acetylpiperazin-1-ium", id="prefix_cited_before_the_centre_marker_alphabetically"),
        pytest.param("CC(=O)N1CC[NH+](C)CC1", "4-acetyl-1-methylpiperazin-1-ium", id="substituted_cationic_nitrogen_beside_an_acyl_group"),
        pytest.param("CC(=O)N1CCC[NH2+]CC1C", "4-acetyl-3-methyl-1,4-diazepan-1-ium", id="diazepane_with_an_acyl_group_and_a_methyl_group"),
    ],
)
def test_cationic_ring_nitrogen_of_a_substituted_diazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-16.3.6(a), P-16.3.6(e): bis/tris before an anion name that 'di' would merge into another structure
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[Ru+4].[O-]C([O-])=O.[O-]C([O-])=O", "ruthenium(4+) bis(carbonate)", id="mononuclear_oxoanion"),
        pytest.param("[Al+3].[I-].[I-].[I-]", "aluminium tris(iodide)", id="iodide_against_triiodide"),
        pytest.param("[Mn+2].[N-]=C=S.[N-]=C=S", "manganese(2+) bis(sulfanylidenemethaniminide)", id="anion_opening_with_a_prefix"),
        pytest.param("[Cu+2].[O-]C(=O)C(Cl)(Cl)Cl.[O-]C(=O)C(Cl)(Cl)Cl", "copper(2+) bis(trichloroacetate)", id="prefix_with_its_own_multiplier"),
    ],
)
def test_bis_before_an_anion_that_di_would_misname(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[N@+](CC)(CCC)c1ccccc1", "(S)-N-ethyl-N-methyl-N-propylanilinium", id="azanium_beside_a_ring"),
        pytest.param("C(C1=CC=CC=C1)[N@+](C1=CC=CC=C1)(CC=C)C.[Br-]", "(R)-N-benzyl-N-methyl-N-(prop-2-en-1-yl)anilinium bromide", id="descriptor_stays_with_the_cation_of_a_salt"),
        pytest.param("C[N@@+](CC)(CCC)CC(C)C", "(R)-N-ethyl-N,2-dimethyl-N-propylpropan-1-aminium", id="acyclic_aminium"),
        pytest.param("C[P@+](CC)(CCC)c1ccccc1.[Br-]", "(S)-ethyl(methyl)(phenyl)(propyl)phosphanium bromide", id="acyclic_phosphanium"),
        pytest.param("C[P@@+](CC)(CCC)c1ccccc1", "(R)-ethyl(methyl)(phenyl)(propyl)phosphanium", id="phosphanium_without_counter_ion"),
        pytest.param("C[P@H+](CC)CCC", "(R)-ethyl(methyl)(propyl)phosphanium", id="hydrogen_is_the_fourth_ligand"),
        pytest.param("C[As@+](CC)(CCC)c1ccccc1", "(S)-ethyl(methyl)(phenyl)(propyl)arsanium", id="arsanium_like_phosphanium"),
        pytest.param("C1CC2=CC=CC=C2C[P@@+]1(c1ccccc1)c1ccc(O)cc1.[Cl-]", "(2S)-2-(4-hydroxyphenyl)-2-phenyl-1,2,3,4-tetrahydroisophosphinolin-2-ium chloride", id="ring_phosphonium_centre_is_a_hydro_position"),
    ],
)
def test_cip_descriptor_of_a_stereogenic_azanium_or_phosphanium_centre(smiles, expected):
    # P-93.2.2, P-93.5.4
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C@@H]([O-])[C@H](O)[C@H](C)O", "(2R,3S,4S)-3,4-dihydroxypentan-2-olate", id="olate_beside_an_alcohol"),
        pytest.param("[O-][C@@H]1[C@@H](O)[C@@H](O)CO[C@H]1C", "(2S,3R,4R,5S)-4,5-dihydroxy-2-methyloxan-3-olate", id="olate_on_a_ring"),
        pytest.param("C[C@H]([NH3+])[C@H](O)[C@H](C)N", "(2S,3R,4S)-4-amino-3-hydroxypentan-2-aminium", id="aminium_beside_an_amine"),
    ],
)
def test_charged_ligand_ranks_below_its_neutral_twin_in_cip_descriptors(smiles, expected):
    # P-92.1.1, P-92.2.1
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC[S+](C)CCOCCOCC[S+](C)CC", "3,12-dimethyl-6,9-dioxa-3,12-dithioniatetradecane", id="onia_on_a_chain_with_a_multiplied_suffix"),
        pytest.param("C[N+]12CCC(CC1)C2.[Cl-]", "1-methyl-1-azoniabicyclo[2.2.1]heptane chloride", id="onia_on_von_baeyer_with_a_counter_ion"),
        pytest.param("C1CCC[As+]12CCCC2", "5-arsoniaspiro[4.4]nonane", id="ylium_on_a_lambda_parent_becomes_onia"),
        pytest.param("C1CC2CCC1[N+]2", "7-azanyliabicyclo[2.2.1]heptane", id="ylia"),
        pytest.param("CCOP(=O)(C1CCCCC1)OCC[S+](C)CCCCCC", "4-cyclohexyl-8-methyl-4-oxo-3,5-dioxa-8-thionia-4λ5-phosphatetradecane", id="neutral_lambda_heteroatom_beside_onia"),
        pytest.param("C[B-]1(C)CC2(CCCCC2)CC1", "2,2-dimethyl-2-boranuidaspiro[4.5]decane", id="uida"),
        pytest.param("C1CCCC[PH2-]12CCCCC2", "6λ5-phosphanuidaspiro[5.5]undecane", id="uida_keeps_the_lambda_parent"),
        pytest.param("C1CC2CCC1[P-]2", "7-phosphanidabicyclo[2.2.1]heptane", id="ida"),
        pytest.param("C[N+]12CCC(CC1)C2O", "1-methyl-1-azoniabicyclo[2.2.1]heptan-7-ol", id="hydroxy_becomes_the_suffix"),
        pytest.param("C[N+]12CC(=O)C(CC1)C2", "1-methyl-1-azoniabicyclo[2.2.1]heptan-3-one", id="oxo_becomes_the_suffix"),
        pytest.param(
            "C[N+]12CCC(CC1)C2C(=O)O", "1-methyl-1-azoniabicyclo[2.2.1]heptane-7-carboxylic acid", id="carboxy_becomes_the_suffix"
        ),
        pytest.param("O[C@H]1C[N+]2(C)CCC1CC2", "(3R)-1-methyl-1-azoniabicyclo[2.2.2]octan-3-ol", id="descriptor_stays_with_the_ring"),
        pytest.param("C[B-]1(C)CC2(CCC(O)CC2)CC1", "2,2-dimethyl-2-boranuidaspiro[4.5]decan-8-ol", id="anionic_prefix_beside_a_suffix"),
        pytest.param("C[N+]12CCN(CC1)C2", "1-methyl-1-azonia-4-azabicyclo[2.2.1]heptane", id="charged_and_neutral_atoms_of_one_element"),
        pytest.param(
            "C[N+]12CCN(CC1)C(O)C2",
            "4-methyl-4-azonia-1-azabicyclo[2.2.2]octan-2-ol",
            id="numbering_by_the_suffix_instead_of_the_ionic_ending",
        ),
        pytest.param(
            "C[N+]12CCC(CC1)C2.[Br-].[Na+]",
            "1-methyl-1-azoniabicyclo[2.2.1]heptane—bromide—sodium(1+) (1/1/1)",
            id="component_of_a_multi_component_adduct",
        ),
        pytest.param(
            "C1=C2C3=C([P+]4(C2=CC=C1)C1=C(C2=CC=CC=C24)C=CC=C1)C=CC=C3",
            "9-phosphonia-9,9′-spirobi[fluorene]",
            id="spiro_atom_of_two_fused_components",
        ),
        pytest.param(
            "C[P+]1(C)c2ccc(O)cc2-c2ccccc21",
            "9,9-dimethyl-9-phosphoniafluoren-3-ol",
            id="fused_system_with_the_indicated_hydrogen_on_the_replaced_atom",
        ),
        pytest.param(
            "C[As+]1(C)c2ccccc2Oc2ccccc21",
            "9,9-dimethyl-10-oxa-9-arsoniaanthracene",
            id="hydro_prefixes_of_two_replaced_atoms_and_the_element_order",
        ),
        pytest.param(
            "c1ccc2c(c1)C[N+]1(C2)Cc2ccccc2C1",
            "1,1′,3,3′-tetrahydro-2-azonia-2,2′-spirobi[indene]",
            id="hydro_prefixes_stay_before_the_replacement_prefix",
        ),
    ],
)
def test_ionic_skeletal_replacement_prefixes(smiles, expected):
    # P-15.4.1.4, P-72.4 method (2), P-73.4 method (2)
    assert [a.name for a in smiles_to_alternative_names(smiles) if a.kind == "ionic-prefix"] == [expected]


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C[N+]12CCC(CC1)C2CO", id="senior_group_on_a_side_chain"),
        pytest.param("C[N+]12CCC(CC1)C2NC", id="substituted_amine_as_the_senior_group"),
        pytest.param("c1ccc2c(c1)[S+](C)c1ccccc12", id="aromatic_heteroatom_of_a_fusion_name"),
        pytest.param("c1ccc2c(c1)[S+]=CCSC=C2", id="multiply_bonded_heteroatom_of_a_fused_ring"),
        pytest.param("C[N+]1(C)CCCCC1", id="no_replacement_prefix_in_the_name"),
        pytest.param("C1CC2CCC1[Si-]2(C)C", id="bonding_number_without_an_ionic_prefix"),
    ],
)
def test_no_ionic_prefix_name_without_a_skeletal_replacement_name(smiles):
    assert [a for a in smiles_to_alternative_names(smiles) if a.kind == "ionic-prefix"] == []


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[CH2]CC(=O)[O-].[Na+]", "sodium 2-carboxylatoethyl", id="radical_on_a_carboxylate_anion"),
        pytest.param("[CH2]C[N+](C)(C)C.[Cl-]", "2-(trimethylazaniumyl)ethyl chloride", id="radical_on_an_ammonium_cation"),
        pytest.param("[CH2]C[B-](C)(C)C.[Na+]", "sodium 2-(trimethylboranuidyl)ethyl", id="radical_on_a_boranuide_anion"),
        pytest.param("[CH2]CC(=O)[O-].[K+]", "potassium 2-carboxylatoethyl", id="other_counter_cation"),
    ],
)
def test_salt_cites_the_radical_centre_of_an_ionic_fragment(smiles, expected):
    # P-75.4: the radical parent carries the ionic centres as prefixes; the Blue Book names no salt of such an ion
    with pytest.warns(NonPreferredNameWarning, match="P-75.4"):
        assert smiles_to_iupac(smiles) == expected


def test_uide_centre_shared_by_two_ring_components_is_refused_not_overflowed():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cc[n+]2c(c1)CO[B-]21OCCO1")


def test_atom_whose_hydride_is_still_a_radical_is_refused_not_overflowed():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C.[Md]")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1C=[Zr]", "benzylidenezirconium", id="aromatic_ylidene_ligand"),
        pytest.param("Cl[CH]=[Zr]", "chloromethylidenezirconium", id="heteroatom_on_the_ylidene_carbon"),
    ],
)
def test_ylidene_ligands_with_aromatic_or_heteroatom_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
