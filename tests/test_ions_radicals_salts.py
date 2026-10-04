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


def test_alkoxide_with_halogen_substituent():
    assert smiles_to_iupac("FC[O-]") == "fluoromethanolate"


def test_alkoxide_with_unsaturation():
    assert smiles_to_iupac("C=CC[O-]") == "prop-2-en-1-olate"


def test_tert_butoxide():
    assert smiles_to_iupac("CC(C)(C)[O-]") == "tert-butoxide"


def test_branched_alkoxide():
    assert smiles_to_iupac("CC(C)C[O-]") == "2-methylpropan-1-olate"


def test_phenoxide():
    assert smiles_to_iupac("[O-]c1ccccc1") == "phenoxide"
    assert smiles_to_iupac("[O-]c1ccc(C)cc1") == "4-methylphenoxide"
    assert smiles_to_iupac("[O-]c1ccccc1Cl") == "2-chlorophenoxide"


def test_phenoxide_non_alkyl_ring_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(O)c1ccccc1[O-].[Cu+]")


def test_phenoxide_unsaturated_ring_substituent_is_named():
    assert smiles_to_iupac("[O-]c1ccc(CC=C)cc1") == "4-(prop-2-en-1-yl)phenoxide"


def test_two_alkoxide_groups_is_named():
    assert smiles_to_iupac("[O-]CC[O-]") == "ethane-1,2-bis(olate)"


def test_ether_oxygen_alongside_alkoxide_is_named():
    assert smiles_to_iupac("[O-]CCOC") == "2-methoxyethanolate"


def test_alkoxide_stereocenter_with_coexisting_halogen():
    assert smiles_to_iupac("CC[C@@H](Cl)[O-]") == "(1R)-1-chloropropan-1-olate"


def test_phenyl_chain_alkoxide():
    assert smiles_to_iupac("c1ccccc1CC[O-]") == "2-phenylethanolate"
    assert smiles_to_iupac("c1ccccc1CCC[O-]") == "3-phenylpropan-1-olate"


def test_phenyl_substituted_benzene_ring_alkoxide_is_named():
    assert smiles_to_iupac("Cc1ccccc1CC[O-]") == "2-(2-methylphenyl)ethanolate"


def test_phenyl_chain_alkoxide_unsaturation_is_named():
    assert smiles_to_iupac("C=Cc1ccccc1CC[O-]") == "2-(2-ethenylphenyl)ethanolate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[NH-]", "propan-2-aminide"),
    ],
)
def test_aminide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aromatic_aminide_is_named():
    assert smiles_to_iupac("[NH-]c1ccccc1") == "benzenaminide"


def test_two_aminide_groups_is_named():
    assert smiles_to_iupac("[NH-]CC[NH-]") == "ethane-1,2-bis(aminide)"


def test_enamine_aminide_is_named():
    assert smiles_to_iupac("C=C[NH-]") == "eth-1-en-1-aminide"


def test_second_nitrogen_is_named():
    assert smiles_to_iupac("NCC[NH-]") == "2-aminoethanaminide"


def test_phenyl_chain_aminide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC[NH-])cc1") == "3-(4-chlorophenyl)propan-1-aminide"


def test_phenyl_chain_aminide_unsaturation_is_named():
    assert smiles_to_iupac("C=Cc1ccccc1CC[NH-]") == "2-(2-ethenylphenyl)ethanaminide"


def test_ring_ammonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[NH2+]C1CCCCC1")


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


def test_benzenide_name():
    assert smiles_to_iupac("[c-]1ccccc1") == "benzenide"


def test_carbanide_with_ring_is_named():
    assert smiles_to_iupac("[CH-]1CCCCC1") == "cyclohexanide"


def test_carbanide_with_halogen_is_named():
    assert smiles_to_iupac("ClC[CH2-]") == "2-chloroethan-1-ide"


def test_carbanide_with_unsaturation_is_named():
    assert smiles_to_iupac("C=CC[CH2-]") == "but-3-en-1-ide"


def test_acetyl_anion_with_oxo_on_anion_carbon():
    assert smiles_to_iupac("C[C-]=O") == "1-oxoethan-1-ide"


def test_oxo_substituent_elsewhere_on_chain():
    assert smiles_to_iupac("CCC(=O)[CH-]C") == "3-oxopentan-2-ide"


def test_carbanide_with_second_ketone_is_named():
    assert smiles_to_iupac("O=CC(=O)[CH-]C") == "3,4-dioxobutan-2-ide"


def test_carbanide_with_aldehyde_shaped_carbonyl_is_named():
    assert smiles_to_iupac("O=C[CH-]C") == "1-oxopropan-2-ide"


def test_methylium_name():
    # Blue Book P-73.2.2.1.1 worked example: "[CH3]+ -> methylium (PIN)".
    assert smiles_to_iupac("[CH3+]") == "methylium"


def test_ethylium_name():
    assert smiles_to_iupac("[CH2+]C") == "ethylium"


def test_cyclobutylium_name():
    assert smiles_to_iupac("C1C[CH+]C1") == "cyclobutylium"


def test_branch_point_carbenium_butan_2_ylium_name():
    assert smiles_to_iupac("C[CH+]CC") == "butan-2-ylium"


def test_branch_point_carbenium_three_branches_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C+](C)C")


def test_branch_point_carbenium_with_further_branching_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[CH+]C(C)C")


def test_branched_chain_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C(C)C")


def test_substituted_ring_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC[CH+]C1")


def test_halogen_substituted_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C(Cl)")


def test_unsaturated_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C=C")


def test_acetylium_acylium_cation_name():
    assert smiles_to_iupac("C[C+]=O") == "acetylium"


def test_branched_acylium_name():
    assert smiles_to_iupac("CC(C)C[C+]=O") == "3-methylbutanoylium"


def test_cyclohexanecarbonylium_acylium_name():
    assert smiles_to_iupac("O=[C+]C1CCCCC1") == "cyclohexanecarbonylium"


def test_benzoylium_acylium_name():
    assert smiles_to_iupac("O=[C+]c1ccccc1") == "benzoylium"


def test_cyclopentadienide_name():
    assert smiles_to_iupac("[CH-]1C=CC=C1") == "cyclopenta-2,4-dien-1-ide"


def test_cyclopentadiene_neutral_parent_unaffected():
    assert smiles_to_iupac("C1=CC=CC1") == "cyclopenta-1,3-diene"


def test_diphosphoric_acid():
    assert smiles_to_iupac("OP(=O)(O)OP(=O)(O)O") == "diphosphoric acid"


def test_nitrone_unsubstituted_nitrogen():
    assert smiles_to_iupac("C=[N+]([H])[O-]") == "methanimine N-oxide"


def test_nitrile_oxide_propane():
    assert smiles_to_iupac("CCC#[N+][O-]") == "propanenitrile oxide"


def test_nitrile_oxide_fulminic_acid():
    assert smiles_to_iupac("C#[N+][O-]") == "formonitrile oxide"


def test_nitrile_imide_still_out_of_scope():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC#[N+][N-]C")


def test_two_halide_fragments_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCCN.Cl.Cl")


def test_plain_mixture_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCO.CCO")


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


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[OH2+]")


def test_ring_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[OH2+]C1CCCCC1")


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


def test_branched_chain_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C(C)C")


def test_branch_point_radical_butan_2_yl_name():
    assert smiles_to_iupac("C[CH]CC") == "butan-2-yl"


def test_branch_point_radical_tert_butyl_name():
    assert smiles_to_iupac("[C](C)(C)C") == "tert-butyl"


def test_branch_point_radical_2_methylbutan_2_yl_name():
    assert smiles_to_iupac("CC[C](C)C") == "2-methylbutan-2-yl"


def test_branch_point_radical_with_further_branching_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH](C(C)C)C")


def test_substituted_ring_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC[CH]C1")


def test_two_radical_centers_diyl_name():
    assert smiles_to_iupac("[CH2][CH2]") == "ethane-1,2-diyl"


def test_three_radical_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2][CH][CH2]")


def test_halogen_substituted_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C(Cl)")


def test_ethylidene_radical_name():
    assert smiles_to_iupac("[CH]C") == "ethylidene"


def test_cyclobutylidene_radical_name():
    assert smiles_to_iupac("[C]1CCC1") == "cyclobutylidene"


def test_branch_point_divalent_radical_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[C](C)C")


def test_branched_acyl_radical_name():
    assert smiles_to_iupac("[C](=O)C(C)CCC") == "2-methylpentanoyl"


def test_cyclohexanecarbonyl_acyl_radical_name():
    assert smiles_to_iupac("[C](=O)C1CCCCC1") == "cyclohexanecarbonyl"


def test_benzoyl_acyl_radical_name():
    # Blue Book P-71.3.1 worked example: the retained 'benzoyl' name.
    assert smiles_to_iupac("[C](=O)c1ccccc1") == "benzoyl"


def test_aminyl_radical_name():
    # Blue Book P-71.3.2 worked example.
    assert smiles_to_iupac("C[NH]") == "methanaminyl"


def test_iminyl_radical_name():
    assert smiles_to_iupac("CC=[N]") == "ethaniminyl"


def test_amidyl_radical_name():
    assert smiles_to_iupac("CC(=O)[NH]") == "ethanamidyl"


def test_vinyl_carbyne_name():
    assert smiles_to_iupac("C#C[CH]") == "prop-2-yn-1-ylidene"


def test_butoxyl_radical_name():
    assert smiles_to_iupac("CCCC[O]") == "butoxyl"


def test_methylselanyl_radical_name():
    assert smiles_to_iupac("C[Se]") == "methylselanyl"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[CH]1CC[CH]CC1", "cyclohexane-1,4-diyl"),
    ],
)
def test_ring_diradical_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diradical_branch_point_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2]C([CH2])C")


def test_ethyloxidaniumyl():
    assert smiles_to_iupac("CC[OH+]") == "ethyloxidaniumyl"


def test_methanaminyliumyl():
    assert smiles_to_iupac("C[N+]") == "methanaminyliumyl"


def test_ethaniminyliumyl():
    assert smiles_to_iupac("CC=[N+]") == "ethaniminyliumyl"


def test_ethanamidyliumyl():
    assert smiles_to_iupac("CC(=O)[N+]") == "ethanamidyliumyl"


def test_calcium_bis_compound_carboxylate():
    assert (
        smiles_to_iupac("[Ca+2].CCCCC(C)(C(=O)[O-])c1ccccc1.CCCCC(C)(C(=O)[O-])c1ccccc1")
        == "calcium bis(2-methyl-2-phenylhexanoate)"
    )


def test_sodium_methoxide():
    assert smiles_to_iupac("[Na+].C[O-]") == "sodium methoxide"


def test_ammonium_ethanethioate():
    assert smiles_to_iupac("[NH4+].CC(=O)[S-]") == "azanium ethanethioate"


def test_calcium_dichloride():
    assert smiles_to_iupac("[Ca+2].[Cl-].[Cl-]") == "calcium dichloride"


def test_mixed_halide_anions_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ca+2].[Cl-].[Br-]")


def test_carbanide_salt_names():
    assert smiles_to_iupac("[CH3-].[Li+]") == "lithium methanide"
    assert smiles_to_iupac("C[CH2-].[Na+]") == "sodium ethanide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Valine zwitterion (branched chain).
        ("CC(C)C(C(=O)[O-])[NH3+]", "2-azaniumyl-3-methylbutanoate"),
        ("C[N+](C)(C)CCS(=O)(=O)[O-]", "2-(N,N-dimethylmethanaminiumyl)ethanesulfonate"),
    ],
)
def test_zwitterion_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_zwitterion_ionic_center_in_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC([NH3+])C1C(=O)[O-]")


def test_zwitterion_ammonium_bonded_directly_to_sulfonate_carbon_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH3+]C(S(=O)(=O)[O-])")


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
