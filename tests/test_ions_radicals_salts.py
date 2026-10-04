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


def test_phenoxide_unsaturated_ring_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]c1ccc(CC=C)cc1")


def test_two_alkoxide_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]CC[O-]")


def test_ether_oxygen_alongside_alkoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]CCOC")


def test_alkoxide_stereocenter_with_coexisting_halogen():
    assert smiles_to_iupac("CC[C@@H](Cl)[O-]") == "(1R)-1-chloropropan-1-olate"


def test_phenyl_chain_alkoxide():
    assert smiles_to_iupac("c1ccccc1CC[O-]") == "2-phenylethanolate"
    assert smiles_to_iupac("c1ccccc1CCC[O-]") == "3-phenylpropan-1-olate"


def test_phenyl_substituted_benzene_ring_alkoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[O-]")


def test_phenyl_chain_alkoxide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[O-]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[NH-]", "propan-2-aminide"),
    ],
)
def test_aminide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aromatic_aminide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH-]c1ccccc1")


def test_two_aminide_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH-]CC[NH-]")


def test_enamine_aminide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[NH-]")


def test_second_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC[NH-]")


def test_phenyl_chain_aminide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC[NH-])cc1") == "3-(4-chlorophenyl)propan-1-aminide"


def test_phenyl_chain_aminide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[NH-]")


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


def test_carbanide_with_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH-]1CCCCC1")


def test_carbanide_with_halogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC[CH2-]")


def test_carbanide_with_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[CH2-]")


def test_acetyl_anion_with_oxo_on_anion_carbon():
    assert smiles_to_iupac("C[C-]=O") == "1-oxoethan-1-ide"


def test_oxo_substituent_elsewhere_on_chain():
    assert smiles_to_iupac("CCC(=O)[CH-]C") == "3-oxopentan-2-ide"


def test_carbanide_with_second_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC(=O)[CH-]C")


def test_carbanide_with_aldehyde_shaped_carbonyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C[CH-]C")


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
        == "butyl 2-[(2-phenylethanoyl)amino]benzoate;hydrochloride"
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
