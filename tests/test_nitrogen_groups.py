import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_n_n_disubstituted_amide_with_locant_leading_parent_name():
    assert (
        smiles_to_iupac("O=C(C(C(F)(F)F)C(F)(F)F)N(CC)CC")
        == "N,N-diethyl-3,3,3-trifluoro-2-(trifluoromethyl)propanamide"
    )


def test_halogenated_n_substituent():
    assert smiles_to_iupac("CC(=O)NCCCl") == "N-(2-chloroethyl)ethanamide"


def test_unsaturated_n_substituent():
    assert smiles_to_iupac("CC(=O)NC=C") == "N-ethenylethanamide"


def test_lactam_is_named_via_ketone_module():
    assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"


def test_diamide():
    assert smiles_to_iupac("NC(=O)CC(N)=O") == "propanediamide"


def test_phenyl_chain_amide_with_hydroxyl():
    assert smiles_to_iupac("OCc1ccccc1CC(N)=O") == "2-[2-(hydroxymethyl)phenyl]ethanamide"


def test_phenyl_chain_amide_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(N)=O") == "2-(2-ethenylphenyl)ethanamide"


def test_n_phenyl_formamide():
    assert smiles_to_iupac("O=CNc1ccccc1") == "N-phenylmethanamide"


def test_n_phenyl_amide_with_second_n_substituent():
    assert smiles_to_iupac("CC(=O)N(C)c1ccccc1") == "N-methyl-N-phenylethanamide"


def test_n_phenyl_amide_substituted_ring():
    assert smiles_to_iupac("CC(=O)Nc1ccc(O)cc1") == "N-(4-hydroxyphenyl)ethanamide"


def test_benzamide():
    assert smiles_to_iupac("NC(=O)c1ccccc1") == "benzamide"
    assert smiles_to_iupac("NC(=O)c1ccc(C)cc1") == "4-methylbenzamide"
    assert smiles_to_iupac("NC(=O)c1ccccc1Cl") == "2-chlorobenzamide"


def test_amide_enol_mix():
    assert smiles_to_iupac("OC=CC(N)=O") == "3-hydroxyprop-2-enamide"


def test_n_substituted_amide_stereocenter():
    assert smiles_to_iupac("CC[C@@H](C)C(=O)NCC") == "(2R)-N-ethyl-2-methylbutanamide"


def test_cyclic_amide_ring_stereocenter():
    assert (
        smiles_to_iupac("NC(=O)[C@H]1CCCC[C@@H]1Cl") == "(1R,2S)-2-chlorocyclohexane-1-carboxamide"
    )


def test_amide_ring_branch_stereocenter():
    assert (
        smiles_to_iupac("NC(=O)[C@H]1CC[C@H](C[C@@H](C)CC)C1")
        == "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carboxamide"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)NO", "N-hydroxyethanamide"),
    ],
)
def test_hydroxamic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hydroxamic_acid_on_ring():
    assert smiles_to_iupac("C1CCC(CC1)C(=O)NO") == "N-hydroxycyclohexanecarboxamide"


def test_hydroxamic_acid_with_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(C)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C(=N)N", "methanimidamide"),
        # -imidamide combined with existing unsaturation support.
        ("C=CCC(=N)N", "but-3-enimidamide"),
    ],
)
def test_amidine_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diamidine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=N)CC(=N)N")


def test_amidine_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(=N)N")


def test_n_prime_phenylamidine():
    assert smiles_to_iupac("CC(=Nc1ccccc1)N") == "N'-phenylethanimidamide"


def test_n_prime_methyl_n_phenylamidine():
    assert smiles_to_iupac("CC(=NC)Nc1ccccc1") == "N'-methyl-N-phenylethanimidamide"


def test_n_phenylamidine_with_second_imino_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=Nc1ccccc1)Nc1ccccc1")


def test_phenyl_directly_attached_amidine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=N)N")


def test_phenyl_chain_amidine_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC(=N)N)cc1") == "3-(4-chlorophenyl)propanimidamide"


def test_phenyl_chain_amidine_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=N)N")


def test_halogen_substituted_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC(=N)N(CC)CC(F)(F)F")


def test_geminal_diamine_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(N)N(C)C")


def test_two_amines_with_unsaturated_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCCN(C)CC=C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCNc1ccc(Cl)cc1", "4-chloro-N-ethylaniline"),
        ("c1ccccc1CNc1ccccc1", "N-benzylaniline"),
    ],
)
def test_n_substituted_aniline_and_multi_ring_chain_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_enamine():
    assert smiles_to_iupac("NC=CC") == "prop-1-en-1-amine"


def test_cyclic_amine_branch_stereocenter():
    assert smiles_to_iupac("NC1(CCCCC1)[C@@H](C)CC") == "1-[(2S)-butan-2-yl]cyclohexan-1-amine"


def test_acyclic_primary_amine_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("C/C=C/[C@@H](C)N") == "(2R,3E)-pent-3-en-2-amine"
    assert smiles_to_iupac("C/C=C/[C@H](C)N") == "(2S,3E)-pent-3-en-2-amine"


def test_ammonium_stereocenter_with_ez_double_bond_coexistence():
    assert (
        smiles_to_iupac("C[N+](C)(C)[C@H](C)/C=C/C") == "(2R,3E)-N,N,N-trimethylpent-3-en-2-aminium"
    )


def test_phenyl_chain_amine_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CCN") == "2-(2-ethenylphenyl)ethanamine"


def test_phenyl_chain_diamine():
    assert smiles_to_iupac("c1ccccc1C(N)CCN") == "1-phenylpropane-1,3-diamine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CNCc1ccccc1", "N-benzyl-1-phenylmethanamine"),
        ("CC(C)NCc1ccccc1", "N-benzylpropan-2-amine"),
    ],
)
def test_phenyl_chain_secondary_tertiary_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_amine_with_substituent():
    assert smiles_to_iupac("NC1CCCC=C1C") == "2-methylcyclohex-2-en-1-amine"


def test_unsaturated_ring_amine_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC#C1")


def test_ring_substituent_chain_amine():
    assert smiles_to_iupac("NCC1CCCCC1") == "cyclohexylmethanamine"


def test_ring_with_amine_chain_amine_tie():
    assert smiles_to_iupac("NC1CCCCC1CN") == "2-(aminomethyl)cyclohexan-1-amine"


def test_ring_with_amine_chain_amine_chain_wins():
    assert smiles_to_iupac("NC1CCCCC1C(N)CN") == "1-(2-aminocyclohexyl)ethane-1,2-diamine"


def test_primary_amine_oxide_not_matched():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[NH2+][O-]")


def test_amine_oxide_ring_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[N+]1([O-])CCCCC1")


def test_amine_oxide_stereocenter():
    assert smiles_to_iupac("CC[C@@H](C)[N+](C)(C)[O-]") == "(2R)-N,N-dimethylbutan-2-amine N-oxide"
    assert smiles_to_iupac("CC[C@H](C)[N+](C)(C)[O-]") == "(2S)-N,N-dimethylbutan-2-amine N-oxide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCN=[N+]=[N-]", "azidoethane"),
    ],
)
def test_azide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CCCCCCCCCCN=[N+]=[N-]", "(10-azidodecyl)benzene"),
    ],
)
def test_azide_benzene_ring_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_benzene_ring_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CN=[N+]=[N-]")


def test_fused_aromatic_ring_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2ccccc2c1CN=[N+]=[N-]")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(N=[N+]=[N-])CC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCN=[N+]=[N-]")


def test_acetone_azine_name():
    assert smiles_to_iupac("CC(C)=NN=C(C)C") == "N-(propan-2-ylideneamino)propan-2-imine"


def test_asymmetric_azine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NN=C(C)C")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1=NN=C1CCCCC1")


def test_aromatic_carbon_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C=NN=Cc1ccccc1")


def test_azine_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(C)=NN=C(C)[C@H](C)CC")


def test_methyl_n_n_dimethylcarbamate():
    assert smiles_to_iupac("COC(=O)N(C)C") == "methyl N,N-dimethylcarbamate"


def test_plain_cyclic_r_group():
    assert smiles_to_iupac("O=C(N)OC1CCCCC1") == "cyclohexyl carbamate"
    assert smiles_to_iupac("NC(=O)OC1CCCC1") == "cyclopentyl carbamate"


def test_n_methyl_cyclic_r_group():
    assert smiles_to_iupac("CNC(=O)OC1CCCCC1") == "cyclohexyl N-methylcarbamate"


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCOC(N)=O")


def test_phenyl_r_group_still_not_supported():
    # The phenyl exception only applies to the amide nitrogen, not R.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)Oc1ccccc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)N(C)c1ccccc1")


def test_unsaturated_r_not_supported__cyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=COC#N")


def test_cyclic_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1OC#N")


def test_two_cyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#COCCOC#N")


def test_phenyl_cyanate():
    assert smiles_to_iupac("c1ccccc1OC#N") == "phenyl cyanate"  # CID 70740


def test_phenyl_cyanate_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1COC#N")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N=N", "diazene"),
        ("CN=NCC", "ethyl(methyl)diazene"),
    ],
)
def test_diazene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1N=Nc1ccccc1", "diphenyldiazene"),
    ],
)
def test_ring_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_phenyl():
    assert smiles_to_iupac("Cc1ccccc1N=N") == "(2-methylphenyl)diazene"


def test_bis_chloroethyldiazene_name():
    assert smiles_to_iupac("ClCCN=NCCCl") == "bis(2-chloroethyl)diazene"


def test_halogen_on_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClN=N")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=[N+]=[N-]", "diazomethane"),
    ],
)
def test_diazo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported__diazo():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1(=[N+]=[N-])CCCCC1")


def test_unsaturated_chain_not_supported__diazo():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC=[N+]=[N-]")


def test_two_diazo_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[N-]=[N+]=CC=[N+]=[N-]")


def test_methanediazonium():
    assert smiles_to_iupac("C[N+]#N") == "methanediazonium"


def test_chloroethane_diazonium():
    assert smiles_to_iupac("ClCC[N+]#N") == "2-chloroethanediazonium"


def test_pent_4_ene_1_diazonium():
    assert smiles_to_iupac("C=CCCC[N+]#N") == "pent-4-ene-1-diazonium"


def test_two_diazonium_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#[N+]CC[N+]#N")


def test_cyclohexanediazonium():
    assert smiles_to_iupac("C1CCCCC1[N+]#N") == "cyclohexanediazonium"


def test_substituted_ring_diazonium_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCCCC1[N+]#N")


def test_benzenediazonium():
    assert smiles_to_iupac("c1ccccc1[N+]#N") == "benzenediazonium"  # CID 9718


def test_diazonium_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[N+]#N")


def test_phenyl_chain_diazonium():
    assert smiles_to_iupac("c1ccccc1CC[N+]#N") == "2-phenylethanediazonium"
    assert smiles_to_iupac("c1ccccc1CCC[N+]#N") == "3-phenylpropane-1-diazonium"


def test_dimethyl_benzenediazonium():
    assert smiles_to_iupac("Cc1ccc(C)cc1[N+]#N") == "2,5-dimethylbenzenediazonium"


def test_phenyl_substituted_benzene_ring_diazonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[N+]#N")


def test_phenyl_chain_diazonium_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[N+]#N")


def test_guanidine():
    assert smiles_to_iupac("NC(=N)N") == "guanidine"


def test_different_substituent_counts_on_amino_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(C)C(=N)NC")


def test_imino_plus_amino_different_names_alphabetized():
    assert smiles_to_iupac("CNC(=NCC)N") == "N''-ethyl-N-methylguanidine"


def test_unsaturated_imino_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=NC=C)N")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=N)N")


def test_n_n_prime_diphenylguanidine():
    assert smiles_to_iupac("c1ccccc1NC(=N)Nc1ccccc1") == "N,N'-diphenylguanidine"


def test_n_double_prime_phenylguanidine():
    assert smiles_to_iupac("NC(=Nc1ccccc1)N") == "N''-phenylguanidine"


def test_substituted_phenyl_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc(NC(=N)N)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported__guanidine():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=N)N")


def test_ring_fused_guanidine_still_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N1CCCCC1=N")


def test_pyrrolidine_amine_2():
    assert smiles_to_iupac("NC1CCCN1") == "pyrrolidin-2-amine"


def test_secondary_exocyclic_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC1CCNCC1")


def test_two_exocyclic_amines():
    assert smiles_to_iupac("NC1CCNC(N)C1") == "piperidine-2,4-diamine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # -hydrazide combined with existing unsaturation support.
        ("CC=CC(=O)NN", "but-2-enehydrazide"),
        ("OCCC(=O)NN", "3-hydroxypropanehydrazide"),
    ],
)
def test_hydrazide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formohydrazide():
    assert smiles_to_iupac("C(=O)NN") == "formohydrazide"


def test_acetohydrazide():
    assert smiles_to_iupac("CC(=O)NN") == "acetohydrazide"


def test_n_substituted_stereocenter_hydrazide():
    assert smiles_to_iupac("CC[C@@H](C)C(=O)N(C)N") == "(2R)-N-methyl-2-methylbutanehydrazide"


def test_n_substituted_formohydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CN(C)N")


def test_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(N)C(C)C")


def test_diacylhydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCC(=O)NNC(=O)C")


def test_phenyl_chain_hydrazide_retained_name():
    assert smiles_to_iupac("c1ccccc1CC(=O)NN") == "2-phenylacetohydrazide"


def test_phenyl_directly_attached_hydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)NN")


def test_n_prime_phenylhydrazide():
    assert smiles_to_iupac("CC(=O)NNc1ccccc1") == "N'-phenylacetohydrazide"
    assert smiles_to_iupac("CCC(=O)NNc1ccccc1") == "N'-phenylpropanehydrazide"


def test_n_phenylhydrazide():
    assert smiles_to_iupac("CC(=O)N(c1ccccc1)N") == "N-phenylacetohydrazide"


def test_n_phenylhydrazide_with_second_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(c1ccccc1)Nc1ccccc1")


def test_phenyl_chain_hydrazide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC(=O)NN)cc1") == "3-(4-chlorophenyl)propanehydrazide"


def test_phenyl_chain_hydrazide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)NN")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)NNC(=O)C", "N'-acetylacetohydrazide"),
        ("CCC(=O)NNC(=O)CC", "N'-propanoylpropanehydrazide"),
    ],
)
def test_symmetric_diacyl_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diacyl_hydrazide_formyl():
    # Mononuclear retained name on both sides.
    assert smiles_to_iupac("O=CNNC=O") == "N'-formylformohydrazide"


def test_diacyl_hydrazide_aromatic_acyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NNC(=O)c1ccccc1")


def test_diacyl_hydrazide_branched_acyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)NNC(=O)C(C)C")


def test_halogen_substituted_n_alkyl_raises__hydrazide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC(=O)NNCC(F)(F)F")


def test_plain_hydrazine_name():
    assert smiles_to_iupac("NN") == "hydrazine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1NNc1ccccc1", "1,2-diphenylhydrazine"),
    ],
)
def test_ring_substituent__hydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_phenyl__hydrazine():
    assert smiles_to_iupac("Cc1ccccc1NN") == "(2-methylphenyl)hydrazine"


def test_halogen_on_nitrogen_raises__hydrazine():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClNN")


def test_n_substituted_hydrazone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NNC")


def test_ring_raises__hydrazone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1=NN")


def test_aromatic_carbon_raises__hydrazone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C=NN")


def test_hydrazone_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)C(C)=NN") == "3-methylpentan-2-ylidenehydrazine"


def test_hydrazone_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(C)=NN")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C/C=N/N", "(E)-ethylidenehydrazine"),
    ],
)
def test_specified_double_bond_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NO", "hydroxylamine"),
    ],
)
def test_hydroxylamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNO")


def test_n_o_disubstituted_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNOC")


def test_ring_not_supported__hydroxylamine():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ONC1CCCCC1")


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CON")


def test_phenyl_o_substituent():
    assert smiles_to_iupac("NOc1ccccc1") == "O-phenylhydroxylamine"


def test_benzyl_o_substituent_still_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NOCc1ccccc1")


def test_substituted_phenyl_o_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NOc1ccccc1C")


def test_branched_symmetric_imide():
    assert smiles_to_iupac("CC(C)C(=O)NC(=O)C(C)C") == "N-(2-methylpropanoyl)-2-methylpropanamide"


def test_halogen_substituent():
    assert smiles_to_iupac("ClCC(=O)NC(=O)CCl") == "N-(2-chloroethanoyl)-2-chloroethanamide"


def test_unsymmetric_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC(=O)CC")


def test_unsaturated_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)NC(=O)C=C")


def test_n_substituted_cyclic_imide_named_via_ketone_module():
    assert smiles_to_iupac("CN1C(=O)CCC1=O") == "1-methylpyrrolidine-2,5-dione"


def test_imide_phenyl_chain_unsymmetric_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)NC(=O)CCc1ccccc1")


def test_imide_phenyl_ring_directly_on_acyl_carbon_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)NC(=O)c1ccccc1")


def test_imide_one_sided_phenyl_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)NC(=O)CC")


def test_imide_substituted_benzene_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)NC(=O)Cc1ccccc1C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC=NOCC", "N-ethoxypropan-1-imine"),
    ],
)
def test_smiles_to_iupac_simple_imine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_raises__imine():
    # A cyclic/aromatic imine (e.g. 'thiolan-2-imine') is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N=C1CCCC1")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC=N")


def test_imine_with_a_senior_alcohol_is_an_imino_prefix():
    assert smiles_to_iupac("OCC=N") == "2-iminoethanol"


def test_branched_oxime_o_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NOC(C)C")


def test_two_oxygens_on_oxime_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NOO")


def test_imine_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(C)=N")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C/C=N/O", "(1E)-N-hydroxyethanimine"),
    ],
)
def test_specified_double_bond_stereo_oxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC/C=N/C", "(1E)-N-methyl-3-phenylpropan-1-imine"),
    ],
)
def test_specified_double_bond_stereo_n_alkyl_imine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_imine():
    assert smiles_to_iupac("c1ccccc1CC=N") == "2-phenylethanimine"
    assert smiles_to_iupac("c1ccccc1CCC=N") == "3-phenylpropan-1-imine"
    assert smiles_to_iupac("c1ccccc1CC(C)=N") == "1-phenylpropan-2-imine"


def test_phenyl_chain_oxime_ether_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCC=NOCC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)N=C=O", "2-isocyanatopropane"),
    ],
)
def test_isocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported__isocyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(N=C=O)CC1")


def test_unsaturated_chain_not_supported__isocyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCN=C=O")


def test_two_isocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C=NCN=C=O")


def test_benzene_ring_direct_bond():
    assert smiles_to_iupac("c1ccccc1N=C=O") == "isocyanatobenzene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[N+]#[C-]", "isocyanoethane"),
    ],
)
def test_isocyanide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported__isocyanide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N+]#[C-])CC1")


def test_phenyl_isocyanide_chain():
    assert smiles_to_iupac("c1ccccc1C[N+]#[C-]") == "(isocyanomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CC[N+]#[C-]") == "(2-isocyanoethyl)benzene"


def test_phenyl_isocyanide_unsaturation_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1C[N+]#[C-]")


def test_unsaturated_chain_not_supported__isocyanide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N+]#[C-]")


def test_two_isocyanide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[C-]#[N+]C[N+]#[C-]")


def test_multiple_fragments_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CO[N+](=O)[O-].C")


def test_two_nitrate_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-][N+](=O)OCCCO[N+](=O)[O-]")


def test_nitro_unaffected():
    assert smiles_to_iupac("C[N+](=O)[O-]") == "nitromethane"


def test_formonitrile_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C#N")


def test_trinitrile():
    assert smiles_to_iupac("N#CC(CC#N)CC#N") == "propane-1,2,3-tricarbonitrile"


def test_dinitrile_alongside_ring():
    assert smiles_to_iupac("N#CC1CCC(C#N)CC1") == "cyclohexane-1,4-dicarbonitrile"


def test_benzonitrile():
    assert smiles_to_iupac("N#Cc1ccccc1") == "benzonitrile"
    assert smiles_to_iupac("N#Cc1ccc(C)cc1") == "4-methylbenzonitrile"
    assert smiles_to_iupac("N#Cc1ccccc1Cl") == "2-chlorobenzonitrile"


def test_phenyl_chain_nitrile_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC#N") == "2-(2-ethenylphenyl)ethanenitrile"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N#CC1CCCCC1c1cccs1", "2-(thiophen-2-yl)cyclohexane-1-carbonitrile"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_nitrile(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_aromatic_substituent_nitrile_substituted_ring():
    assert (
        smiles_to_iupac("N#CC1CCCC(c2ccc(C)cc2)C1")
        == "3-(4-methylphenyl)cyclohexane-1-carbonitrile"
    )


def test_two_ring_aromatic_substituent_nitrile_chain_nitrile():
    assert smiles_to_iupac("N#CCC1CCCCC1c1ccccc1") == "2-(2-phenylcyclohexyl)ethanenitrile"


def test_heteroaromatic_ring_directly_attached_nitrile():
    assert smiles_to_iupac("N#Cc1cccnc1") == "pyridine-3-carbonitrile"


def test_acyclic_nitrile_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("N#C[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enenitrile"
    assert smiles_to_iupac("N#C[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enenitrile"


def test_cyclic_nitrile_ring_stereocenter():
    assert smiles_to_iupac("N#C[C@H]1CCCC[C@@H]1Cl") == "(1R,2S)-2-chlorocyclohexane-1-carbonitrile"


def test_nitrile_ring_branch_stereocenter():
    assert (
        smiles_to_iupac("N#C[C@H]1CC[C@H](C[C@@H](C)CC)C1")
        == "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carbonitrile"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCON=O", "ethyl nitrite"),  # CID 8026
    ],
)
def test_nitrite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_fragments_raises__nitrite_ester():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CON=O.C")


def test_two_nitrite_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=NOCCCON=O")


def test_nitrate_ester_unaffected():
    assert smiles_to_iupac("CO[N+](=O)[O-]") == "methyl nitrate"


def test_ring():
    assert smiles_to_iupac("C1CCC([N+](=O)[O-])CC1") == "nitrocyclohexane"


def test_unsaturated_chain():
    assert smiles_to_iupac("C=CC[N+](=O)[O-]") == "3-nitroprop-1-ene"


def test_phenyl_nitro_chain():
    assert smiles_to_iupac("c1ccccc1CC[N+](=O)[O-]") == "(2-nitroethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CCCCCCC[N+](=O)[O-]") == "(7-nitroheptyl)benzene"


def test_phenyl_nitro_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1C[N+](=O)[O-]") == "1-ethenyl-2-(nitromethyl)benzene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCN=O", "nitrosoethane"),
    ],
)
def test_nitroso(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring__nitroso():
    assert smiles_to_iupac("C1CCCCC1N=O") == "nitrosocyclohexane"


def test_unsaturated_chain__nitroso():
    assert smiles_to_iupac("C=CN=O") == "nitrosoethene"


def test_phenyl_nitroso_chain():
    assert smiles_to_iupac("c1ccccc1CN=O") == "(nitrosomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CCN=O") == "(2-nitrosoethyl)benzene"


def test_phenyl_nitroso_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CN=O") == "1-ethenyl-2-(nitrosomethyl)benzene"


def test_methylazepane():
    # 7-membered ring analogue of the piperidine/pyrrolidine cases above.
    assert smiles_to_iupac("CN1CCCCCC1") == "1-methylazepane"


def test_piperazine():
    assert smiles_to_iupac("CN1CCNCC1") == "1-methylpiperazine"


def test_unsaturated_n_substituent__ring_amine():
    assert smiles_to_iupac("C=CN1CCCCC1") == "1-ethenylpiperidine"


def test_ring_amine_two_ring_carbon_substituents():
    assert smiles_to_iupac("CN1CC(C)C(C)C1") == "1,3,4-trimethylpyrrolidine"


def test_ring_amine_sulfonyl_with_ring_carbon_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)(=O)N1CCC(C)CC1")


def test_methylsulfonylmorpholine():
    assert smiles_to_iupac("CS(=O)(=O)N1CCOCC1") == "4-methylsulfonylmorpholine"


def test_urea():
    assert smiles_to_iupac("NC(=O)N") == "urea"


def test_n_methylurea():
    assert smiles_to_iupac("CNC(=O)N") == "N-methylurea"


def test_n_n_dimethylurea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=O)N") == "N,N-dimethylurea"


def test_n_ethyl_n_methylurea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=O)N") == "N-ethyl-N-methylurea"


def test_different_substituent_counts_on_different_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCN(C)C(=O)NC")


def test_unsaturated_n_substituent_not_supported__urea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=O)N")


def test_n_methyl_n_prime_phenylurea_different_nitrogens():
    assert smiles_to_iupac("CNC(=O)Nc1ccccc1") == "N-methyl-N'-phenylurea"


def test_n_n_prime_diphenylurea():
    assert smiles_to_iupac("c1ccc(NC(=O)Nc2ccccc2)cc1") == "N,N'-diphenylurea"


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported__urea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=O)N")


def test_semicarbazide():
    assert smiles_to_iupac("NC(=O)NN") == "aminourea"


def test_semicarbazide_with_n_alkyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=O)NN")


def test_double_amino_substituted_urea_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NNC(=O)NN")
