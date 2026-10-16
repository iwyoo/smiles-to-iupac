import pytest
from smiles_to_iupac import smiles_to_alternative_names, smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C(C(C(F)(F)F)C(F)(F)F)N(CC)CC", "N,N-diethyl-3,3,3-trifluoro-2-(trifluoromethyl)propanamide", id="n_n_disubstituted_amide_with_locant_leading_parent_name"),
        pytest.param("CC(=O)NCCCl", "N-(2-chloroethyl)acetamide", id="halogenated_n_substituent"),
        pytest.param("CC(=O)NC=C", "N-ethenylacetamide", id="unsaturated_n_substituent"),
        pytest.param("O=C1CCCN1", "pyrrolidin-2-one", id="lactam_is_named_via_ketone_module"),
        pytest.param("NC(=O)CC(N)=O", "propanediamide", id="diamide"),
        pytest.param("OCc1ccccc1CC(N)=O", "2-[2-(hydroxymethyl)phenyl]acetamide", id="phenyl_chain_amide_with_hydroxyl"),
        pytest.param("C=Cc1ccccc1CC(N)=O", "2-(2-ethenylphenyl)acetamide", id="phenyl_chain_amide_unsaturation"),
        pytest.param("O=CNc1ccccc1", "N-phenylformamide", id="n_phenyl_formamide"),
        pytest.param("CC(=O)N(C)c1ccccc1", "N-methyl-N-phenylacetamide", id="n_phenyl_amide_with_second_n_substituent"),
        pytest.param("CC(=O)Nc1ccc(O)cc1", "N-(4-hydroxyphenyl)acetamide", id="n_phenyl_amide_substituted_ring"),
        pytest.param("CC(=O)NC1=NCCS1", "N-(4,5-dihydro-1,3-thiazol-2-yl)acetamide", id="ring_amidine_carbon_is_a_ring_atom_of_the_n_substituent"),
        pytest.param(
            "CN(C)CCNC1=NCCC1",
            "N2-(3,4-dihydro-2H-pyrrol-5-yl)-N1,N1-dimethylethane-1,2-diamine",
            id="ring_amidine_substituent_of_a_diamine",
        ),
    ],
)
def test_n_n_disubstituted_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzamide():
    assert smiles_to_iupac("NC(=O)c1ccccc1") == "benzamide"
    assert smiles_to_iupac("NC(=O)c1ccc(C)cc1") == "4-methylbenzamide"
    assert smiles_to_iupac("NC(=O)c1ccccc1Cl") == "2-chlorobenzamide"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OC=CC(N)=O", "3-hydroxyprop-2-enamide", id="amide_enol_mix"),
        pytest.param("CC[C@@H](C)C(=O)NCC", "(2R)-N-ethyl-2-methylbutanamide", id="n_substituted_amide_stereocenter"),
        pytest.param("NC(=O)[C@H]1CCCC[C@@H]1Cl", "(1R,2S)-2-chlorocyclohexane-1-carboxamide", id="cyclic_amide_ring_stereocenter"),
        pytest.param("NC(=O)[C@H]1CC[C@H](C[C@@H](C)CC)C1", "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carboxamide", id="amide_ring_branch_stereocenter"),
    ],
)
def test_amide_enol_mix_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)NO", "N-hydroxyacetamide"),
        ("ONC(=O)c1ccc(Cl)cc1", "4-chloro-N-hydroxybenzamide"),
    ],
)
def test_hydroxamic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CN(C)O", "N-hydroxy-N-methylformamide"),
        ("NC(=O)C(N)=O", "oxamide"),
        ("N#CC#N", "oxalonitrile"),
    ],
)
def test_retained_stems_for_one_and_two_carbon_acyl_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1CCC(CC1)C(=O)NO", "N-hydroxycyclohexanecarboxamide", id="on_ring"),
        pytest.param("CC(=O)N(C)O", "N-hydroxy-N-methylacetamide", id="with_n_alkyl"),
    ],
)
def test_hydroxamic_acid_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=Nc1ccccc1)N", "N'-phenylethanimidamide", id="phenylamidine"),
        pytest.param("CC(=NC)Nc1ccccc1", "N'-methyl-N-phenylethanimidamide", id="methyl_n_phenylamidine"),
    ],
)
def test_n_prime_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_amidine_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC(=N)N)cc1") == "3-(4-chlorophenyl)propanimidamide"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NCCN(C)CC=C", "N1-methyl-N1-(prop-2-en-1-yl)ethane-1,2-diamine", id="chain_diamine_with_an_unsaturated_n_substituent"),
        pytest.param(
            "C=CN(C)CCNCCN(C)C",
            "N2-[2-(dimethylamino)ethyl]-N1-ethenyl-N1-methylethane-1,2-diamine",
            id="triamine_with_an_ethenyl_group_and_lowest_n_locants_for_all_substituents",
        ),
        pytest.param("CNCc1ncccc1N", "2-[(methylamino)methyl]pyridin-3-amine", id="secondary_amine_prefix_on_a_pyridine_amine"),
        pytest.param("CNCc1nc(Br)ccc1N", "6-bromo-2-[(methylamino)methyl]pyridin-3-amine", id="halogenated_pyridine_amine_with_a_secondary_amine_prefix"),
        pytest.param("CNc1cccnc1N", "N3-methylpyridine-2,3-diamine", id="secondary_amine_nitrogen_on_the_ring_is_a_principal_group"),
        pytest.param("NCCNCc1ccccn1", "N1-[(pyridin-2-yl)methyl]ethane-1,2-diamine", id="chain_diamine_beats_a_ring_with_one_amine"),
        pytest.param("Nc1ccc2ccccc2c1CNC", "1-[(methylamino)methyl]naphthalen-2-amine", id="fused_ring_amine_with_a_secondary_amine_prefix"),
    ],
)
def test_polyamines_of_a_ring_and_acyclic_amine_nitrogens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCCC(CC)(NCC)NC", "N3-ethyl-N'3-methylhexane-3,3-diamine", id="geminal_primed_nitrogen_locants"),
        pytest.param("CNCCC(NCC)(NC)CCCN", "N3-ethyl-N1,N'3-dimethylhexane-1,3,3,6-tetramine", id="geminal_among_other_amines"),
        pytest.param("CCNC(C)NCC", "N1,N'1-diethylethane-1,1-diamine", id="geminal_same_substituents"),
        pytest.param("NCCNc1ccccc1", "N1-phenylethane-1,2-diamine", id="chain_parent_with_aryl_on_nitrogen"),
        pytest.param("CNc1ccc(N)cc1", "N1-methylbenzene-1,4-diamine", id="ring_parent_n_substituted"),
        pytest.param("CNC1CCC(NC)CC1", "N1,N4-dimethylcyclohexane-1,4-diamine", id="saturated_ring_parent_with_both_nitrogens_substituted"),
        pytest.param("CNC1(NC)CCOCC1", "N4,N'4-dimethyloxane-4,4-diamine", id="geminal_primed_nitrogen_locants_on_a_ring"),
        pytest.param(
            "Nc1ccc(Nc2ccc(Nc3ccccc3)cc2)cc1",
            "N1-(4-aminophenyl)-N4-phenylbenzene-1,4-diamine",
            id="ring_parent_aminophenyl_branch",
        ),
        pytest.param("NCc1ccc(NC)cc1", "4-(aminomethyl)-N-methylaniline", id="ring_with_one_amine_beats_chain"),
        pytest.param(
            "NC1=C(OC2=C(NC)C=CC=C2)C=CC(=C1)C",
            "2-(2-amino-4-methylphenoxy)-N-methylaniline",
            id="parent_with_the_lower_prefix_locant_set_counts_the_n_prefix",
        ),
    ],
)
def test_polyamine_n_locants_and_ring_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCNCNCC", "N,N'-methylenediethanamine", id="methylene_joins_identical_monoamine_nitrogens"),
        pytest.param(
            "Cc1ccc(NCNc2ccc(C)cc2)cc1", "N,N'-methylenebis(4-methylaniline)", id="substituted_monoamine_units_take_bis"
        ),
        pytest.param(
            "CC(=O)N(C)CN(C)C(C)=O", "N,N'-methylenebis(N-methylacetamide)", id="substituted_nitrogens_stay_in_the_unit"
        ),
        pytest.param(
            "OCCNCNCCO", "2,2'-[methylenebis(azanediyl)]di(ethan-1-ol)", id="nitrogen_cited_as_prefix_is_not_a_parent_locant"
        ),
        pytest.param(
            "CCN(C)CNCC", "N,N'-diethyl-N-methylmethanediamine", id="geminal_diamine_unprimed_nitrogen_has_more_substituents"
        ),
        pytest.param(
            "C=CCNCNC(C)C",
            "N-(propan-2-yl)-N'-(prop-2-en-1-yl)methanediamine",
            id="geminal_diamine_unsaturated_substituent_cited_in_alphanumerical_order",
        ),
        pytest.param(
            "C1CCCCC1NCNCC", "N-cyclohexyl-N'-ethylmethanediamine", id="geminal_diamine_ring_substituent_on_nitrogen"
        ),
        pytest.param("NCCNCNCCN", "N1,N1'-methylenedi(ethane-1,2-diamine)", id="methylene_joins_diamine_nitrogens"),
        pytest.param(
            "Nc1ccc(NCNc2ccc(N)cc2)cc1", "N1,N1'-methylenedi(benzene-1,4-diamine)", id="methylene_joins_aryl_diamine"
        ),
        pytest.param("CC(=O)NNCNNC(C)=O", "N',N'''-methylenediacetohydrazide", id="methylene_joins_hydrazide_nitrogens"),
        pytest.param(
            "CC(=O)NNCCNNC(C)=O", "N',N'''-(ethane-1,2-diyl)diacetohydrazide", id="longer_chain_joins_hydrazide_nitrogens"
        ),
        pytest.param("CNCCNC", "N1,N2-dimethylethane-1,2-diamine", id="longer_chain_between_amines_is_a_diamine"),
        pytest.param("NCCNCCN", "N1-(2-aminoethyl)ethane-1,2-diamine", id="azanediyl_between_amine_units_is_substitutive"),
        pytest.param("Nc1ccc(Nc2ccc(N)cc2)cc1", "N1-(4-aminophenyl)benzene-1,4-diamine", id="azanediyl_between_aniline_units"),
    ],
)
def test_nitrogen_multiplicative_amines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCNc1ccc(Cl)cc1", "4-chloro-N-ethylaniline"),
        ("c1ccccc1CNc1ccccc1", "N-benzylaniline"),
        ("Brc1cc(Cl)ccc1Nc1ccc(Br)cc1Br", "2-bromo-4-chloro-N-(2,4-dibromophenyl)aniline"),
        ("Brc1ccc(N(C)c2ccc(Cl)cc2Br)c(Br)c1", "2-bromo-4-chloro-N-(2,4-dibromophenyl)-N-methylaniline"),
        ("Brc1ccc(Cl)cc1Nc1ccc(Br)cc1Br", "2,4-dibromo-N-(2-bromo-5-chlorophenyl)aniline"),
        ("CCN(CC)C1=NCCN1", "N,N-diethyl-4,5-dihydro-1H-imidazol-2-amine"),
        ("ONC1=NCCCO1", "N-hydroxy-5,6-dihydro-4H-1,3-oxazin-2-amine"),
    ],
)
def test_n_substituted_aniline_and_multi_ring_chain_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC=CC", "prop-1-en-1-amine", id="enamine"),
        pytest.param("CC(N)=CC(C)=NC", "4-(methylimino)pent-2-en-2-amine", id="enamine_outranks_an_imine_cited_as_a_prefix"),
        pytest.param("C/C=C\\N(/C=C\\C)C(C)(C)C", "2-methyl-N,N-bis[(1Z)-prop-1-en-1-yl]propan-2-amine", id="stereodefined_substituents_outside_the_parent_arm"),
        pytest.param("C=C(C)N(/C=C/N)CC(C)C", "(1E)-N1-(2-methylpropyl)-N1-(prop-1-en-2-yl)ethene-1,2-diamine", id="stereodefined_bond_inside_the_parent_arm"),
        pytest.param("NC1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexan-1-amine", id="cyclic_amine_branch_stereocenter"),
        # P-44.1.1, P-62.2.4.1.2, P-14.3.4.4: every amine group of the chain is in the suffix, also on a carbon of the double bond
        pytest.param("NC(N)=C", "ethene-1,1-diamine", id="geminal_groups_on_the_double_bond_of_ethene"),
        pytest.param("NC(N)=CC", "prop-1-ene-1,1-diamine", id="geminal_groups_on_the_double_bond_of_a_longer_chain"),
        pytest.param("NC(N)=C1CCCC1", "cyclopentylidenemethanediamine", id="geminal_groups_on_the_carbon_of_a_ylidene"),
        pytest.param("ClC=CN", "2-chloroethen-1-amine", id="locant_of_the_double_bond_of_a_substituted_ethene_omitted"),
        pytest.param("C[Si](C)(C)C/C=C/N(C)C", "(1E)-N,N-dimethyl-3-(trimethylsilyl)prop-1-en-1-amine", id="silyl_group_beside_a_configured_tertiary_enamine"),
        pytest.param("C[Si](C)(C)C/C=C(/C)N(C)C", "(2Z)-N,N-dimethyl-4-(trimethylsilyl)but-2-en-2-amine", id="silyl_group_beside_a_configured_enamine_on_a_branched_carbon"),
    ],
)
def test_enamine_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC=C[C@H]1CCCN1", "(2R)-2-(prop-1-en-1-yl)pyrrolidine", id="nitrogen_ring"),
        pytest.param("CC=C[C@H]1CCCO1", "(2R)-2-(prop-1-en-1-yl)oxolane", id="oxygen_ring"),
        pytest.param("CC=C[C@H](C)N1CCCC1", "1-[(2S)-pent-3-en-2-yl]pyrrolidine", id="centre_on_a_chain_bonded_to_a_ring_nitrogen"),
    ],
)
def test_unspecified_double_bond_adds_no_descriptor_beside_a_ring_stereocentre(smiles, expected):
    # P-93.4.2.1.1
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_primary_amine_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("C/C=C/[C@@H](C)N") == "(2R,3E)-pent-3-en-2-amine"
    assert smiles_to_iupac("C/C=C/[C@H](C)N") == "(2S,3E)-pent-3-en-2-amine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[N+](C)(C)[C@H](C)/C=C/C", "(2R,3E)-N,N,N-trimethylpent-3-en-2-aminium", id="ammonium_stereocenter_with_ez_double_bond_coexistence"),
        pytest.param("C[C@H](CC(=O)OC)[N+](C)(C)C", "(2R)-4-methoxy-N,N,N-trimethyl-4-oxobutan-2-aminium", id="quaternary_aminium_stereocentre_beside_an_ester"),
        pytest.param("C=Cc1ccccc1CCN", "2-(2-ethenylphenyl)ethan-1-amine", id="phenyl_chain_amine_unsaturation"),
        pytest.param("c1ccccc1C(N)CCN", "1-phenylpropane-1,3-diamine", id="phenyl_chain_diamine"),
    ],
)
def test_ammonium_stereocenter_with_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NCC1CCCCC1", "1-cyclohexylmethanamine", id="substituent_chain_amine"),
        pytest.param("NCNN", "1-hydrazinylmethanamine", id="locant_kept_beside_a_prefix_that_could_stand_on_nitrogen"),
        pytest.param("NC(Cl)Cl", "1,1-dichloromethanamine", id="every_prefix_locant_of_a_substituted_methanamine"),
        pytest.param("NC1CCCCC1CN", "2-(aminomethyl)cyclohexan-1-amine", id="with_amine_chain_amine_tie"),
        pytest.param("NC1CCCCC1C(N)CN", "1-(2-aminocyclohexyl)ethane-1,2-diamine", id="with_amine_chain_amine_chain_wins"),
    ],
)
def test_ring_substituent_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_amine_oxide_ring_nitrogen():
    assert smiles_to_iupac("C[N+]1([O-])CCCCC1") == "1-methylpiperidine 1-oxide"


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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Cc1ccccc1CN=[N+]=[N-]", "1-(azidomethyl)-2-methylbenzene"),
        ("C1CCC(N=[N+]=[N-])CC1", "azidocyclohexane"),
        ("C=CCN=[N+]=[N-]", "3-azidoprop-1-ene"),
    ],
)
def test_azide_prefix_on_substituted_ring_and_unsaturated_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)=NN=C(C)C", "di(propan-2-ylidene)hydrazine", id="symmetric_azine"),
        pytest.param("CCC(C)=NN=C1CCCCC1", "(butan-2-ylidene)(cyclohexylidene)hydrazine", id="unsymmetrical_azine_with_a_ring"),
        pytest.param("OC(=O)C1CCC(CC1)=NN=C(C)C", "4-[(propan-2-ylidene)hydrazinylidene]cyclohexane-1-carboxylic acid", id="azine_as_a_prefix_on_an_acid"),
        pytest.param("NC(=S)NN=C(C)C", "2-(propan-2-ylidene)hydrazine-1-carbothioamide", id="thiosemicarbazone"),
        pytest.param("NNC(N)=[Se]", "hydrazinecarboselenoamide", id="selenosemicarbazide"),
        pytest.param("CC[C@@H](C)C(C)=NN=C(C)[C@H](C)CC", "bis[(3R)-3-methylpentan-2-ylidene]hydrazine", id="azine_with_specified_chain_stereocenters"),
    ],
)
def test_azines_and_chalcogen_semicarbazones_are_hydrazine_derivatives(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_methyl_n_n_dimethylcarbamate():
    assert smiles_to_iupac("COC(=O)N(C)C") == "methyl dimethylcarbamate"


def test_plain_cyclic_r_group():
    assert smiles_to_iupac("O=C(N)OC1CCCCC1") == "cyclohexyl carbamate"
    assert smiles_to_iupac("NC(=O)OC1CCCC1") == "cyclopentyl carbamate"


def test_n_methyl_cyclic_r_group():
    assert smiles_to_iupac("CNC(=O)OC1CCCCC1") == "cyclohexyl methylcarbamate"


def test_phenyl_cyanate():
    assert smiles_to_iupac("c1ccccc1OC#N") == "phenyl cyanate"  # CID 70740


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Cc1ccccc1N=N", "(2-methylphenyl)diazene", id="substituted_phenyl"),
        pytest.param("ClCCN=NCCCl", "bis(2-chloroethyl)diazene", id="bis_chloroethyldiazene_name"),
    ],
)
def test_substituted_phenyl_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_on_diazene_nitrogen():
    assert smiles_to_iupac("ClN=N") == "chlorodiazene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=[N+]=[N-]", "diazomethane"),
    ],
)
def test_diazo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1(=[N+]=[N-])CCCCC1", "diazocyclohexane", id="diazo_on_a_ring"),
        pytest.param("C=CC=[N+]=[N-]", "3-diazoprop-1-ene", id="diazo_beside_a_double_bond"),
        pytest.param("[N-]=[N+]=CC=[N+]=[N-]", "1,2-didiazoethane", id="two_diazo_groups"),
        pytest.param("[N-]=[N+]=CC(=O)O", "diazoacetic acid", id="diazo_beside_a_carboxylic_acid"),
        pytest.param("[N-]=[N+]=CC(=O)OCC", "ethyl diazoacetate", id="diazo_acetate_ester"),
        pytest.param("[N-]=[N+]=CC(C)=O", "1-diazopropan-2-one", id="diazo_beside_a_ketone"),
        pytest.param("[N-]=[N+]=C(C(C)=O)[Si](C)(C)C", "1-diazo-1-(trimethylsilyl)propan-2-one", id="diazo_and_silyl_on_one_carbon"),
        pytest.param("[N-]=[N+]=C1C=CC=C1", "5-diazocyclopenta-1,3-diene", id="diazo_on_an_unsaturated_ring"),
        pytest.param("[N-]=[N+]=C1C=CC(C(O)=O)=C1", "3-diazocyclopenta-1,4-diene-1-carboxylic acid", id="diazo_on_an_unsaturated_ring_with_acid"),
    ],
)
def test_diazo_prefix_on_rings_and_beside_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Si](C)(C)[N+](=O)[O-]", "trimethyl(nitro)silane", id="nitro_on_silicon"),
        pytest.param("CB[N+](=O)[O-]", "methyl(nitro)borane", id="nitro_on_boron"),
        pytest.param("C[Si](C)(C)N=O", "trimethyl(nitroso)silane", id="nitroso_on_silicon"),
        pytest.param("C[Si](C)(C)N=[N+]=[N-]", "azidotri(methyl)silane", id="azido_on_silicon"),
        pytest.param("C[Ge](C)(C)[N+](=O)[O-]", "trimethyl(nitro)germane", id="nitro_on_germanium"),
        pytest.param("CB(C)N=O", "dimethyl(nitroso)borane", id="nitroso_on_boron"),
    ],
)
def test_nitrogen_group_prefixes_on_mononuclear_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[N+]#N", "methanediazonium", id="methanediazonium"),
        pytest.param("ClCC[N+]#N", "2-chloroethane-1-diazonium", id="chloroethane_diazonium"),
        pytest.param("C=CCCC[N+]#N", "pent-4-ene-1-diazonium", id="pent_4_ene_1_diazonium"),
        pytest.param("N#[N+]Cc1ccc(N)cc1", "1-(4-aminophenyl)methanediazonium", id="amino_group_on_the_substituent_skeleton"),
        pytest.param("N#[N+]c1ccc(CN)cc1", "4-(aminomethyl)benzenediazonium", id="amino_group_on_a_chain_of_the_diazonium_ring"),
        pytest.param("N#[N+]CC(N)CC(N)C", "2,4-diaminopentane-1-diazonium", id="amino_groups_on_the_diazonium_skeleton"),
    ],
)
def test_methanediazonium_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_diazonium_groups_are_bis_diazonium():
    assert smiles_to_iupac("N#[N+]CC[N+]#N") == "ethane-1,2-bis(diazonium)"
    assert smiles_to_iupac("N#[N+]c1ccc(cc1)[N+]#N") == "benzene-1,4-bis(diazonium)"


def test_cyclohexanediazonium():
    assert smiles_to_iupac("C1CCCCC1[N+]#N") == "cyclohexanediazonium"


def test_substituted_ring_diazonium():
    assert smiles_to_iupac("CC1CCCCC1[N+]#N") == "2-methylcyclohexane-1-diazonium"


def test_benzenediazonium():
    assert smiles_to_iupac("c1ccccc1[N+]#N") == "benzenediazonium"  # CID 9718


def test_diazonium_outranks_alcohol_and_ketone():
    assert smiles_to_iupac("OCC[N+]#N") == "2-hydroxyethane-1-diazonium"
    assert smiles_to_iupac("CC(=O)C(C(C)=O)[N+]#N") == "2,4-dioxopentane-3-diazonium"


def test_phenyl_chain_diazonium():
    assert smiles_to_iupac("c1ccccc1CC[N+]#N") == "2-phenylethane-1-diazonium"
    assert smiles_to_iupac("c1ccccc1CCC[N+]#N") == "3-phenylpropane-1-diazonium"


def test_dimethyl_benzenediazonium():
    assert smiles_to_iupac("Cc1ccc(C)cc1[N+]#N") == "2,5-dimethylbenzenediazonium"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("Cc1ccccc1CC[N+]#N", "2-(2-methylphenyl)ethane-1-diazonium", id="substituted_benzene_ring_diazonium"),
        pytest.param("C=Cc1ccccc1CC[N+]#N", "2-(2-ethenylphenyl)ethane-1-diazonium", id="chain_diazonium_beside_alkenyl_ring"),
    ],
)
def test_phenyl_substituted_benzene_diazonium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_guanidine():
    assert smiles_to_iupac("NC(=N)N") == "guanidine"


def test_imino_plus_amino_different_names_alphabetized():
    assert smiles_to_iupac("CNC(=NCC)N") == "N''-ethyl-N-methylguanidine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1NC(=N)Nc1ccccc1", "N,N'-diphenylguanidine", id="n_prime_diphenylguanidine"),
        pytest.param("NC(=Nc1ccccc1)N", "N''-phenylguanidine", id="double_prime_phenylguanidine"),
    ],
)
def test_n_n_prime_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_pyrrolidine_amine_2():
    assert smiles_to_iupac("NC1CCCN1") == "pyrrolidin-2-amine"


def test_secondary_exocyclic_amine():
    assert smiles_to_iupac("CNC1CCNCC1") == "N-methylpiperidin-4-amine"


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C(=O)NN", "formohydrazide", id="formohydrazide"),
        pytest.param("CC(=O)NN", "acetohydrazide", id="acetohydrazide"),
        pytest.param("NNC(=O)C(=O)NN", "oxalohydrazide", id="oxalohydrazide"),
        pytest.param("CNNC(=O)C(=O)NNC", "N'1,N'2-dimethyloxalohydrazide", id="n_substituted_oxalohydrazide"),
        pytest.param("ClC(C(=O)N(N(C)C)C)C", "2-chloro-N,N',N'-trimethylpropanehydrazide", id="n_locants_ordered_with_other_prefixes"),
        pytest.param("CC[C@@H](C)C(=O)N(C)N", "(2R)-N,2-dimethylbutanehydrazide", id="n_substituted_stereocenter_hydrazide"),
        pytest.param("CN(C)NC=O", "N',N'-dimethylformohydrazide", id="n_substituted_formohydrazide"),
        pytest.param("O=CNN=CCC", "N'-propylideneformohydrazide", id="ylidene_on_the_terminal_nitrogen_of_formohydrazide"),
        pytest.param("CC(C)=NN(C)C(=S)CC", "N-methyl-N'-(propan-2-ylidene)propanethiohydrazide", id="ylidene_on_a_thiohydrazide"),
    ],
)
def test_formohydrazide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_hydrazide_retained_name():
    assert smiles_to_iupac("c1ccccc1CC(=O)NN") == "2-phenylacetohydrazide"


def test_n_prime_phenylhydrazide():
    assert smiles_to_iupac("CC(=O)NNc1ccccc1") == "N'-phenylacetohydrazide"
    assert smiles_to_iupac("CCC(=O)NNc1ccccc1") == "N'-phenylpropanehydrazide"


def test_n_phenylhydrazide():
    assert smiles_to_iupac("CC(=O)N(c1ccccc1)N") == "N-phenylacetohydrazide"


def test_phenyl_chain_hydrazide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC(=O)NN)cc1") == "3-(4-chlorophenyl)propanehydrazide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)NNC(=O)C", "N'-acetylacetohydrazide"),
        ("CCC(=O)NNC(=O)CC", "N'-propanoylpropanehydrazide"),
    ],
)
def test_symmetric_diacyl_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)NNC(=O)NC", "N'-(methylcarbamoyl)acetohydrazide", id="carbamoyl_on_terminal_hydrazide_nitrogen"),
        pytest.param("c1ccccc1C(=O)NNC(=S)Nc1ccccc1", "N'-(phenylcarbamothioyl)benzohydrazide", id="carbamothioyl_keeps_the_carboxylic_hydrazide_parent"),
        pytest.param("CC(=S)NNC(N)=O", "N'-carbamoylethanethiohydrazide", id="thiohydrazide_with_carbamoyl"),
    ],
)
def test_hydrazide_with_carbamoyl_or_carbamothioyl_on_the_terminal_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)N(C)C#N", "N-cyano-N-methylacetamide", id="cyano_on_amide_nitrogen"),
        pytest.param("N#CN(N)C(=O)CO", "N-cyano-2-hydroxyacetohydrazide", id="cyano_on_hydrazide_nitrogen"),
        pytest.param("COC(=O)NC#N", "methyl cyanocarbamate", id="cyano_on_carbamate_nitrogen"),
        pytest.param("CSC(=Nc1ccc(F)cc1)NC#N", "methyl N-cyano-N'-(4-fluorophenyl)carbamimidothioate", id="cyano_on_amino_nitrogen_of_carbamimidothioate"),
        pytest.param("N#CN=C(N)SC", "methyl N'-cyanocarbamimidothioate", id="cyano_on_imino_nitrogen_of_carbamimidothioate"),
        pytest.param("N#CNC(=N)N", "N-cyanoguanidine", id="cyanoguanidine_unchanged"),
    ],
)
def test_cyano_on_nitrogen_is_a_prefix_of_the_senior_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)NNCCNN", "N'-(2-hydrazinylethyl)acetohydrazide", id="hydrazine_at_the_chain_end"),
        pytest.param("CCC(=O)N(CNN)N", "N-(hydrazinylmethyl)propanehydrazide", id="hydrazine_on_the_carbon_bonded_to_the_hydrazide_nitrogen"),
    ],
)
def test_hydrazide_nitrogen_substituent_keeps_its_own_hydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diacyl_hydrazide_formyl():
    # Mononuclear retained name on both sides.
    assert smiles_to_iupac("O=CNNC=O") == "N'-formylformohydrazide"


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("NN1CCCCC1", "piperidin-1-amine"),
        ("NN1CCCC1", "pyrrolidin-1-amine"),
        ("NN1CCCCC1C", "2-methylpiperidin-1-amine"),
        ("NN1CCOCC1", "morpholin-4-amine"),
        ("NN1C=CC=C1", "1H-pyrrol-1-amine"),
        ("CCNn1cccc1", "N-ethyl-1H-pyrrol-1-amine"),
        ("CCCC(C)CNn1c(C)ccc1C", "2,5-dimethyl-N-(2-methylpentyl)-1H-pyrrol-1-amine"),
        ("CCNn1ccnc1", "N-ethyl-1H-imidazol-1-amine"),
        ("CN(C)N1CCCC1", "N,N-dimethylpyrrolidin-1-amine"),
        ("c1ccncc1Nn1cccc1", "N-(1H-pyrrol-1-yl)pyridin-3-amine"),
        ("CCNn1cccc1C(=O)O", "1-(ethylamino)-1H-pyrrole-2-carboxylic acid"),
    ],
)
def test_ring_nitrogen_amino_is_amine_not_hydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("ClNN", "chlorohydrazine"),
        ("F[N](F)N(F)F", "tetrafluorohydrazine"),
        ("CN(C)N(C)C", "tetramethylhydrazine"),
        ("CN(C)N(C)CC", "1-ethyl-1,2,2-trimethylhydrazine"),
        ("c1ccccc1NC(=O)NN", "N-phenylhydrazinecarboxamide"),
    ],
)
def test_halogens_on_hydrazine_nitrogen_and_fully_substituted_hydrazines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hydrazone_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)C(C)=NN") == "3-methylpentan-2-ylidenehydrazine"


def test_hydrazone_specified_chain_stereocenter_is_cited():
    assert smiles_to_iupac("CC[C@@H](C)C(C)=NN") == "[(3R)-3-methylpentan-2-ylidene]hydrazine"


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


def test_phenyl_o_substituent():
    assert smiles_to_iupac("NOc1ccccc1") == "O-phenylhydroxylamine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)C(=O)NC(=O)C(C)C", "2-methyl-N-(2-methylpropanoyl)propanamide", id="branched_symmetric_imide"),
        pytest.param("ClCC(=O)NC(=O)CCl", "2-chloro-N-(chloroacetyl)acetamide", id="halogen_substituent"),
        pytest.param("CC(=O)NC(=O)C", "N-acetylacetamide", id="diacetylamine"),
        pytest.param("CC(=O)N(C(C)=O)C(C)=O", "N,N-diacetylacetamide", id="triacetylamine"),
        pytest.param("O=CNC=O", "N-formylformamide", id="diformylamine"),
        pytest.param("CC(=O)NC(=O)CC", "N-acetylpropanamide", id="unsymmetric_imide"),
        pytest.param("C=CC(=O)NC(=O)C=C", "N-(prop-2-enoyl)prop-2-enamide", id="unsaturated_imide"),
        pytest.param("O=C(c1ccccc1)NC(=O)c1ccccc1", "N-benzoylbenzamide", id="dibenzoylamine"),
        pytest.param("O=C(c1ccccc1)NC(C)=O", "N-acetylbenzamide", id="acetyl_benzoyl_amine"),
        pytest.param("O=C(c1ccco1)NC(=O)c1ccco1", "N-(furan-2-carbonyl)furan-2-carboxamide", id="difuroylamine"),
        pytest.param(
            "O=C(C1CCCCC1)N(C(=O)C1CCCCC1)C(=O)C1CCCCC1",
            "N,N-di(cyclohexanecarbonyl)cyclohexanecarboxamide",
            id="tri_acylamine_ring_acyl",
        ),
        pytest.param("CC(=O)N(C(=O)CCCl)C(=O)c1ccccc1", "N-acetyl-N-(3-chloropropanoyl)benzamide", id="three_different_acyls"),
        pytest.param("CC(=O)N(C1CCCC1)C(C)=O", "N-acetyl-N-cyclopentylacetamide", id="acylamine_with_alkyl"),
        pytest.param("CC(=O)N(c1ccc2ccccc2c1)C(=O)c1ccccc1", "N-acetyl-N-(naphthalen-2-yl)benzamide", id="acylamine_with_aryl"),
        pytest.param("CS(=O)(=O)NS(C)(=O)=O", "N-(methanesulfonyl)methanesulfonamide", id="disulfonylamine"),
        pytest.param("CS(=O)(=O)NC(C)=O", "N-(methanesulfonyl)acetamide", id="carboxamide_senior_to_sulfonamide"),
        pytest.param(
            "O=C(c1ccc2ccccc2c1)N(C(=O)c1ccc2ncccc2c1)C",
            "N-methyl-N-(naphthalene-2-carbonyl)quinoline-6-carboxamide",
            id="heterocyclic_ring_system_before_carbocyclic",
        ),
        pytest.param(
            "O=C(c1ccsc1)N(C(=O)c1ccoc1)C",
            "N-methyl-N-(thiophene-3-carbonyl)furan-3-carboxamide",
            id="oxygen_ring_before_sulfur_ring",
        ),
    ],
)
def test_branched_symmetric_imide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)CCN(C(C)=O)C(C)=O", "3-(N-acetylacetamido)propanoic acid", id="diacylamino_prefix"),
        pytest.param(
            "OC(=O)CCN(S(C)(=O)=O)S(C)(=O)=O",
            "3-[N-(methanesulfonyl)methanesulfonamido]propanoic acid",
            id="disulfonylamino_prefix",
        ),
    ],
)
def test_diacylamino_prefix_is_the_n_acylamido_of_the_senior_acyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_cyclic_imide_named_via_ketone_module():
    assert smiles_to_iupac("CN1C(=O)CCC1=O") == "1-methylpyrrolidine-2,5-dione"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC=NOCC", "N-ethoxypropan-1-imine"),
    ],
)
def test_smiles_to_iupac_simple_imine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_raises_imine_is_named():
    assert smiles_to_iupac("N=C1CCCC1") == "cyclopentanimine"


def test_imine_with_a_senior_alcohol_is_an_imino_prefix():
    assert smiles_to_iupac("OCC=N") == "2-iminoethan-1-ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("ON=C1CCc2ccccc21", "N-hydroxy-2,3-dihydro-1H-inden-1-imine", id="oxime_on_a_fused_ring"),
        pytest.param(
            "ON=C1CC(c2cc(Cl)cc(Cl)c2)c2ccccc21",
            "3-(3,5-dichlorophenyl)-N-hydroxy-2,3-dihydro-1H-inden-1-imine",
            id="n_prefix_between_ring_prefixes",
        ),
        pytest.param("N=C1CCCN1CCOC", "1-(2-methoxyethyl)pyrrolidin-2-imine", id="heterocyclic_imine_with_an_ether_on_the_ring_nitrogen"),
        pytest.param("COc1ccc(N2CCCC2=N)c(C)n1", "1-(6-methoxy-2-methylpyridin-3-yl)pyrrolidin-2-imine", id="heterocyclic_imine_with_a_heteroaryl_on_the_ring_nitrogen"),
    ],
)
def test_imines_on_fused_and_heterocyclic_rings_with_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC=NOO", id="two_oxygens_on_oxime_nitrogen_raises"),
    ],
)
def test_branched_oxime_o_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=N/C=N\\N(C)C", "N-[(Z)-N-(dimethylamino)methanimidoyl]methanimine"),
    ],
)
def test_imidoyl_prefix_with_a_disubstituted_amino_nitrogen_cites_its_descriptor(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)/C(C#N)=N/O", "(1E)-N-hydroxy-2-oxopropanimidoyl cyanide"),
        ("CC(=O)/C(C#N)=N\\O", "(1Z)-N-hydroxy-2-oxopropanimidoyl cyanide"),
        ("O=C/C(C#N)=N/O", "(1E)-N-hydroxy-2-oxoethanimidoyl cyanide"),
        ("O=CC(=NO)O", "N-hydroxy-2-oxoethanimidic acid"),
    ],
)
def test_imidoyl_cyanide_descriptor_ranks_the_acyl_ligand_above_cyano(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_imine():
    assert smiles_to_iupac("c1ccccc1CC=N") == "2-phenylethan-1-imine"
    assert smiles_to_iupac("c1ccccc1CCC=N") == "3-phenylpropan-1-imine"
    assert smiles_to_iupac("c1ccccc1CC(C)=N") == "1-phenylpropan-2-imine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)N=C=O", "2-isocyanatopropane"),
    ],
)
def test_isocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC(N=C=O)CC1", "isocyanatocyclohexane"),
        ("C=CCN=C=O", "3-isocyanatoprop-1-ene"),
        ("O=C=NCN=C=O", "diisocyanatomethane"),
        pytest.param(
            "O=C=NC[C@@H]1CCCCO1", "(2S)-2-(isocyanatomethyl)oxane", id="cumulated_bond_is_not_an_unspecified_stereo_element"
        ),
    ],
)
def test_isocyanate_on_ring_unsaturated_and_repeated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


def test_phenyl_isocyanide_chain():
    assert smiles_to_iupac("c1ccccc1C[N+]#[C-]") == "(isocyanomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CC[N+]#[C-]") == "(2-isocyanoethyl)benzene"


def test_nitro_unaffected():
    assert smiles_to_iupac("C[N+](=O)[O-]") == "nitromethane"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("N#CC(CC#N)CC#N", "propane-1,2,3-tricarbonitrile", id="trinitrile"),
        pytest.param("N#CC1CCC(C#N)CC1", "cyclohexane-1,4-dicarbonitrile", id="dinitrile_alongside_ring"),
    ],
)
def test_trinitrile_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzonitrile():
    assert smiles_to_iupac("N#Cc1ccccc1") == "benzonitrile"
    assert smiles_to_iupac("N#Cc1ccc(C)cc1") == "4-methylbenzonitrile"
    assert smiles_to_iupac("N#Cc1ccccc1Cl") == "2-chlorobenzonitrile"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("N#Cc1ccc(Cc2ccccc2)c(C#N)c1", "4-benzylbenzene-1,3-dicarbonitrile", id="benzyl"),
        pytest.param("N#Cc1ccc(C#N)c(Cc2ccccc2)c1", "2-benzylbenzene-1,4-dicarbonitrile", id="locant_set_before_prefix"),
        pytest.param("N#Cc1ccc(C=Cc2ccccc2)c(C#N)c1", "4-(2-phenylethen-1-yl)benzene-1,3-dicarbonitrile", id="styryl"),
        pytest.param("N#Cc1ccc(Cc2ccncc2)c(C#N)c1", "4-[(pyridin-4-yl)methyl]benzene-1,3-dicarbonitrile", id="heteroarylmethyl"),
        pytest.param("N#Cc1ccc(C#N)c(C#N)c1Cc1ccccc1", "3-benzylbenzene-1,2,4-tricarbonitrile", id="trinitrile"),
    ],
)
def test_benzene_with_several_nitriles_cites_an_aromatic_ring_on_a_chain_as_a_prefix(smiles, expected):
    # P-66.5.1.1, P-44.1.1, P-31.1.4
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_nitrile_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC#N") == "(2-ethenylphenyl)acetonitrile"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N#CC1CCCCC1c1cccs1", "2-(thiophen-2-yl)cyclohexane-1-carbonitrile"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_nitrile(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("N#CC1CCCC(c2ccc(C)cc2)C1", "3-(4-methylphenyl)cyclohexane-1-carbonitrile", id="two_ring_aromatic_substituent_nitrile_substituted_ring"),
        pytest.param("N#CCC1CCCCC1c1ccccc1", "(2-phenylcyclohexyl)acetonitrile", id="two_ring_aromatic_substituent_nitrile_chain_nitrile"),
        pytest.param("N#Cc1cccnc1", "pyridine-3-carbonitrile", id="heteroaromatic_ring_directly_attached_nitrile"),
    ],
)
def test_two_ring_aromatic_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_nitrile_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("N#C[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enenitrile"
    assert smiles_to_iupac("N#C[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enenitrile"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("N#C[C@H]1CCCC[C@@H]1Cl", "(1R,2S)-2-chlorocyclohexane-1-carbonitrile", id="cyclic_nitrile_ring_stereocenter"),
        pytest.param("N#C[C@H]1CC[C@H](C[C@@H](C)CC)C1", "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carbonitrile", id="nitrile_ring_branch_stereocenter"),
    ],
)
def test_cyclic_nitrile_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCON=O", "ethyl nitrite"),  # CID 8026
    ],
)
def test_nitrite_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-65.6.3.2.2, P-67.1.3.2: identical nitrate or nitrite groups on one organyl group are multiplied like acetate groups
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=[N+]([O-])OCCO[N+](=O)[O-]", "ethane-1,2-diyl dinitrate", id="dinitrate"),
        pytest.param("CC(CO[N+](=O)[O-])O[N+](=O)[O-]", "propane-1,2-diyl dinitrate", id="vicinal_dinitrate"),
        pytest.param("O=[N+]([O-])OCC(CO[N+](=O)[O-])O[N+](=O)[O-]", "propane-1,2,3-triyl trinitrate", id="trinitrate"),
        pytest.param(
            "O=[N+]([O-])OCC(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]",
            "2,2-bis[(nitrooxy)methyl]propane-1,3-diyl dinitrate",
            id="fourth_nitrate_group_is_a_nitrooxy_prefix",
        ),
        pytest.param("O=[N+]([O-])OC1CCCCC1O[N+](=O)[O-]", "cyclohexane-1,2-diyl dinitrate", id="ring_diyl"),
        pytest.param("O=[N+]([O-])Oc1ccccc1O[N+](=O)[O-]", "1,2-phenylene dinitrate", id="arylene"),
        pytest.param("OCC(CO[N+](=O)[O-])O[N+](=O)[O-]", "3-hydroxypropane-1,2-diyl dinitrate", id="alcohol_beside_the_dinitrate"),
        pytest.param("O=NOCCCCCCON=O", "hexane-1,6-diyl dinitrite", id="dinitrite"),
        pytest.param("O=NOCC(O)CON=O", "2-hydroxypropane-1,3-diyl dinitrite", id="alcohol_beside_the_dinitrite"),
        pytest.param("O=NOC[C@H](C)ON=O", "(2S)-propane-1,2-diyl dinitrite", id="stereocentre_of_the_organyl_group"),
    ],
)
def test_identical_nitrate_and_nitrite_groups_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nitrate_beside_a_nitrite_is_not_multiplied():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[N+]([O-])OCCON=O")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CO[N+](=O)[O-]", "methyl nitrate", id="nitrate_ester_unaffected"),
        pytest.param("C1CCC([N+](=O)[O-])CC1", "nitrocyclohexane", id="ring"),
        pytest.param("C=CC[N+](=O)[O-]", "3-nitroprop-1-ene", id="unsaturated_chain"),
    ],
)
def test_nitrate_ester_unaffected_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-41, P-67.1.3.2, P-65.6.3.2: a nitrate or nitrite ester outranks an alcohol, amine, nitrile, aldehyde, ketone, amide, ether and ring,
# which are prefixes of the organyl group; an acid or another ester outranks it
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=[N+]([O-])OCCO", "2-hydroxyethyl nitrate", id="alcohol"),
        pytest.param("O=[N+]([O-])OCCCO", "3-hydroxypropyl nitrate", id="alcohol_three_carbons"),
        pytest.param("O=[N+]([O-])OC(C)CO", "1-hydroxypropan-2-yl nitrate", id="branched_alcohol_group"),
        pytest.param("O=[N+]([O-])OCCN", "2-aminoethyl nitrate", id="primary_amine"),
        pytest.param("O=[N+]([O-])OCCNC", "2-(methylamino)ethyl nitrate", id="secondary_amine"),
        pytest.param("O=[N+]([O-])OCCN1CCCC1", "2-(pyrrolidin-1-yl)ethyl nitrate", id="ring_amine"),
        pytest.param("O=[N+]([O-])OCCN1CCOCC1", "2-(morpholin-4-yl)ethyl nitrate", id="morpholine"),
        pytest.param("O=[N+]([O-])OCCC#N", "2-cyanoethyl nitrate", id="nitrile"),
        pytest.param("O=[N+]([O-])OCCC=O", "3-oxopropyl nitrate", id="aldehyde"),
        pytest.param("O=[N+]([O-])OCCC(C)=O", "3-oxobutyl nitrate", id="ketone"),
        pytest.param("O=[N+]([O-])OCCOC", "2-methoxyethyl nitrate", id="ether"),
        pytest.param("O=[N+]([O-])OCCNC(=O)c1ccccc1", "2-benzamidoethyl nitrate", id="amide_ranks_below_the_ester"),
        pytest.param("O=[N+]([O-])OCc1ccccn1", "(pyridin-2-yl)methyl nitrate", id="pyridine"),
        pytest.param("O=[N+]([O-])OCc1ccncc1", "(pyridin-4-yl)methyl nitrate", id="pyridine_position_4"),
        pytest.param("O=[N+]([O-])OCCn1ccnc1", "2-(1H-imidazol-1-yl)ethyl nitrate", id="imidazole_through_nitrogen"),
        pytest.param("O=[N+]([O-])OCc1nccs1", "(1,3-thiazol-2-yl)methyl nitrate", id="thiazole"),
        pytest.param("O=[N+]([O-])OCCc1cccs1", "2-(thiophen-2-yl)ethyl nitrate", id="thiophene"),
        pytest.param("O=[N+]([O-])OCCc1ccco1", "2-(furan-2-yl)ethyl nitrate", id="furan"),
        pytest.param("O=[N+]([O-])OCc1ccc2ccccc2n1", "(quinolin-2-yl)methyl nitrate", id="fused_nitrogen_ring"),
        pytest.param("OCCON=O", "2-hydroxyethyl nitrite", id="nitrite_with_alcohol"),
        pytest.param("NCCON=O", "2-aminoethyl nitrite", id="nitrite_with_amine"),
        pytest.param("O=[N+]([O-])OCC(=O)O", "(nitrooxy)acetic acid", id="acid_outranks_the_nitrate"),
        pytest.param("O=[N+]([O-])OCCOC(C)=O", "2-(nitrooxy)ethyl acetate", id="carboxylic_ester_outranks_the_nitrate"),
    ],
)
def test_nitrate_and_nitrite_esters_outrank_the_groups_of_their_organyl_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1CCCCC1N=O", "nitrosocyclohexane", id="ring__nitroso"),
        pytest.param("C=CN=O", "nitrosoethene", id="unsaturated_chain__nitroso"),
    ],
)
def test_ring__nitroso_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_nitroso_chain():
    assert smiles_to_iupac("c1ccccc1CN=O") == "(nitrosomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CCN=O") == "(2-nitrosoethyl)benzene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=Cc1ccccc1CN=O", "1-ethenyl-2-(nitrosomethyl)benzene", id="phenyl_nitroso_unsaturation"),
        pytest.param("CN1CCCCCC1", "1-methylazepane", id="methylazepane"),
        pytest.param("CN1CCNCC1", "1-methylpiperazine", id="piperazine"),
        pytest.param("C=CN1CCCCC1", "1-ethenylpiperidine", id="unsaturated_n_substituent__ring_amine"),
        pytest.param("CN1CC(C)C(C)C1", "1,3,4-trimethylpyrrolidine", id="ring_amine_two_ring_carbon_substituents"),
        pytest.param("CS(=O)(=O)N1CCC(C)CC1", "1-(methanesulfonyl)-4-methylpiperidine", id="ring_amine_sulfonyl_with_ring_carbon_substituent"),
        pytest.param("CS(=O)(=O)N1CCOCC1", "4-methylsulfonylmorpholine", id="methylsulfonylmorpholine"),
        pytest.param("NC(=O)N", "urea", id="urea"),
        pytest.param("CNC(=O)N", "methylurea", id="n_methylurea"),
        pytest.param("CN(C)C(=O)N", "N,N-dimethylurea", id="n_n_dimethylurea_same_nitrogen"),
        pytest.param("CCN(C)C(=O)N", "N-ethyl-N-methylurea", id="n_ethyl_n_methylurea_same_nitrogen"),
    ],
)
def test_phenyl_nitroso_unsaturation_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CNC(=O)Nc1ccccc1", "N-methyl-N'-phenylurea", id="methyl_n_prime_phenylurea_different_nitrogens"),
        pytest.param("c1ccc(NC(=O)Nc2ccccc2)cc1", "N,N'-diphenylurea", id="n_prime_diphenylurea"),
        pytest.param("CCN(C)C(=O)NC", "N-ethyl-N,N'-dimethylurea", id="more_substituents_take_the_unprimed_locant"),
        pytest.param("O=C(Nc1ccccc1)Nc1ccccn1", "N-phenyl-N'-(pyridin-2-yl)urea", id="heteroaromatic_ring_n_substituent"),
        pytest.param("C=CNC(=O)Nc1ccc(Cl)cc1", "N-(4-chlorophenyl)-N'-ethenylurea", id="halogenated_ring_and_unsaturated_n_substituents"),
        pytest.param("CCNC(=O)N[Si](C)(C)C", "N-ethyl-N'-(trimethylsilyl)urea", id="silyl_group_on_a_urea_nitrogen"),
        pytest.param("O=C(N[SiH3])N[SiH3]", "N,N'-disilylurea", id="same_element_group_on_both_urea_nitrogens"),
        pytest.param("NC(=O)N[GeH3]", "germylurea", id="element_hydride_atom_is_a_substituent_of_the_urea"),
        pytest.param("NC(=N)N[Ge](C)(C)C", "N-(trimethylgermyl)guanidine", id="germyl_group_on_a_guanidine_nitrogen"),
        pytest.param("NC(=N)NP", "N-phosphanylguanidine", id="phosphanyl_group_outranked_by_the_guanidine"),
        pytest.param("CNC(=O)NP(C)C", "N-(dimethylphosphanyl)-N'-methylurea", id="phosphanyl_group_outranked_by_the_urea"),
        pytest.param("NC(=N)N(C)[SiH3]", "N-methyl-N-silylguanidine", id="silyl_and_carbon_group_on_one_guanidine_nitrogen"),
    ],
)
def test_n_methyl_n_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_semicarbazide():
    assert smiles_to_iupac("NC(=O)NN") == "hydrazinecarboxamide"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCCCC1=NN", "cyclohexylidenehydrazine"),
        ("C1CCCCC1=NNC", "cyclohexylidene(methyl)hydrazine"),
        ("CC(C)=NN(C)C", "1,1-dimethyl-2-(propan-2-ylidene)hydrazine"),
        ("CC=N[NH]Cl", "chloro(ethylidene)hydrazine"),
        ("c1ccccc1C=NN", "benzylidenehydrazine"),
        ("OC1CCCCC1C(=NN)C1CCCCC1O", "2,2'-(hydrazinylidenemethylene)di(cyclohexan-1-ol)"),
    ],
)
def test_hydrazone_with_ring_or_substituted_nitrogen_is_ylidene_hydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCCCC1=NO", "N-hydroxycyclohexanimine"),
        ("ClC1CCCCC1=NC", "2-chloro-N-methylcyclohexan-1-imine"),
        ("OC1CCCCC1=N", "2-iminocyclohexan-1-ol"),
        pytest.param("CN=CC1CCCCC1=N", "2-(N-methylmethanimidoyl)cyclohexan-1-imine", id="n_substituent_stays_on_the_prefix"),
        pytest.param(
            "CN=CC1CCCCC1=NCC", "N-ethyl-2-(N-methylmethanimidoyl)cyclohexan-1-imine", id="parent_and_prefix_both_substituted"
        ),
    ],
)
def test_ring_imine_and_oxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-62.3.1.1, P-62.2.1.3, P-41: a ring nitrogen on the imine nitrogen is the group 'piperidin-1-yl'; an amino group on a ring nitrogen is the amine suffix
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC1CCC(=NN2CCCCC2)CC1", "4-methyl-N-(piperidin-1-yl)cyclohexan-1-imine", id="alkyl_on_the_ring"),
        pytest.param("CC1CCCCC1=NN1CCCCC1", "2-methyl-N-(piperidin-1-yl)cyclohexan-1-imine", id="alkyl_beside_the_imine_carbon"),
        pytest.param("CC1CCC(=NN2CCCC2)C1", "3-methyl-N-(pyrrolidin-1-yl)cyclopentan-1-imine", id="five_membered_rings"),
        pytest.param("ClC1CCC(=NN2CCCCC2)CC1", "4-chloro-N-(piperidin-1-yl)cyclohexan-1-imine", id="halogen_on_the_ring"),
        pytest.param("CC1CCC(=NN2CCOCC2)CC1", "4-methyl-N-(morpholin-4-yl)cyclohexan-1-imine", id="second_heteroatom_in_the_nitrogen_ring"),
        pytest.param(
            "CC1CCC(=NN2CCN(C)CC2)CC1", "4-methyl-N-(4-methylpiperazin-1-yl)cyclohexan-1-imine", id="substituted_nitrogen_ring"
        ),
        pytest.param("N=C1CCN(N)CC1", "4-iminopiperidin-1-amine", id="amine_on_a_ring_nitrogen_outranks_the_imine"),
        pytest.param("CN=C1CCN(N)CC1", "4-(methylimino)piperidin-1-amine", id="substituted_imino_prefix"),
        pytest.param("ON=C1CCN(N)CC1", "4-(hydroxyimino)piperidin-1-amine", id="oxime_prefix"),
        pytest.param("N=C1CCCN(N)C1", "3-iminopiperidin-1-amine", id="imine_at_position_3"),
        pytest.param("N=C1CCN(N)C1", "3-iminopyrrolidin-1-amine", id="five_membered_ring_amine"),
        pytest.param("O=C1CCN(N)CC1", "1-aminopiperidin-4-one", id="ketone_outranks_the_amine"),
    ],
)
def test_ring_imine_with_a_ring_nitrogen_group_or_an_amino_ring_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)=Nc1ccccc1", "N-phenylpropan-2-imine", id="n_aryl_imine"),
        pytest.param("c1ccccc1C=Nc1ccccc1", "N,1-diphenylmethanimine", id="n_aryl_aldimine_of_an_aryl_aldehyde"),
        pytest.param("CC=NOCc1ccccc1", "N-(benzyloxy)ethanimine", id="oxime_ether_with_an_aromatic_chain"),
        pytest.param("CC=NOC(C)C", "N-[(propan-2-yl)oxy]ethanimine", id="branched_oxime_o_substituent"),
        pytest.param("CSN=CC", "N-(methylsulfanyl)ethanimine", id="sulfanyl_on_the_imine_nitrogen"),
        pytest.param("CSN=C1CCCCC1", "N-(methylsulfanyl)cyclohexanimine", id="sulfanyl_on_the_nitrogen_of_a_ring_imine"),
        pytest.param("C[Se]N=CC", "N-(methylselanyl)ethanimine", id="selanyl_on_the_imine_nitrogen"),
        pytest.param("CSSN=CC", "N-(methyldisulfanyl)ethanimine", id="disulfanyl_on_the_imine_nitrogen"),
        pytest.param("N=CC=N", "ethane-1,2-diimine", id="polyimine"),
        pytest.param("CC(=NO)N", "N'-hydroxyethanimidamide", id="amidoxime"),
        pytest.param("NC(=NO)c1ccccc1", "N'-hydroxybenzenecarboximidamide", id="aryl_amidoxime"),
        pytest.param("CCC(=NO)NOCC", "N-ethoxy-N'-hydroxypropanimidamide", id="oxygen_groups_on_both_amidine_nitrogens"),
        pytest.param("CC(N)=NOCC(=O)OC", "methyl {[(1-aminoethylidene)amino]oxy}acetate", id="amidoxime_ether_beside_an_ester"),
        pytest.param("CC(N)=NOCC(=O)O", "{[(1-aminoethylidene)amino]oxy}acetic acid", id="amidoxime_ether_beside_an_acid"),
        pytest.param("CC(NC)=NOCC(=O)OC", "methyl ({[1-(methylamino)ethylidene]amino}oxy)acetate", id="substituted_amidoxime_ether_beside_an_ester"),
        pytest.param("NC(=NOCC(=O)OCC)c1ccccc1", "ethyl ({[amino(phenyl)methylidene]amino}oxy)acetate", id="aryl_amidoxime_ether_beside_an_ester"),
        pytest.param("NC(=N)NOCC", "N-ethoxyguanidine", id="alkoxy_group_on_a_guanidine_nitrogen"),
        pytest.param("N=C(N)c1ccccc1C(=N)N", "benzene-1,2-dicarboximidamide", id="ring_diamidine"),
        pytest.param("CC(=NC)N(C)C", "N,N,N'-trimethylethanimidamide", id="n_and_n_prime_substituents"),
        pytest.param("CC(=N[Si](C)(C)C)N[Si](C)(C)C", "N,N'-bis(trimethylsilyl)ethanimidamide", id="hydride_groups_on_both_amidine_nitrogens"),
        pytest.param("CC(=N[Si](C)(C)C)NC", "N-methyl-N'-(trimethylsilyl)ethanimidamide", id="hydride_group_beside_a_carbon_group"),
        pytest.param("NC(=N)CC(=N)N", "propanediimidamide", id="chain_diamidine"),
        pytest.param("CC(=N)CN=CC", "1-(ethylideneamino)propan-2-imine", id="second_imine_as_ylideneamino_prefix"),
        pytest.param(
            "CC(=N)CN=Cc1c(Cl)cccc1Cl",
            "1-{[(2,6-dichlorophenyl)methylidene]amino}propan-2-imine",
            id="aryl_methylideneamino_prefix",
        ),
        pytest.param("N=C1CCCC1CN=CC", "2-[(ethylideneamino)methyl]cyclopentan-1-imine", id="ring_imine_parent_with_ylideneamino_prefix"),
        pytest.param("CC(=N)CN=CCCC", "N-(2-iminopropyl)butan-1-imine", id="longer_imine_chain_is_the_parent"),
        pytest.param("CCC=NCCN=CCC", "N,N'-(ethane-1,2-diyl)di(propan-1-imine)", id="identical_imines_joined_multiplicatively"),
        pytest.param("c1ccccc1C=NCN=Cc1ccccc1", "N,N'-methylenebis(phenylmethanimine)", id="substituted_imine_units_use_bis"),
        pytest.param("C=NC1CCC(N=C)CC1", "N,N'-(cyclohexane-1,4-diyl)dimethanimine", id="imine_units_on_a_ring_linker"),
        pytest.param("CC=NC(C)N=CC", "N,N'-(ethane-1,1-diyl)diethanimine", id="imine_units_on_a_geminal_linker"),
        pytest.param("C=NCC(C)N=C", "N,N'-(propane-1,2-diyl)dimethanimine", id="imine_units_on_an_unsymmetrical_linker"),
        pytest.param("CC(C)=NCC(C)N=C(C)C", "N,N'-(propane-1,2-diyl)di(propan-2-imine)", id="imine_unit_with_a_locant_uses_di_and_enclosing_marks"),
        pytest.param("c1ccccc1C=Nc1ccc(N=Cc2ccccc2)cc1", "N,N'-(1,4-phenylene)bis(phenylmethanimine)", id="substituted_imine_units_on_an_arene_linker"),
        pytest.param("CC=NN1CCCC1", "N-(pyrrolidin-1-yl)ethanimine", id="imine_nitrogen_on_a_ring_nitrogen"),
        pytest.param("O=CC=NN1CCCC1", "[(pyrrolidin-1-yl)imino]acetaldehyde", id="ring_nitrogen_imino_prefix_under_an_aldehyde"),
    ],
)
def test_imine_and_amidine_n_substituents_and_polyfunctional_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC=NNC(N)=O", "2-ethylidenehydrazine-1-carboxamide", id="semicarbazone"),
        # P-93.4.2.1.3, P-93.5.1.4.2.1: the C=N descriptor of a ring ylidene is cited with the locant of the attachment, or alone when the hydrazine has no locant
        pytest.param("N/N=C1/CCCCC1C", "(Z)-(2-methylcyclohexylidene)hydrazine", id="ring_ylidene_z"),
        pytest.param("N/N=C1\\CCCCC1C", "(E)-(2-methylcyclohexylidene)hydrazine", id="ring_ylidene_e"),
        pytest.param("N/N=C1/CCCC1C", "(Z)-(2-methylcyclopentylidene)hydrazine", id="five_membered_ring_ylidene"),
        pytest.param("CN/N=C1/CCCCC1C", "methyl[(1Z)-2-methylcyclohexylidene]hydrazine", id="ring_ylidene_beside_a_substituent"),
        pytest.param("CN(C)/N=C1/CCCCC1C", "1,1-dimethyl-2-[(1Z)-2-methylcyclohexylidene]hydrazine", id="ring_ylidene_with_locants"),
        pytest.param("NC(=O)N/N=C1/CCCCC1C", "2-[(1Z)-2-methylcyclohexylidene]hydrazine-1-carboxamide", id="ring_semicarbazone_z"),
        pytest.param("NC(=O)N/N=C1\\CCCCC1C", "2-[(1E)-2-methylcyclohexylidene]hydrazine-1-carboxamide", id="ring_semicarbazone_e"),
        pytest.param("NC(=S)N/N=C1/CCCCC1C", "2-[(1Z)-2-methylcyclohexylidene]hydrazine-1-carbothioamide", id="ring_thiosemicarbazone"),
        pytest.param("CC(=O)N/N=C1/CCCCC1C", "N'-[(1Z)-2-methylcyclohexylidene]acetohydrazide", id="ring_acylhydrazone"),
        pytest.param("c1ccccc1N/N=C1/CCCCC1C", "[(1Z)-2-methylcyclohexylidene](phenyl)hydrazine", id="ring_phenylhydrazone"),
        pytest.param("CNC(=O)NNC", "N,2-dimethylhydrazine-1-carboxamide", id="substituted_semicarbazide"),
        pytest.param("CC=NNC(=O)c1ccccc1", "N'-ethylidenebenzohydrazide", id="acylhydrazone"),
        pytest.param("O=C(NN=Cc1ccccc1)c1ccncc1", "N'-benzylidenepyridine-4-carbohydrazide", id="ring_carbohydrazide_hydrazone"),
        pytest.param("CC(=O)N(c1ccccc1)Nc1ccccc1", "N,N'-diphenylacetohydrazide", id="n_and_n_prime_substituents"),
        pytest.param("NC(=Nc1ccccn1)N", "N''-(pyridin-2-yl)guanidine", id="guanidine_with_heteroaromatic_ring"),
        pytest.param("N=C(N)NC1=NCCN1", "N-(4,5-dihydro-1H-imidazol-2-yl)guanidine", id="guanidine_on_a_ring_amidine_carbon_is_not_a_condensed_guanidine"),
        pytest.param("CCCCCN=C(N)NC1=NCCCN1", "N''-pentyl-N-(1,4,5,6-tetrahydropyrimidin-2-yl)guanidine", id="ring_amidine_group_with_a_substituent_on_the_imino_nitrogen"),
        pytest.param("CN(C)C(=NC)N(C)C", "N,N,N',N',N''-pentamethylguanidine", id="guanidine_with_substituents_on_every_nitrogen"),
        # P-66.4.1.1, P-68.3.1.2.4: an amidine on a hydrazine nitrogen is a suffix of the hydrazine; N and N' are the amidine nitrogens
        pytest.param("CC=NN(C)C(N)=N", "2-ethylidene-1-methylhydrazine-1-carboximidamide", id="aminoguanidine_hydrazone_with_a_substituent"),
        pytest.param("CNC(=N)NN", "N-methylhydrazinecarboximidamide", id="substituent_on_the_amidine_amino_nitrogen"),
        pytest.param("CN=C(N)NN", "N'-methylhydrazinecarboximidamide", id="substituent_on_the_amidine_imino_nitrogen"),
        # P-41: an amidine outranks a hydrazide, so an acyl or sulfonyl group on the hydrazine nitrogens is a prefix
        pytest.param("NC(=N)NNC(C)=O", "2-acetylhydrazine-1-carboximidamide", id="acyl_on_the_terminal_hydrazine_nitrogen"),
        pytest.param("CNC(=N)NNC(C)=O", "2-acetyl-N-methylhydrazine-1-carboximidamide", id="acyl_with_a_substituent_on_the_amidine"),
        pytest.param("NC(=N)N(C)NC(C)=O", "2-acetyl-1-methylhydrazine-1-carboximidamide", id="acyl_beside_a_substituent_on_the_amidine_hydrazine_nitrogen"),
        pytest.param("NC(=N)NNC(=O)c1ccccc1", "2-benzoylhydrazine-1-carboximidamide", id="aroyl_on_the_terminal_hydrazine_nitrogen"),
        pytest.param("NC(=N)NNS(C)(=O)=O", "2-(methanesulfonyl)hydrazine-1-carboximidamide", id="sulfonyl_on_the_terminal_hydrazine_nitrogen"),
        pytest.param("NC(=O)NNC(C)=O", "N'-carbamoylacetohydrazide", id="urea_keeps_the_acyl_hydrazide_as_the_parent"),
        # P-66.4.1.4.1, P-66.4.1.4.2: several amidines on a hydrazine; N and N' carry the position of the group, geminal groups continue N'', N'''
        pytest.param("NC(=N)NNC(=N)N", "hydrazine-1,2-dicarboximidamide", id="amidine_on_each_hydrazine_nitrogen"),
        pytest.param("CCNC(=N)NN(C)C(N)=N", "N2-ethyl-1-methylhydrazine-1,2-dicarboximidamide", id="hydrazine_substituent_fixes_the_numbering"),
        pytest.param("CN=C(N)NNC(N)=N", "N'1-methylhydrazine-1,2-dicarboximidamide", id="imino_nitrogen_of_a_group_by_position"),
        pytest.param("CC=NN(C(N)=N)C(N)=N", "2-ethylidenehydrazine-1,1-dicarboximidamide", id="geminal_groups_with_a_hydrazone"),
        pytest.param("CCN=C(N)N(N)C(=N)N(C)C", "N'''1-ethyl-N1,N1-dimethylhydrazine-1,1-dicarboximidamide", id="geminal_more_substituted_group_unprimed"),
        pytest.param("CNC(=N)NN(C(N)=N)C(N)=N", "N2-methylhydrazine-1,1,2-tricarboximidamide", id="geminal_pair_and_a_single_group"),
    ],
)
def test_hydrazide_semicarbazone_and_guanidine_n_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param(
            "O=C(NN=Cc1ccccc1)N(C)N=CC",
            "N'-benzylidene-2-ethylidene-1-methylhydrazine-1-carbohydrazide",
            id="twice_substituted_unit_is_the_parent_hydrazine",
        ),
        pytest.param("NNC(=O)N(C)N", "1-methylhydrazine-1-carbohydrazide", id="single_substituent_on_the_parent_unit"),
        pytest.param("CNNC(=O)NNC", "N',2-dimethylhydrazine-1-carbohydrazide", id="equal_units_cite_n_prime_locant"),
        pytest.param("NNC(=S)NNC", "2-methylhydrazine-1-carbothiohydrazide", id="chalcogen_analogue"),
    ],
)
def test_substituted_hydrazinecarbohydrazides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)NC(=O)OC", "methyl acetylcarbamate", id="n_acyl_carbamate"),
        pytest.param("O=C(NC(=O)c1ccccc1)OC", "methyl benzoylcarbamate", id="n_aroyl_carbamate"),
        pytest.param("CC(=O)N(C)C(=O)OC", "methyl acetyl(methyl)carbamate", id="n_acyl_n_alkyl_carbamate"),
    ],
)
def test_n_acyl_carbamates(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C(NCCC(N)=N)C1CCC(N)C1", "3-amino-N-(3-amino-3-iminopropyl)cyclopentane-1-carboxamide", id="amidine_ends_a_substituent_chain"),
        pytest.param("OC(=O)c1ccccc1CCC(NC)=NCC", "2-[3-(ethylimino)-3-(methylamino)propyl]benzoic acid", id="substituted_amidine_ends_a_substituent_chain"),
        pytest.param("OC(=O)c1ccc(cc1)C(N)=N", "4-carbamimidoylbenzoic acid", id="amidine_on_a_ring_stays_carbamimidoyl"),
    ],
)
def test_amidine_carbon_in_a_substituent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acylated_nitrogen_is_not_a_hetero_parent_for_a_carboxylic_acid():
    assert smiles_to_iupac("CC(=O)NC(=O)O") == "acetylcarbamic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CNC(=O)CCCC(=O)NC", "N1,N5-dimethylpentanediamide", id="n_locants_of_a_chain_diamide"),
        pytest.param("CNC(=O)CCC(N)=O", "N1-methylbutanediamide", id="one_substituted_group_of_a_diamide"),
        pytest.param("CCNC(=O)c1cccc(C(=O)NC)c1", "N1-ethyl-N3-methylbenzene-1,3-dicarboxamide", id="n_locants_of_a_ring_diamide"),
        pytest.param("CC(=O)NCC(N)=O", "2-acetamidoacetamide", id="second_amide_cited_as_a_prefix"),
        pytest.param("CNC(=O)CCc1cc(CC(C)C(N)=O)ccc1C", "3-[5-(3-amino-2-methyl-3-oxopropyl)-2-methylphenyl]-N-methylpropanamide", id="parent_amide_carries_the_n_prefix"),
        pytest.param("CN(C)C(=O)C(C)Cc1ccc(C)c(CC(C)C(=O)NC)c1", "N,N,2-trimethyl-3-{4-methyl-3-[2-methyl-3-(methylamino)-3-oxopropyl]phenyl}propanamide", id="parent_amide_with_more_prefixes"),
        pytest.param("CNC(CCC1=CC(=C(C=C1)C)CCC(NCCC)=O)=O", "N-methyl-3-{4-methyl-3-[3-oxo-3-(propylamino)propyl]phenyl}propanamide", id="equal_locant_sets_compared_in_citation_order_with_n_first"),
        pytest.param("C(C)NC(=O)C(CCC(=O)NC)C(=O)N(C)C", "N'1-ethyl-N1,N1,N3-trimethylpropane-1,1,3-tricarboxamide", id="primed_n_locants_on_a_geminal_carboxamide"),
        pytest.param("CNC(=O)C(CC(=O)NC)CC(=O)NC", "N1,N2,N3-trimethylpropane-1,2,3-tricarboxamide", id="n_substituents_of_carboxamides_on_a_chain"),
        pytest.param("O=CCNSCNON=C(C)S", "N-{[({[(2-oxoethyl)amino]sulfanyl}methyl)amino]oxy}ethanimidothioic acid", id="n_substituent_of_an_imidothioic_acid"),
    ],
)
def test_n_locants_of_polyamides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(OC(=O)c1ccccc1)C(=O)Nc1ccccc1", "1-anilino-1-oxopropan-2-yl benzoate", id="n_phenyl_carbamoyl_prefix"),
        pytest.param("COC(=O)NNC(=S)Nc1ccccc1", "methyl 2-(phenylcarbamothioyl)hydrazine-1-carboxylate", id="n_phenyl_thiocarbamoyl_prefix"),
        pytest.param("O=C(CS(=O)(=O)Nc1ccccc1)C(=O)OC", "methyl 2-oxo-3-(phenylsulfamoyl)propanoate", id="n_phenyl_sulfamoyl_prefix"),
        pytest.param("CC(=O)OCON(O)O", "[(dihydroxyamino)oxy]methyl acetate", id="n_dihydroxy_is_not_nitro"),
        pytest.param("CC(=O)OCO[N+](=O)[O-]", "(nitrooxy)methyl acetate", id="nitrate_ester_prefix"),
    ],
)
def test_nitrogen_prefix_keeps_every_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=NN)N", "ethanehydrazonamide", id="hydrazonamide"),
        pytest.param("NC(=NN)c1ccccc1", "benzenecarbohydrazonamide", id="carbohydrazonamide_on_a_ring"),
        pytest.param("C(=N)NN", "methanimidohydrazide", id="imidohydrazide_of_a_one_carbon_parent"),
        pytest.param("NNC=NN", "methanehydrazonohydrazide", id="hydrazonohydrazide_of_a_one_carbon_parent"),
        pytest.param("CNC=NN", "N-methylmethanehydrazonamide", id="hydrazonamide_of_a_one_carbon_parent"),
        pytest.param("N=C(NN)C1CCCCC1", "cyclohexanecarboximidohydrazide", id="carboximidohydrazide_on_a_ring"),
        pytest.param("CN(C)C=NN=C(C)C", "N,N-dimethyl-N'-(propan-2-ylidene)methanehydrazonamide", id="ylidene_on_the_terminal_nitrogen"),
        pytest.param("CCN=C(c1ccccc1)N(C)Nc1ccccc1", "N''-ethyl-N-methyl-N'-phenylbenzenecarboximidohydrazide", id="three_nitrogen_locants"),
        pytest.param("NC(=NN)C(N)=NN", "ethanedihydrazonamide", id="two_hydrazonamide_groups"),
        pytest.param("N=C(NN)N1CCCC1", "pyrrolidine-1-carboximidohydrazide", id="imidohydrazide_on_a_ring_nitrogen"),
        pytest.param("N=C(NN)N1CCOCC1", "morpholine-4-carboximidohydrazide", id="imidohydrazide_on_a_second_heteroatom_ring"),
        pytest.param("N=C(NN)N1CCCC1C", "2-methylpyrrolidine-1-carboximidohydrazide", id="substituent_on_the_nitrogen_ring"),
        pytest.param("CN=C(NN)N1CCCC1", "N''-methylpyrrolidine-1-carboximidohydrazide", id="imino_nitrogen_of_a_ring_imidohydrazide"),
        pytest.param("N=C(NNC)N1CCCC1", "N'-methylpyrrolidine-1-carboximidohydrazide", id="terminal_nitrogen_of_a_ring_imidohydrazide"),
        pytest.param("NN=C(N)N1CCCC1", "pyrrolidine-1-carbohydrazonamide", id="hydrazonamide_on_a_ring_nitrogen"),
        pytest.param("NN=C(N)N1CCCCC1", "piperidine-1-carbohydrazonamide", id="hydrazonamide_on_a_six_membered_ring"),
        pytest.param("NNC(=N)C(=N)NN", "ethanediimidohydrazide", id="two_imidohydrazide_groups"),
        pytest.param("NC=N", "methanimidamide", id="amidine_of_a_one_carbon_parent"),
        pytest.param("CN=CNC1CCOC1", "N'-methyl-N-(oxolan-3-yl)methanimidamide", id="one_carbon_amidine_with_a_ring_on_nitrogen"),
        pytest.param(
            "CN(C)C(=NN)c1cc(C(=NN)N(CC)CC)cc2ccccc12",
            "N3,N3-diethyl-N1,N1-dimethylnaphthalene-1,3-dicarbohydrazonamide",
            id="substituted_hydrazonamides_on_two_ring_positions",
        ),
        pytest.param(
            "CCN(CC)NC(=N)C1(CCCCC1)C(=N)NN(C)C",
            "N'1,N'1-diethyl-N''' '1,N''' '1-dimethylcyclohexane-1,1-dicarboximidohydrazide",
            id="geminal_imidohydrazides_continue_the_primes",
        ),
        pytest.param("NN=C(N)C(=N)NN", "2-hydrazinyl-2-iminoethanehydrazonamide", id="imidohydrazide_at_a_chain_end_is_prefixes"),
        pytest.param("NC(=N)CC(=N)NN", "3-hydrazinyl-3-iminopropanimidamide", id="imidohydrazide_beside_an_amidine"),
        pytest.param("NN=S(N)c1ccccc1", "benzenesulfinohydrazonamide", id="sulfinic_amidrazone"),
        pytest.param("CS(=S)(=NN)N", "methanesulfonohydrazonothioamide", id="replaced_oxygen_in_a_sulfonohydrazonamide"),
        pytest.param("CS(=O)(=NNC)N", "N'-methylmethanesulfonohydrazonamide", id="sulfur_centred_hydrazone_terminal_is_primed"),
        pytest.param("CS(=O)(=NNC)NNC", "N',N''-dimethylmethanesulfonohydrazonohydrazide", id="sulfur_centred_hydrazonohydrazide_locants"),
        pytest.param("CS(=NN)(=NN)N", "methanesulfonodihydrazonamide", id="two_hydrazono_positions_on_sulfur"),
        pytest.param("OC(=O)CCS(=O)(=NN)N", "3-(S-aminosulfonohydrazonoyl)propanoic acid", id="sulfur_centred_hydrazonamide_as_a_prefix"),
    ],
)
def test_amidrazones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CNO", "N-hydroxymethanamine", id="n_hydroxy_amine"),
        pytest.param("CNCl", "methylhypochlorous amide", id="amide_of_hypochlorous_acid"),
        pytest.param("ClN(Cl)C1CCC1", "chloro(cyclobutyl)hypochlorous amide", id="second_halogen_on_nitrogen_is_a_prefix"),
        pytest.param("FNN1CCCC1", "(pyrrolidin-1-yl)hypofluorous amide", id="halogen_amide_of_a_ring_nitrogen_group"),
        pytest.param("ClNC1CCCCO1", "(oxan-2-yl)hypochlorous amide", id="halogen_amide_of_a_ring_ether_group"),
        pytest.param("CCON(C)Cl", "ethoxy(methyl)hypochlorous amide", id="alkoxy_group_on_the_nitrogen_of_a_halogen_amide"),
        pytest.param("CN(Br)SC", "methyl(methylsulfanyl)hypobromous amide", id="sulfanyl_group_on_the_nitrogen_of_a_halogen_amide"),
        pytest.param("O=C(N(F)F)N(F)F", "tetrafluorourea", id="urea_all_positions_halogenated"),
        pytest.param("CC(=O)N(Cl)C", "N-chloro-N-methylacetamide", id="n_halogen_amide"),
        pytest.param("CCCN(N=O)C(=N)N[N+](=O)[O-]", "N'-nitro-N-nitroso-N-propylguanidine", id="guanidine_nitro_nitroso"),
    ],
)
def test_heteroatom_substituents_on_nitrogen_and_pseudohalide_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(N)=S", "ethanethioamide", id="chain_thioamide"),
        pytest.param("NC=S", "methanethioamide", id="one_carbon_thioamide_has_no_retained_name"),
        pytest.param("NC(=S)CC(N)=S", "propanedithioamide", id="dithioamide_on_a_chain"),
        pytest.param("NC(=S)c1ccccn1", "pyridine-2-carbothioamide", id="carbothioamide_on_a_heteroarene"),
        pytest.param("NC(=S)c1ccc(Cl)cc1", "4-chlorobenzene-1-carbothioamide", id="carbothioamide_with_substituent_locants"),
        pytest.param("NC(=S)c1ccc(cc1)C(N)=S", "benzene-1,4-dicarbothioamide", id="carbothioamides_on_one_ring_are_not_multiplicative"),
        pytest.param("NC(=S)CC(CC(N)=S)C(N)=S", "propane-1,2,3-tricarbothioamide", id="carbothioamides_on_a_chain"),
        pytest.param("CCC(=S)N(C)C", "N,N-dimethylpropanethioamide", id="n_substituted_thioamide"),
        pytest.param("ON(C)C(C)=S", "N-hydroxy-N-methylethanethioamide", id="n_hydroxy_thioamide"),
        pytest.param("CC(N)=[Se]", "ethaneselenoamide", id="selenoamide"),
        pytest.param("CC(N)=[Te]", "ethanetelluroamide", id="telluroamide"),
        pytest.param("CC(=S)N1CCCC1", "1-(pyrrolidin-1-yl)ethane-1-thione", id="hidden_amide_is_a_thione"),
        pytest.param("CCC(=S)NC(C)=O", "N-(propanethioyl)acetamide", id="n_acyl_amide_keeps_the_oxygen_amide_as_parent"),
        pytest.param("CC(=S)NC(C)=S", "N-(ethanethioyl)ethanethioamide", id="identical_thioacyl_groups_on_nitrogen"),
        pytest.param("CC(=S)N(C1CCCCC1)C(C)=S", "N-cyclohexyl-N-(ethanethioyl)ethanethioamide", id="tertiary_thioamide_with_a_third_substituent"),
        pytest.param("CC(=S)NC(=S)c1ccccc1", "N-(ethanethioyl)benzenecarbothioamide", id="imide_parent_is_the_ring_acyl_group"),
        pytest.param("O=C(O)c1ccc(C=[Se])cc1", "4-(methaneselenoyl)benzoic acid", id="methane_chalcogenoyl_prefix_is_compound"),
        pytest.param("NC(=O)CC(N)=S", "3-amino-3-sulfanylidenepropanamide", id="thioamide_under_an_amide_joins_the_chain"),
        pytest.param("NC(=[Se])CC(N)=S", "3-amino-3-selanylidenepropanethioamide", id="selenoamide_under_a_thioamide"),
        pytest.param("OC(=O)c1ccc(cc1)C(N)=S", "4-carbamothioylbenzoic acid", id="thioamide_prefix_under_an_acid"),
        pytest.param("NC(=[Se])c1ccc(cc1)C(=O)O", "4-carbamoselenoylbenzoic acid", id="selenoamide_prefix_under_an_acid"),
        pytest.param("NC(=O)c1ccc(NC(C)=S)cc1", "4-(ethanethioamido)benzamide", id="thioacylamino_prefix_is_an_amido_prefix"),
        pytest.param("OC(=O)CCN(C)C(=S)C(Cl)C", "3-(2-chloro-N-methylpropanethioamido)propanoic acid", id="n_substituted_amido_prefix_merges_with_the_acyl_prefixes"),
        pytest.param("NC(=S)NCCC(=O)O", "3-(carbamothioylamino)propanoic acid", id="thiourea_group_stays_an_acylamino_prefix"),
        pytest.param("OC(=O)c1ccccc1C(=S)C(N)=S", "2-[amino(sulfanylidene)ethanethioyl]benzoic acid", id="completely_substituted_acyl_group_omits_locants"),
    ],
)
def test_chalcogen_analogues_of_amides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CS(=S)(=O)N", "methanesulfonothioamide", id="sulfonothioamide"),
        pytest.param("NS(=S)(=S)c1ccc2ccccc2c1", "naphthalene-2-sulfonodithioamide", id="sulfonodithioamide_on_a_ring"),
        pytest.param("CS(=S)N", "methanesulfinothioamide", id="sulfinothioamide"),
        pytest.param("C[Te](=O)(=S)N", "methanetelluronothioamide", id="replaced_oxygen_on_a_tellurium_centre"),
        pytest.param("CS(=S)(=[Se])N", "methanesulfonoselenothioamide", id="infixes_in_alphanumerical_order"),
        pytest.param("CS(=S)(=O)NC", "N-methylmethanesulfonothioamide", id="n_substituted_sulfonothioamide"),
        pytest.param("OC(=O)CNS(C)=S", "[(methanesulfinothioyl)amino]acetic acid", id="sulfinothioyl_amino_prefix"),
    ],
)
def test_chalcogen_analogues_of_sulfonamides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NC(=O)NCCO", "(2-hydroxyethyl)urea", id="hydroxy_group_in_an_n_substituent"),
        pytest.param("COCCNC(=O)Nc1cccs1", "N-(2-methoxyethyl)-N'-(thiophen-2-yl)urea", id="ether_and_heteroaryl_on_the_two_nitrogens"),
        pytest.param("NC(=O)N(O)C", "N-hydroxy-N-methylurea", id="hydroxy_on_the_nitrogen_itself"),
        pytest.param("N#CC(CCSC)NC(=O)NC", "N-[1-cyano-3-(methylsulfanyl)propyl]-N'-methylurea", id="nitrile_and_sulfanyl_in_a_substituent"),
        pytest.param("NC(=S)NCCO", "(2-hydroxyethyl)thiourea", id="thiourea_with_a_hydroxy_substituent"),
        pytest.param("CC(C)=NC(=S)N(C)C", "N,N-dimethyl-N'-(propan-2-ylidene)thiourea", id="ylidene_on_a_thiourea_nitrogen"),
        pytest.param("CNC(=S)NCCS(C)=O", "N-[2-(methanesulfinyl)ethyl]-N'-methylthiourea", id="sulfinyl_in_a_thiourea_substituent"),
        pytest.param("CNC(=S)NCCS(C)(=O)=O", "N-[2-(methanesulfonyl)ethyl]-N'-methylthiourea", id="sulfonyl_in_a_thiourea_substituent"),
        pytest.param("NC(=S)NCCS(C)=O", "[2-(methanesulfinyl)ethyl]thiourea", id="sulfinyl_on_an_otherwise_free_thiourea"),
        pytest.param("CNC(=S)N(C)CCS(C)=O", "N-[2-(methanesulfinyl)ethyl]-N,N'-dimethylthiourea", id="sulfinyl_beside_two_methyl_groups"),
        pytest.param("CNC(=S)NCCS(=O)c1ccccc1", "N-[2-(benzenesulfinyl)ethyl]-N'-methylthiourea", id="aryl_sulfinyl_in_a_thiourea_substituent"),
        pytest.param("CNC(=[Se])NCCS(C)=O", "N-[2-(methanesulfinyl)ethyl]-N'-methylselenourea", id="sulfinyl_in_a_selenourea_substituent"),
        pytest.param("NC(=O)NC(=O)c1ccccc1", "N-carbamoylbenzamide", id="urea_beneath_a_carboxamide"),
        pytest.param("NC(=O)NS(=O)(=O)c1ccccc1", "N-carbamoylbenzenesulfonamide", id="urea_beneath_a_sulfonamide"),
        pytest.param("CNC(=S)NC=S", "N-(methylcarbamothioyl)methanethioamide", id="thiourea_beneath_a_thioamide"),
        pytest.param("NC(=S)NS(C)(=O)=O", "N-carbamothioylmethanesulfonamide", id="thiourea_beneath_a_sulfonamide"),
        pytest.param("CC(=O)NCCNC(N)=O", "N-[2-(carbamoylamino)ethyl]acetamide", id="carbamoylamino_prefix_beneath_an_amide"),
        pytest.param("O=CNCCCNC(N)=O", "N-[3-(carbamoylamino)propyl]formamide", id="formamide_outranks_urea"),
        pytest.param("NC(=O)NC(=O)O", "carbamoylcarbamic acid", id="carbamoyl_on_the_nitrogen_of_carbamic_acid"),
        pytest.param("NC(=O)NC(=O)OC", "methyl carbamoylcarbamate", id="ester_of_carbamoylcarbamic_acid"),
        pytest.param("NC(=O)NC(=O)Cl", "carbamoylcarbamoyl chloride", id="acyl_halide_of_carbamoylcarbamic_acid"),
        pytest.param("CNC(=O)NC(=O)O", "(methylcarbamoyl)carbamic acid", id="substituted_carbamoyl_on_carbamic_acid"),
    ],
)
def test_ureas_with_further_substituents_and_ureas_beneath_senior_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC(=O)NCCNC(=O)NCC", "N-[2-(carbamoylamino)ethyl]-N'-ethylurea", id="urea_with_more_substituents_is_the_parent"),
        pytest.param("c1ccccc1NC(=O)NCCNC(=O)NC", "N-{2-[(methylcarbamoyl)amino]ethyl}-N'-phenylurea", id="equal_substituent_counts_by_alphanumerical_order"),
        pytest.param("NC(=O)NCCNC(=S)NCC", "{2-[(ethylcarbamothioyl)amino]ethyl}urea", id="urea_outranks_thiourea"),
        pytest.param("CCNC(=O)NCC(C)NC(=O)NCC", "N,N''-(propane-1,2-diyl)bis(N'-ethylurea)", id="identical_units_stay_multiplicative_on_an_unsymmetrical_linker"),
        pytest.param("OC(=O)NCCN(C)C(=O)O", "[2-(carboxyamino)ethyl](methyl)carbamic acid", id="carbamic_acid_beside_carbamic_acid"),
        pytest.param("OC(=O)NCCC(=O)O", "3-(carboxyamino)propanoic acid", id="carboxylic_acid_outranks_carbamic_acid"),
    ],
)
def test_unequal_units_of_one_class_choose_the_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CN(C)N=O", "dimethylnitrous amide", id="nitroso_on_amine_nitrogen"),
        pytest.param("CN(N=O)N1CCCC1", "methyl(pyrrolidin-1-yl)nitrous amide", id="nitrous_amide_of_a_ring_nitrogen_group"),
        pytest.param("CN(C)[N+](=O)[O-]", "dimethylnitramide", id="nitro_on_amine_nitrogen"),
        pytest.param("NN[N+](=O)[O-]", "nitric hydrazide", id="nitro_on_hydrazine_nitrogen"),
        pytest.param("NNN=O", "nitrous hydrazide", id="nitroso_on_hydrazine_nitrogen"),
        pytest.param("CN([N+](=O)[O-])[N+](=O)[O-]", "methyl(nitro)nitramide", id="second_nitro_group_is_a_prefix_of_nitramide"),
        pytest.param("CCCCN(CC)N=O", "butyl(ethyl)nitrous amide", id="nitrous_amide_prefixes_without_locants"),
        pytest.param("O=C1CCCCC1N=O", "2-nitrosocyclohexan-1-one", id="nitroso_on_a_ring_that_carries_a_ketone"),
        pytest.param("CN(N=O)c1ccc(cc1)C(=O)O", "4-[methyl(nitroso)amino]benzoic acid", id="nitrosoamino_prefix_under_an_acid"),
        pytest.param("OC(=O)c1ccc(cc1)N[N+](=O)[O-]", "4-nitramidobenzoic acid", id="nitramido_prefix_under_an_acid"),
        pytest.param("O=C(C=C)NCN([N+](=O)[O-])CN([N+](=O)[O-])CNC(C)=O", "N-[({[(acetamidomethyl)nitramido]methyl}nitramido)methyl]prop-2-enamide", id="substituted_nitramido_prefix"),
        pytest.param("OC(=O)c1ccc(cc1)N(C)[N+](=O)[O-]", "4-(methylnitramido)benzoic acid", id="substituted_nitramide_nitrogen_is_a_substituted_nitramido"),
        pytest.param("O=NN1CCCC1", "1-nitrosopyrrolidine", id="nitroso_on_a_ring_nitrogen"),
        pytest.param("OC(=O)CN(N=O)CC(=O)O", "2,2'-(nitrosoazanediyl)diacetic acid", id="nitroso_on_a_linking_nitrogen"),
        pytest.param("OC(=O)CCNN[N+](=O)[O-]", "3-(2-nitrohydrazin-1-yl)propanoic acid", id="nitrohydrazinyl_prefix"),
        pytest.param("OC(=O)CCN(N)N(=O)=O", "3-(1-nitrohydrazin-1-yl)propanoic acid", id="nitro_on_the_attached_hydrazine_nitrogen"),
        pytest.param("OC(=O)CCNNN=O", "3-(2-nitrosohydrazin-1-yl)propanoic acid", id="nitrosohydrazinyl_prefix"),
        pytest.param("OC(=O)CCNNN", "3-(triazan-1-yl)propanoic acid", id="triazanyl_prefix_unchanged"),
    ],
)
def test_nitroso_and_nitro_groups_on_nitrogen_and_on_rings_with_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OB(O)c1ccc([N+](=O)[O-])cc1", "(4-nitrophenyl)boronic acid", id="boronic_acid"),
        pytest.param("CP(=O)(O)c1ccccc1[N+](=O)[O-]", "methyl(2-nitrophenyl)phosphinic acid", id="phosphinic_acid"),
        pytest.param("N#CNc1cccc([N+](=O)[O-])c1", "(3-nitrophenyl)cyanamide", id="cyanamide"),
        pytest.param("CNC(=S)Nc1cccc([N+](=O)[O-])c1", "N-methyl-N'-(3-nitrophenyl)thiourea", id="thiourea"),
        pytest.param("CC(=O)NNC(=O)c1cccc([N+](=O)[O-])c1", "N'-acetyl-3-nitrobenzohydrazide", id="diacylhydrazine"),
        pytest.param("CS(=O)(=O)NC(=O)c1cccc([N+](=O)[O-])c1", "N-(methanesulfonyl)-3-nitrobenzamide", id="acylsulfonamide"),
        pytest.param("CC(=O)N(C)[N+](=O)[O-]", "N-methyl-N-nitroacetamide", id="nitro_on_amide_nitrogen"),
        pytest.param("COC(=O)N[N+](=O)[O-]", "methyl nitrocarbamate", id="nitro_on_carbamate_nitrogen"),
        pytest.param("O=[N+]([O-])N=C1NCCS1", "N-nitro-1,3-thiazolidin-2-imine", id="nitro_on_ring_imine_nitrogen"),
        pytest.param("C[NH2+]c1ccc([N+](=O)[O-])cc1", "N-methyl-4-nitroanilinium", id="anilinium"),
        pytest.param("C[N+](C)(C)c1cccc([N+](=O)[O-])c1", "N,N,N-trimethyl-3-nitroanilinium", id="quaternary_anilinium"),
        pytest.param("O=NN1CCc2ccccc2C1", "2-nitroso-1,2,3,4-tetrahydroisoquinoline", id="nitroso_on_a_hydrogenated_ring_amine"),
    ],
)
def test_nitro_and_nitroso_prefixes_in_parent_families_with_their_own_module(smiles, expected):
    # P-61.5.1, P-61.5.2: a nitro or nitroso group is only ever a prefix, whatever the parent family
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NC(=O)N=O", "nitrosoformamide", id="amide_of_nitrosoformic_acid"),
        pytest.param("CNC(=O)N=O", "N-methylnitrosoformamide", id="n_substituted_amide_of_nitrosoformic_acid"),
        pytest.param("ClC(=O)N=O", "nitrosoformyl chloride", id="acyl_halide_of_nitrosoformic_acid"),
        pytest.param("OC(=O)c1ccc(C(=O)N=O)cc1", "4-(nitrosoformyl)benzoic acid", id="nitrosoformyl_prefix"),
    ],
)
def test_nitroso_on_a_carbonyl_carbon_is_a_substituent_of_formic_acid(smiles, expected):
    # P-65.1.8.2, P-65.1.8.3: the hydrogen of formic acid and of formyl is substitutable by nitroso; it is not a halide
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(=O)N=O", id="acyl_nitroso"),
        pytest.param("CC(=O)[N+](=O)[O-]", id="acyl_nitro"),
        pytest.param("CCC(=O)C(=O)N=O", id="acyl_nitroso_beside_a_ketone"),
    ],
)
def test_acyl_group_on_nitroso_or_nitro_nitrogen_is_not_an_acyl_halide_or_aldehyde(smiles):
    # P-65.5.1, P-65.5.2.1: neither N=O nor NO2 is a halide or pseudohalide; the carbonyl carbon has no hydrogen
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)NOC", "N-methoxyacetamide", id="n_alkoxy_amide"),
        pytest.param("CC(=O)NOCC(C)C", "N-(2-methylpropoxy)acetamide", id="branched_alkoxy_on_amide_nitrogen"),
        pytest.param("CC(=O)NOc1ccccc1", "N-phenoxyacetamide", id="n_aryloxy_amide"),
        pytest.param("CC(=O)N(C)OC", "N-methoxy-N-methylacetamide", id="alkoxy_and_alkyl_on_one_amide_nitrogen"),
        pytest.param("O=C(NOC)CCC(=O)NOC", "N1,N4-dimethoxybutanediamide", id="n_alkoxy_groups_of_a_diamide"),
        pytest.param("CC(=O)NOC(C)C", "N-[(propan-2-yl)oxy]acetamide", id="compound_organyloxy_prefix_is_enclosed"),
        pytest.param("CC(=O)N(C)SC", "N-methyl-N-(methylsulfanyl)acetamide", id="sulfanyl_and_alkyl_on_one_amide_nitrogen"),
        pytest.param("O=C(c1ccccc1)NSSc1ccccn1", "N-[(pyridin-2-yl)disulfanyl]benzamide", id="chalcogen_chain_on_amide_nitrogen"),
        pytest.param("CS(=O)(=O)NOC", "N-methoxymethanesulfonamide", id="n_alkoxy_sulfonamide"),
    ],
)
def test_n_alkoxy_and_n_aryloxy_amides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)NC(=NC)SC", "methyl N-acetyl-N'-methylcarbamimidothioate", id="acyl_isothiourea_is_an_ester"),
        pytest.param("CC(=O)NC(=NC)NC", "N-(N,N'-dimethylcarbamimidoyl)acetamide", id="substituted_carbamimidoyl_on_amide_nitrogen"),
        pytest.param("CC(=O)NC(=N)NC(C)=O", "N-(N-acetylcarbamimidoyl)acetamide", id="second_acyl_group_on_the_carbamimidoyl_prefix"),
        pytest.param("ClCC(=O)NC(=N)NC(=O)CCl", "2-chloro-N-[N-(chloroacetyl)carbamimidoyl]acetamide", id="second_acyl_group_with_a_substituent"),
        pytest.param("O=C(NC(=N)NC(=O)c1ccccc1)c1ccccc1", "N-(N-benzoylcarbamimidoyl)benzamide", id="identical_ring_amides_on_a_guanidine"),
        pytest.param("CC(=O)NC(=N)C", "N-(ethanimidoyl)acetamide", id="imidoyl_on_amide_nitrogen"),
        pytest.param("CC(=O)N=C(N)N", "N-(diaminomethylidene)acetamide", id="acyl_guanidine_is_an_amide"),
        pytest.param("CC(=O)N=C(C)N", "N-(1-aminoethylidene)acetamide", id="amidine_ylidene_chain_end_on_amide_nitrogen"),
        pytest.param("CC(N)=NS(C)(=O)=O", "N-(1-aminoethylidene)methanesulfonamide", id="amidine_ylidene_chain_end_on_sulfonamide_nitrogen"),
        pytest.param("CC(N)=NCC=N", "N'-(2-iminoethyl)ethanimidamide", id="imino_prefix_beneath_an_amidine"),
        pytest.param("CC(=O)N=C1SC=CN1C", "N-(3-methyl-1,3-thiazol-2(3H)-ylidene)acetamide", id="ring_ylidene_on_amide_nitrogen"),
    ],
)
def test_amide_nitrogen_bonded_to_an_amidine_or_ylidene_carbon(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("COC(=O)NC(C)C(=O)C", "methyl (3-oxobutan-2-yl)carbamate", id="ketone_in_the_n_substituent_is_a_prefix"),
        pytest.param("COC(=O)NC(C)(C(=O)CCl)c1ccccc1", "methyl (4-chloro-3-oxo-2-phenylbutan-2-yl)carbamate", id="ketone_halogen_and_phenyl_in_the_n_substituent"),
        pytest.param(
            "CC(C)(C)OC(=O)N[C@@H](CO)C1(O)CC1",
            "tert-butyl [(1S)-2-hydroxy-1-(1-hydroxycyclopropyl)ethyl]carbamate",
            id="stereocentre_of_the_n_substituent_is_cited_inside_it",
        ),
        pytest.param(
            "CCC[C@H](C[C@@H]1CCCO1)NC(=O)OC",
            "methyl {(2R)-1-[(2S)-oxolan-2-yl]pentan-2-yl}carbamate",
            id="stereodescriptors_in_a_substituent_with_a_ring_prefix",
        ),
    ],
)
def test_carbamate_esters_with_groups_in_the_n_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)n1cccc1", "1-(1H-pyrrol-1-yl)ethan-1-one", id="acyl_on_the_nitrogen_of_a_mancude_monocycle"),
        pytest.param("CC(=O)n1ccc2ccccc21", "1-(1H-indol-1-yl)ethan-1-one", id="acyl_on_the_nitrogen_of_a_fused_mancude_system"),
        pytest.param("CC(=O)n1nnc2ccccc21", "1-(1H-1,2,3-benzotriazol-1-yl)ethan-1-one", id="acyl_on_a_nitrogen_beside_another_ring_nitrogen"),
        pytest.param("CC(=O)n1nnc2ncccc21", "1-(1H-[1,2,3]triazolo[4,5-b]pyridin-1-yl)ethan-1-one", id="acyl_on_a_fused_heterocycle_with_bracketed_locants"),
        pytest.param("O=CN1CCc2ccccc2C1", "3,4-dihydroisoquinoline-2(1H)-carbaldehyde", id="formyl_on_a_ring_nitrogen_is_a_carbaldehyde"),
    ],
)
def test_n_acyl_groups_on_mancude_and_fused_ring_nitrogens_are_pseudoketones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C(n1ccnc1)N1CCCC1", "(1H-imidazol-1-yl)(pyrrolidin-1-yl)methanone", id="carbonyl_between_two_ring_nitrogens_is_a_methanone"),
        pytest.param("S=C(n1cncc1)n1cncc1", "di(1H-imidazol-1-yl)methanethione", id="thiocarbonyl_between_two_ring_nitrogens_is_a_methanethione"),
        pytest.param("O=C1C=CC=CN1C(=O)N1CCCC1", "1-(pyrrolidine-1-carbonyl)pyridin-2(1H)-one", id="ring_ketone_beside_a_hidden_amide_pair_takes_the_ring_as_parent"),
    ],
)
def test_carbonyl_and_thiocarbonyl_between_ring_nitrogens_are_not_ureas(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N#CCc1ccccc1", "phenylacetonitrile"),
        ("OCC=O", "hydroxyacetaldehyde"),
        ("C1OC1c1ccccc1", "phenyloxirane"),
        ("CC1CN1", "2-methylaziridine"),
        ("ClCC(N)=O", "2-chloroacetamide"),
    ],
)
def test_locant_is_omitted_only_when_the_parent_has_one_kind_of_substitutable_hydrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H](CCc1ccccc1)NC(=S)NC1CC1", "N-cyclopropyl-N'-[(2R)-4-phenylbutan-2-yl]thiourea"),
        ("C[C@H](CCc1ccccc1)NC(=O)NC", "N-methyl-N'-[(2R)-4-phenylbutan-2-yl]urea"),
        ("C[C@@H](CCc1ccccc1)NC(=S)NC", "N-methyl-N'-[(2S)-4-phenylbutan-2-yl]thiourea"),
        (
            "C[C@H](CCc1ccccc1)NC(=S)N[C@@H]1C[C@H]2CC[C@H]1C2",
            "N-[(1S,2R,4S)-bicyclo[2.2.1]heptan-2-yl]-N'-[(2R)-4-phenylbutan-2-yl]thiourea",
        ),
    ],
)
def test_chiral_substituent_groups_of_ureas_and_thioureas_cite_their_descriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("NC(=O)ONC", "N-(carbamoyloxy)methanamine"),
        ("NC(=O)ON", "hydroxylamine-O-carboxamide"),
        ("NC(=O)ON(C)C", "N-(carbamoyloxy)-N-methylmethanamine"),
        ("CNC(=O)ONC", "N-[(methylcarbamoyl)oxy]methanamine"),
    ],
)
def test_carbamic_acid_esters_of_amino_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("ON[SiH3]", "N-hydroxysilanamine", id="n_hydroxy_silanamine"),
        pytest.param("ONB", "N-hydroxyboranamine", id="n_hydroxy_boranamine"),
        pytest.param("NOCCl", "O-(chloromethyl)hydroxylamine", id="o_halogenated_alkyl_hydroxylamine"),
        pytest.param("NOCC1CCCCC1", "O-(cyclohexylmethyl)hydroxylamine", id="o_ring_alkyl_hydroxylamine"),
        pytest.param("NOc1ccccc1", "O-phenylhydroxylamine", id="o_phenyl_hydroxylamine"),
        pytest.param("C=CON", "O-ethenylhydroxylamine", id="o_ethenyl_hydroxylamine"),
        pytest.param("NOCc1ccccc1", "O-benzylhydroxylamine", id="o_benzyl_hydroxylamine"),
        pytest.param("NOc1ccccc1C", "O-(2-methylphenyl)hydroxylamine", id="o_substituted_phenyl_hydroxylamine"),
        pytest.param("NOC(C)C", "O-(propan-2-yl)hydroxylamine", id="o_branched_alkyl_hydroxylamine"),
        pytest.param("c1ccccc1/C=C/CON", "O-[(2E)-3-phenylprop-2-en-1-yl]hydroxylamine", id="o_alkenyl_with_a_ring_cites_its_configuration"),
        pytest.param("NOC[C@H](C)CC", "O-[(2R)-2-methylbutyl]hydroxylamine", id="o_alkyl_cites_its_stereocentre"),
    ],
)
def test_hydroxylamines_with_a_hydride_or_ring_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C(C)N(C(=NC)C1=CC(=CC2=CC=CC=C12)C(N(C)C)=NCC)CC", "N1,N1,N'3-triethyl-N'1,N3,N3-trimethylnaphthalene-1,3-dicarboximidamide"),
        ("C(C)N=C(N)C1(CCCCC1)C(N(C)C)=N", "N'''1-ethyl-N1,N1-dimethylcyclohexane-1,1-dicarboximidamide"),
        ("C(C)NC(=N)C1(CCCCC1)C(N(C)C)=N", "N''1-ethyl-N1,N1-dimethylcyclohexane-1,1-dicarboximidamide"),
        ("CN(C)C(=N)CCC(=N)NCC", "N4-ethyl-N1,N1-dimethylbutanediimidamide"),
        ("CCN=C(N(C)C)c1ccc(C(=O)O)cc1", "4-(N'-ethyl-N,N-dimethylcarbamimidoyl)benzoic acid"),
        ("CN(c1ccc(C(=O)O)cc1)C(=N)N(C)c1ccccc1", "4-[methyl(N-methyl-N-phenylcarbamimidoyl)amino]benzoic acid"),
        ("N=CNc1ccc(C(=O)O)cc1", "4-methanimidamidobenzoic acid"),
        ("CC(=N)Nc1ccc(C(=O)O)cc1", "4-ethanimidamidobenzoic acid"),
        ("NS(=N)CCC(=O)O", "3-(S-aminosulfinimidoyl)propanoic acid"),
        ("N[Se](=N)(=O)CCC(=O)O", "3-(Se-aminoselenonimidoyl)propanoic acid"),        ("NC(N)=NCCCC(=O)O", "4-[(diaminomethylidene)amino]butanoic acid"),
        ("C(N)(=N)NC=O", "N-carbamimidoylformamide"),
        ("C(N)(=N)NC(C)=O", "N-carbamimidoylacetamide"),
        ("C(N)(=N)NC(N)=O", "N-carbamimidoylurea"),
        ("C(N)(=N)NC(N)=N", "imidodicarbonimidic diamide"),
        ("C(N)(=N)NC(=N)NC(N)=N", "diimidotricarbonimidic diamide"),
        ("C(C)N=C(N(C1=CC=CC=C1)C1=CC=CC=C1)NC(N)=N", "N'1-ethyl-N1,N1-diphenylimidodicarbonimidic diamide"),
        ("CN=C(N)N=C(N)N", "N'1-methylimidodicarbonimidic diamide"),
        ("N=C(N)N=C(N)N(C)C", "N1,N1-dimethylimidodicarbonimidic diamide"),
        ("CN=C(N)NC(N)=NC", "N'1,N'3-dimethylimidodicarbonimidic diamide"),
        ("N=C(N)N/C(N)=N/C", "(1E)-N'1-methylimidodicarbonimidic diamide"),
        ("C/N=C(\\N)N/C(N)=N\\C", "(1Z,3E)-N'1,N'3-dimethylimidodicarbonimidic diamide"),
        ("N=C(N)NC(=NC)NC(N)=N", "N3-methyldiimidotricarbonimidic diamide"),
        ("CN=C(N)NC(=N)NC(N)=NC", "N'1,N'5-dimethyldiimidotricarbonimidic diamide"),
        ("N=C(N)NC(=N)NC(C)=O", "N-(N-carbamimidoylcarbamimidoyl)acetamide"),
        ("CN=C(N)NC(N)=NCC(=O)O", "({amino[(N'-methylcarbamimidoyl)amino]methylidene}amino)acetic acid"),
        ("N=C(NC(N)=N)NC(NC(NC(N)=N)=N)=N", "3,5,7-triimino-2,4,6,8-tetraazanonane-1,9-diimidamide"),
        ("NN=C(C)NCCC(=O)O", "3-(ethanehydrazonamido)propanoic acid"),
        ("NN=S(c1ccccc1)Nc1ccc(C(=O)O)cc1", "4-(benzenesulfinohydrazonamido)benzoic acid"),
        ("NN=CNc1ccc(C(=O)O)cc1", "4-(methanehydrazonamido)benzoic acid"),
        ("NNC=Nc1ccc(C(=O)O)cc1", "4-[(hydrazinylmethylidene)amino]benzoic acid"),
        ("N=CNNCCC(=O)O", "3-(methanimidohydrazido)propanoic acid"),
        ("CC(=N)NNCCC(=O)OC", "methyl 3-(ethanimidohydrazido)propanoate"),
        ("N=C(NNc1ccc(C(=O)O)cc1)c1ccccc1", "4-(benzenecarboximidohydrazido)benzoic acid"),
        ("NNC(=N)c1cccc(C(=O)O)c1", "3-(hydrazinecarboximidoyl)benzoic acid"),
        ("OC(=O)CCNNC=NN", "3-[2-(methanehydrazonoyl)hydrazin-1-yl]propanoic acid"),
        ("OC(=O)CCNNC(=NN)c1ccccc1", "3-[2-(benzenecarbohydrazonoyl)hydrazin-1-yl]propanoic acid"),
        ("OC(=O)c1ccc(cc1)NNC(=N)OC(N)=N", "4-[C-(carbamimidoyloxy)methanimidohydrazido]benzoic acid"),
    ],
)
def test_amidine_locants_prefixes_and_condensed_guanidines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N1(CCCCC1)C(=O)NN", "piperidine-1-carbohydrazide"),
        ("C1COCCN1C(=O)NN", "morpholine-4-carbohydrazide"),
        ("CN(N)C(=O)C1=CC=C(C2=CC=CC=C12)C(=O)NNC", "N1,N'4-dimethylnaphthalene-1,4-dicarbohydrazide"),
        ("CN(N)C(=O)C1=CC=C(C2=CC(=CC=C12)C)C(=O)NNC", "N1,N'4,6-trimethylnaphthalene-1,4-dicarbohydrazide"),
        ("C(C1=CC=CC=C1)(=O)NNC(C1=CC=CC=C1)=O", "N'-benzoylbenzohydrazide"),
        ("C(C)(=O)N(N(C(CC)=O)C)CC", "N'-acetyl-N'-ethyl-N-methylpropanehydrazide"),
        ("CCCC(=O)NNC(=O)C", "N'-acetylbutanehydrazide"),
        ("CC(=O)NNC(=O)c1ccccc1", "N'-acetylbenzohydrazide"),
        ("CC(C)C(=O)NNC(=O)C(C)C", "2-methyl-N'-(2-methylpropanoyl)propanehydrazide"),
        ("CC(=O)NNC(=O)CCC(=O)NNC(C)=O", "N'1,N'4-diacetylbutanedihydrazide"),
        ("CC(=O)NNC(=O)c1ccc(cc1)C(=O)NNC(C)=O", "N'1,N'4-diacetylbenzene-1,4-dicarbohydrazide"),
        ("CC(=O)NNC(=O)CCC(=O)NNC(=O)c1ccccc1", "N'1-acetyl-N'4-benzoylbutanedihydrazide"),
        ("C(CC)(NN)=S", "propanethiohydrazide"),
        ("C1(=CC=CC=C1)C(NN)=S", "benzenecarbothiohydrazide"),
        ("NNC(=S)CCC(=S)NN", "butanedithiohydrazide"),
        ("NNC(=S)C(=S)NN", "ethanedithiohydrazide"),
        ("NNC(=S)N1CCCCC1", "piperidine-1-carbothiohydrazide"),
        ("CC(=[Se])NN", "ethaneselenohydrazide"),
        ("NNC(=O)NN", "hydrazinecarbohydrazide"),
        ("C(=O)(OC(=O)NN)NN", "dicarbonic dihydrazide"),
        ("C(=O)(NN)NC(=O)NN", "2-imidodicarbonic dihydrazide"),
        ("NNC(=O)NNCCC(=O)O", "3-(hydrazinecarbohydrazido)propanoic acid"),
        ("O=CNNCCC(=O)O", "3-formohydrazidopropanoic acid"),
        ("N(N)S(=O)C=1C(=CC2=CC=CC=C2C1)C(=O)O", "3-(hydrazinesulfinyl)naphthalene-2-carboxylic acid"),
        ("C(C1=CC=CC=C1)=NNC(C(=O)N)=O", "2-(benzylidenehydrazinyl)-2-oxoacetamide"),
    ],
)
def test_hydrazide_variants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CN=C(NC)N1CCCCC1", "N,N'-dimethylpiperidine-1-carboximidamide", id="amidine_group_on_a_ring_nitrogen"),
        pytest.param("O=C(NNC)N1CCCCC1", "N'-methylpiperidine-1-carbohydrazide", id="hydrazide_group_on_a_ring_nitrogen"),
        pytest.param("C1CCCCN1NC(C)=O", "N-(piperidin-1-yl)acetamide", id="amide_nitrogen_on_a_ring_nitrogen"),
        pytest.param("C1CCCCN1NS(=O)(=O)C", "N-(piperidin-1-yl)methanesulfonamide", id="sulfonamide_nitrogen_on_a_ring_nitrogen"),
        pytest.param("CC(=O)Nn1cccn1", "N-(1H-pyrazol-1-yl)acetamide", id="amide_nitrogen_on_a_pyrazole_nitrogen"),
        pytest.param("CC(=O)Nn1ccc(C)n1", "N-(3-methyl-1H-pyrazol-1-yl)acetamide", id="substituted_azole_ring_group"),
        pytest.param("CC(=O)Nn1cncn1", "N-(1H-1,2,4-triazol-1-yl)acetamide", id="amide_nitrogen_on_a_triazole_nitrogen"),
        pytest.param("CC(=O)Nn1nnnc1", "N-(1H-tetrazol-1-yl)acetamide", id="amide_nitrogen_on_a_tetrazole_nitrogen"),
        pytest.param("CC(=O)N(C)n1cccn1", "N-methyl-N-(1H-pyrazol-1-yl)acetamide", id="n_methyl_beside_the_azole_group"),
        pytest.param("CS(=O)(=O)Nn1cccn1", "N-(1H-pyrazol-1-yl)methanesulfonamide", id="sulfonamide_nitrogen_on_a_pyrazole_nitrogen"),
        pytest.param("CC(=O)NN1CCCO1", "N-(1,2-oxazolidin-2-yl)acetamide", id="amide_nitrogen_on_a_saturated_n_o_ring_nitrogen"),
        pytest.param("CC(=N)NN1CCCCC1", "N-(piperidin-1-yl)ethanimidamide", id="amidine_nitrogen_on_a_ring_nitrogen"),
        pytest.param("N=C(N)NN1CCCCC1", "N-(piperidin-1-yl)guanidine", id="guanidine_nitrogen_on_a_ring_nitrogen"),
        pytest.param(
            "OC(=O)CC(=O)NN1CCCCC1", "3-oxo-3-[(piperidin-1-yl)amino]propanoic acid", id="amido_prefix_with_a_ring_nitrogen_beside_an_acid"
        ),
    ],
)
def test_group_nitrogen_bonded_to_a_ring_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NC1=CC=C(C=C1)S(=O)(=O)NC=1SC(=CN1)C(=O)O", "2-(4-aminobenzene-1-sulfonamido)-1,3-thiazole-5-carboxylic acid"),
        ("OC(=O)c1cnc(NS(=O)(=O)C)s1", "2-(methanesulfonamido)-1,3-thiazole-5-carboxylic acid"),
        ("OC(=O)c1cnc(NC(C)=O)s1", "2-acetamido-1,3-thiazole-5-carboxylic acid"),
    ],
)
def test_acylamino_prefix_on_heteroaromatic_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCOC=NCC(=O)O", "[(ethoxymethylidene)amino]acetic acid"),
        ("OC(=O)c1ccc(cc1)N=C(OC)SC", "4-{[methoxy(methylsulfanyl)methylidene]amino}benzoic acid"),
        (
            "CN(C)CCOC=NC1=CC=C2C=CC(=CC2=C1)C(=O)O",
            "7-({[2-(dimethylamino)ethoxy]methylidene}amino)naphthalene-2-carboxylic acid",
        ),
    ],
)
def test_imidate_prefix_on_ring_and_chain_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NNC(=S)C(=O)NN", "2-hydrazinyl-2-sulfanylideneacetohydrazide"),
        ("NNC(=S)CC(=O)NN", "3-hydrazinyl-3-sulfanylidenepropanehydrazide"),
        ("NNC(=S)CCCC(=O)NN", "5-hydrazinyl-5-sulfanylidenepentanehydrazide"),
    ],
)
def test_chalcogen_hydrazide_end_group_beside_a_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NNC(=NN)C", "ethanehydrazonohydrazide"),
        ("CNNC(=NNC)c1cccs1", "N',N''-dimethylthiophene-2-carbohydrazonohydrazide"),
        ("NNC(=NN)c1cccc(c1)C(=O)O", "3-(hydrazinecarbohydrazonoyl)benzoic acid"),
        ("NNC(=NN)CC(=O)O", "3-hydrazinyl-3-hydrazinylidenepropanoic acid"),
        ("NNS(=NN)c1ccccc1", "benzenesulfinohydrazonohydrazide"),
        ("NNS(=O)(=NN)C", "methanesulfonohydrazonohydrazide"),
        ("NNC(=NN)NN", "hydrazinecarbohydrazonohydrazide"),
        ("NNC(=N)NN", "hydrazinecarboximidohydrazide"),
        ("NC(=NN)N", "carbonohydrazonic diamide"),
        ("NC(=NN)OC(=NN)N", "dicarbonohydrazonic diamide"),
    ],
)
def test_hydrazidines_and_carbonic_amidrazones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=N)NC(=N)C", "N-ethanimidoylethanimidamide"),
        ("CC(N)=NC(C)=N", "N'-ethanimidoylethanimidamide"),
        ("CCC(=N)NC(=N)C", "N-ethanimidoylpropanimidamide"),
        ("N=C(c1ccccc1)NC(=N)c1ccccc1", "N-benzenecarboximidoylbenzenecarboximidamide"),
        ("CCNC(=N)SSC(=N)NC", "N-ethylcarbamimidic N-methylcarbamimidic dithioperoxyanhydride"),
        ("[SiH2]([SiH2]C(N)=N)C(N)=N", "disilane-1,2-dicarboximidamide"),
        ("NC(=O)[SiH2][SiH2]C(N)=O", "disilane-1,2-dicarboxamide"),
    ],
)
def test_diamidides_formamidine_disulfides_and_amidines_on_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CN(N=NN(C)C)C=O", "1,4,4-trimethyltetraaz-2-ene-1-carbaldehyde"),
        ("O=CN=NN=NC", "4-methyltetraaza-1,3-diene-1-carbaldehyde"),
    ],
)
def test_carbaldehyde_on_a_polyaza_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C(CON)ON", "O,O'-(ethane-1,2-diyl)bis(hydroxylamine)"),
        ("NSC", "S-methyl(thiohydroxylamine)"),
        ("NS", "thiohydroxylamine"),
        ("CS(=O)ON", "O-(methanesulfinyl)hydroxylamine"),
        ("NOS(=O)(=O)O", "hydroxylamine-O-sulfonic acid"),
        ("NOC(=O)O", "hydroxylamine-O-carboxylic acid"),
        ("NONC", "N-(aminooxy)methanamine"),
        ("NOC(=O)c1ccccc1", "aminooxy(phenyl)methanone"),
        ("CNS", "N-sulfanylmethanamine"),
        ("CC(=O)NS", "N-sulfanylacetamide"),
        ("C[NH][SeH]", "N-selanylmethanamine"),
        ("CN([TeH])C", "N-methyl-N-tellanylmethanamine"),
        ("ON(C=1C=C(C(=O)O)C=CC1)C=1C=C(C(=O)O)C=CC1", "3,3'-(hydroxyazanediyl)dibenzoic acid"),
    ],
)
def test_hydroxylamine_functional_parent_and_chalcogen_analogues(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CN=NC=O", "diazenedicarbaldehyde"),
        ("N#CN=Nc1ccccc1", "phenyldiazenecarbonitrile"),
        ("OC(=O)N=NC(=O)O", "diazenedicarboxylic acid"),
        ("N=Nc1ccc(N=N)c(c1)C(=O)O", "2,5-bis(diazenyl)benzoic acid"),
        ("OC(=O)CCNNC(N)=O", "3-(2-carbamoylhydrazin-1-yl)propanoic acid"),
        ("CC(=NO)C(C)=NO", "N2,N3-dihydroxybutane-2,3-diimine"),
        ("CC(=NOc1ccccc1)CC(C)=NOc1ccccc1", "N2,N4-diphenoxypentane-2,4-diimine"),
        ("CCC(=NO)[N+](=O)[O-]", "N-hydroxy-1-nitropropan-1-imine"),
        ("CC(=NO)N=O", "N-hydroxy-1-nitrosoethan-1-imine"),
        ("CCON=Cc1ccc(cc1)S(=O)(=O)O", "4-[(ethoxyimino)methyl]benzene-1-sulfonic acid"),
        ("CC(=O)C(N=Nc1ccccc1)=NNc1ccccc1", "1-(phenyldiazenyl)-1-(phenylhydrazinylidene)propan-2-one"),
        ("CC(=O)NN=CN=Nc1ccccc1", "N'-[(phenyldiazenyl)methylidene]acetohydrazide"),
        ("CC(=O)N(C)/N=N\\C", "N-methyl-N-[(Z)-methyldiazenyl]acetamide"),
        ("CC(=O)CN(C)/N=N/C", "1-[(2E)-1,3-dimethyltriaz-2-en-1-yl]propan-2-one"),
        ("NC(=O)N=N", "diazenecarboxamide"),
        ("N=NC(=O)NN", "diazenecarbohydrazide"),
        ("S=C(NNc1ccccc1)N=Nc1ccccc1", "N',2-diphenyldiazenecarbothiohydrazide"),
        ("O=C(NNc1ccccc1)N=Nc1ccccc1", "N',2-diphenyldiazenecarbohydrazide"),
        ("CCOC(=O)CCNNC(=O)N=N", "ethyl 3-(diazenecarbohydrazido)propanoate"),
        ("c1ccc(cc1)N=Nc1ccc2ccc(N=Nc3ccc4cc5ccccc5cc4c3)cc2c1", "{7-[(anthracen-2-yl)diazenyl]naphthalen-2-yl}(phenyl)diazene"),
        ("Oc1c(N=Nc2ccccc2)ccc2ccc(N=Nc3ccccc3)c(O)c12", "2,7-bis(phenyldiazenyl)naphthalene-1,8-diol"),
        ("CCN=NC(CC)=NC", "1-(ethyldiazenyl)-N-methylpropan-1-imine"),
        ("CC(=N)N=NC(C)=N", "1,1'-diazenediyldi(ethan-1-imine)"),
    ],
)
def test_diazene_carbo_suffixes_and_substituted_hydrazinyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("NOC(N)=O", "hydroxylamine-O-carboxamide"),
        ("NS", "thiohydroxylamine"),
        ("CSN", "S-methyl(thiohydroxylamine)"),
    ],
)
def test_hydroxylamine_acid_forms_and_chalcogen_analogues(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereogenic_nitrogen_of_an_amine_oxide():
    assert smiles_to_iupac("C(C)[N@+](C1=CC=CC=C1)(C)[O-]") == "(S)-N-ethyl-N-methylaniline N-oxide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N=NC=NN", "formazan"),
        ("c1ccccc1N=NC(c1ccccc1)=NNc1ccccc1", "1,3,5-triphenylformazan"),
        ("N=NC(c1ccccc1)=NNc1ccccc1", "3,5-diphenylformazan"),
        ("OC(=O)C(N=Nc1ccccc1)=NNc1ccccc1", "1,5-diphenylformazan-3-carboxylic acid"),
    ],
)
def test_formazan_and_its_hydrocarbyl_derivatives(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)c1ccc(cc1)N=NC=NN", "4-(formazan-1-yl)benzoic acid", id="1-yl"),
        pytest.param("OC(=O)c1ccc(cc1)NN=CN=N", "4-(formazan-5-yl)benzoic acid", id="5-yl"),
        pytest.param("OC(=O)c1ccc(cc1)C(N=N)=NN", "4-(formazan-3-yl)benzoic acid", id="3-yl"),
        pytest.param("OC(=O)c1ccc(cc1)C(=NNc1ccccc1)N=Nc1ccccc1", "4-(1,5-diphenylformazan-3-yl)benzoic acid", id="substituted_3-yl"),
        pytest.param("Nc1ccc(cc1)N=NC=NN", "4-(formazan-1-yl)aniline", id="on_an_amine"),
        pytest.param("OC(=O)CC(=NNc1ccccc1)N=Nc1ccccc1", "3-(phenyldiazenyl)-3-(phenylhydrazinylidene)propanoic acid", id="chain_parent_includes_C3"),
        pytest.param("OC(=O)c1ccc(cc1)N=NC=NNc1ccc(cc1)C(=O)O", "4,4'-(formazan-1,5-diyl)dibenzoic acid", id="1,5-diyl"),
        pytest.param("OC(=O)c1ccc(cc1)C(N=N)=NNc1ccc(cc1)C(=O)O", "4,4'-(formazan-3,5-diyl)dibenzoic acid", id="3,5-diyl"),
        pytest.param(
            "N#CC(=NNc1cc(ccc1O)S(O)(=O)=O)N=Nc1cc(ccc1O)S(O)(=O)=O",
            "3,3'-(3-cyanoformazan-1,5-diyl)bis(4-hydroxybenzene-1-sulfonic acid)",
            id="substituted_1,5-diyl",
        ),
        pytest.param(
            "OC(=O)c1ccc(cc1)N=NC(c1ccc(cc1)C(=O)O)=NNc1ccc(cc1)C(=O)O", "4,4',4''-(formazan-1,3,5-triyl)tribenzoic acid", id="1,3,5-triyl"
        ),
    ],
)
def test_formazan_substituent_groups_and_linkers_are_retained_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCC([N+]#[C-])CC1", "isocyanocyclohexane"),
        ("OC(=O)c1ccc(cc1)[N+]#[C-]", "4-isocyanobenzoic acid"),
        ("CC(C)(C)N=C(C(=O)O)C(C)(O[N+](=O)[O-])C", "2-(tert-butylimino)-3-methyl-3-(nitrooxy)butanoic acid"),
        ("C1CCC(CC1)N=C=NC1CCCCC1", "dicyclohexylmethanediimine"),
        ("CN=C=NCC", "N-ethyl-N'-methylmethanediimine"),
        ("N=C=NC", "N-methylmethanediimine"),
        ("C[Si](C)(C)N=C=NC", "N-methyl-N'-(trimethylsilyl)methanediimine"),
        ("C[Si](C)(C)N=C=N[Si](C)(C)C", "bis(trimethylsilyl)methanediimine"),
        ("ON=C=NO", "dihydroxymethanediimine"),
        ("CC(=O)N=C=NC(C)=O", "N,N'-(methanediylidene)diacetamide"),
        ("CN=C=NC(=O)c1ccccc1", "N-[(methylimino)methylidene]benzamide"),
        ("CN=C=NCC(=O)O", "{[(methylimino)methylidene]amino}acetic acid"),
        ("CN=C=NC#N", "[(methylimino)methylidene]cyanamide"),
        ("CP=N", "1-methylphosphanimine"),
        ("CN=P", "N-methylphosphanimine"),
        ("N=[AsH]", "arsanimine"),
        ("N=S", "sulfanimine"),
        ("N=B", "boranimine"),
        ("CP=[SiH2]", "methyl(silylidene)phosphane"),
        ("N#P", "phosphanylidyneazane"),
        ("c1ccccc1N=[Se]", "N-phenylselanimine"),
        ("OCCN=S", "2-(sulfanylideneamino)ethan-1-ol"),
        ("C[Si](C)=Nc1ccccc1", "1,1-dimethyl-N-phenylsilanimine"),
        ("S1C(CCC1)=N", "thiolan-2-imine"),
        ("CN=C1C=CC(C2=CC=CC=C12)=NC", "N1,N4-dimethylnaphthalene-1,4-diimine"),
        ("[SiH]#[SiH]", "disilyne"),
        ("NN=NNN", "pentaaz-2-ene"),
        ("c1ccccc1N=NNc1ccccc1", "1,3-diphenyltriaz-1-ene"),
        ("[SnH2]=P[SnH]=P[SnH3]", "tristannaphospha-1,3-diene"),
        ("[SiH2]([SiH3])SS[SiH2][Si](C)(C)C", "2-(disilanyldisulfanyl)-1,1,1-trimethyldisilane"),
        ("CPPN=O", "1-methyl-2-nitrosodiphosphane"),
        ("O=Nc1ccc(cc1)N=O", "1,4-dinitrosobenzene"),
        ("N1=CC#CC=C1", "3,4-didehydropyridine"),
        ("[C-]#[O+]", "carbon monoxide"),
        ("ON[CH]", "N-hydroxy-λ2-methanamine"),
    ],
)
def test_noncarbon_parents_and_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCCCC=NNN=O", "N'-hexylidenenitrous hydrazide"),
        ("CN(C)NN=O", "N',N'-dimethylnitrous hydrazide"),
        ("S=NO", "thionitrous O-acid"),
    ],
)
def test_nitrous_hydrazides_and_chalcogen_nitrous_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NOCNN", "[(aminooxy)methyl]hydrazine"),
        ("c1ccccc1NN=CC=NNc1ccccc1", "1,1′-ethanediylidenebis(2-phenylhydrazine)"),
        ("NNCCNN", "1,1′-ethane-1,2-diyldihydrazine"),
        ("OC(=O)C1CCC(=NN=C2CCC(CC2)C(=O)O)CC1", "4,4'-hydrazinediylidenedi(cyclohexane-1-carboxylic acid)"),
    ],
)
def test_hydrazine_with_aminooxy_prefix_and_multiplied_hydrazine_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[NH][Br+][O-]", "methylbromous amide"),
        ("C[NH][Cl+3]([O-])([O-])[O-]", "methylperchloric amide"),
        ("OC(=O)C1CCC(CC1)=N(=O)O", "4-[hydroxy(oxo)-λ5-azanylidene]cyclohexane-1-carboxylic acid"),
        ("CS[N+](=O)[O-]", "S-methyl thionitrate"),
        ("S=C=N[SiH3]", "isothiocyanatosilane"),
    ],
)
def test_halogen_oxoacid_amides_aci_nitro_prefix_thionitrates_and_pseudohalide_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[N+](C)([O-])CC#N", "cyano-N,N-dimethylmethanamine N-oxide"),
        ("C[N+](C)([O-])CCC#N", "2-cyano-N,N-dimethylethan-1-amine N-oxide"),
    ],
)
def test_amine_oxide_outranks_a_nitrile_cited_as_cyano(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[NH2+][O-]", "ethanamine N-oxide"),
        ("CC[N+](CC)(CC)[S-]", "N,N-diethylethanamine N-sulfide"),
        ("NC(=O)c1cc(CC[NH2+][O-])cc(C(=O)O)c1", "2-(3-carbamoyl-5-carboxyphenyl)ethan-1-amine N-oxide"),
    ],
)
def test_primary_amine_oxides_chalcogen_analogues_and_oxides_beside_senior_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NC(=O)NC(=O)N", "2-imidodicarbonic diamide"),
        ("NC(=O)NC(=O)NC(=O)N", "2,4-diimidotricarbonic diamide"),
        ("CNC(=O)NC(=O)N", "N1-methyl-2-imidodicarbonic diamide"),
        ("CN(C(=O)N)C(=O)NC", "N1,N2-dimethyl-2-imidodicarbonic diamide"),
        ("NC(=O)N(CCO)C(=O)NCC#N", "N1-(cyanomethyl)-N2-(2-hydroxyethyl)-2-imidodicarbonic diamide"),
        ("CN(C(N)=S)C(N)=S", "N2-methyl-2-imido-1,3-dithiodicarbonic diamide"),
        ("NC(=O)N(C)C(=O)N(C)C(N)=O", "N2,N4-dimethyl-2,4-diimidotricarbonic diamide"),
        ("CNC(=S)NC(=O)N", "N1-methyl-2-imido-1-thiodicarbonic diamide"),
        ("CNC(=O)NC(=S)NC(=O)N", "N1-methyl-2,4-diimido-3-thiotricarbonic diamide"),
        ("NC(=O)OC(N)=O", "dicarbonic diamide"),
        ("NC(=O)OC(=O)OC(N)=O", "tricarbonic diamide"),
        ("CC(C)NC(=O)OC(N)=O", "N1-(propan-2-yl)dicarbonic diamide"),
        ("NC(=S)SC(N)=S", "1,2,3-trithiodicarbonic diamide"),
        ("NC(=S)SSC(N)=S", "2-dithioperoxy-1,3-dithiodicarbonic diamide"),
        ("N#CNC(C)C", "(propan-2-yl)cyanamide"),
        ("N#CN(CC)CC", "diethylcyanamide"),
        ("CSCCN(C)C#N", "methyl[2-(methylsulfanyl)ethyl]cyanamide"),
        ("N#CNC1CCCN1", "(pyrrolidin-2-yl)cyanamide"),
        ("C[C@H](Cl)CNC#N", "[(2S)-2-chloropropyl]cyanamide"),
        ("N#CNc1ccco1", "(furan-2-yl)cyanamide"),
        ("N#CNN1CCCCC1", "(piperidin-1-yl)cyanamide"),
        ("N#CNN1CCOCC1", "(morpholin-4-yl)cyanamide"),
        ("N#CN(C)N1CCCC1", "methyl(pyrrolidin-1-yl)cyanamide"),
        ("N#CN(CCO)CCN", "(2-aminoethyl)(2-hydroxyethyl)cyanamide"),
        ("NC(=O)NNC(N)=O", "hydrazine-1,2-dicarboxamide"),
        ("C/C=C\\N(C(N)=O)N(CC)C(N)=O", "1-ethyl-2-[(1Z)-prop-1-en-1-yl]hydrazine-1,2-dicarboxamide"),
        ("CNC(=O)NNC(=O)NC", "N1,N2-dimethylhydrazine-1,2-dicarboxamide"),
        ("NC(=S)NNC(N)=S", "hydrazine-1,2-dicarbothioamide"),
        ("CNC(=S)NNC(=S)NC", "N1,N2-dimethylhydrazine-1,2-dicarbothioamide"),
        ("NC(=S)N(C)NC(N)=S", "1-methylhydrazine-1,2-dicarbothioamide"),
        ("CC(C)NC(=S)NNC(=S)N", "N1-(propan-2-yl)hydrazine-1,2-dicarbothioamide"),
        ("NC(=[Se])NNC(=[Se])N", "hydrazine-1,2-dicarboselenoamide"),
        ("NC(=O)NNC(N)=S", "2-carbamothioylhydrazine-1-carboxamide"),
        ("NC(=S)NNC(N)=O", "2-carbamothioylhydrazine-1-carboxamide"),
        ("NC(=O)NNC(=S)NC", "2-(methylcarbamothioyl)hydrazine-1-carboxamide"),
        ("NC(=S)NNC(=O)NC", "2-carbamothioyl-N-methylhydrazine-1-carboxamide"),
        ("NC(=O)N(C)NC(=S)NC", "1-methyl-2-(methylcarbamothioyl)hydrazine-1-carboxamide"),
        ("NC(=O)NNC(=[Se])N", "2-carbamoselenoylhydrazine-1-carboxamide"),
        ("NC(=O)NN(C(=O)N)C([N+](=O)[O-])[N+](=O)[O-]", "1-(dinitromethyl)hydrazine-1,2-dicarboxamide"),
        ("N#CNN", "cyanohydrazide"),
        ("NC(N)=NN", "carbonohydrazonic diamide"),
        ("CNC(N)=NN", "N-methylcarbonohydrazonic diamide"),
        ("CN(C)C(N)=NN", "N,N-dimethylcarbonohydrazonic diamide"),
        ("CNC(=NN)NC", "N,N'-dimethylcarbonohydrazonic diamide"),
        ("CN(C)C(=NN)N(C)C", "N,N,N',N'-tetramethylcarbonohydrazonic diamide"),
        ("NC(=NN)NC1CC1", "N-cyclopropylcarbonohydrazonic diamide"),
        ("NC(=NN)Nc1ccccc1", "N-phenylcarbonohydrazonic diamide"),
        ("NC(=NN)NOC", "N-methoxycarbonohydrazonic diamide"),
        ("NC(=NN)N[Si](C)(C)C", "N-(trimethylsilyl)carbonohydrazonic diamide"),
        ("CN/C(N)=N/N", "(1E)-N-methylcarbonohydrazonic diamide"),
        ("NC(=NN)NC(C)=O", "N-carbamohydrazonoylacetamide"),
        ("NC(=NNC)N", "N''-methylcarbonohydrazonic diamide"),
        ("NC(=NN(C)C)N", "N'',N''-dimethylcarbonohydrazonic diamide"),
        ("CNC(=NNC)N", "N,N''-dimethylcarbonohydrazonic diamide"),
        ("N#CNNC", "N'-methylcyanohydrazide"),
        ("NNC(=O)OC(=S)NN", "1-thiodicarbonic dihydrazide"),
        ("CNNC(=O)OC(=O)NN", "N'1-methyldicarbonic dihydrazide"),
        ("NNC(=S)SC(=S)NN", "1,2,3-trithiodicarbonic dihydrazide"),
    ],
)
def test_amides_of_cyanic_polycarbonic_and_hydrazinedicarboxylic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CN(C)CCc1cccc(c1)CN(C)C", "2-{3-[(dimethylamino)methyl]phenyl}-N,N-dimethylethan-1-amine"),
        ("C[N+](C)([O-])CCc1cccc(c1)CN(C)C", "2-{3-[(dimethylamino)methyl]phenyl}-N,N-dimethylethan-1-amine N-oxide"),
        ("CN(C)CCc1cccc(c1)C[N+](C)(C)[O-]", "1-{3-[2-(dimethylamino)ethyl]phenyl}-N,N-dimethylmethanamine N-oxide"),
    ],
)
def test_substituted_amines_on_separate_chains_one_is_the_parent_and_the_oxidized_one_wins(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C[N+](C)([O-])CCc1cccc(C[N+](C)(C)[O-])c1",
            "2-(3-{[dimethyl(oxo)-\u03bb5-azanyl]methyl}phenyl)-N,N-dimethylethan-1-amine N-oxide",
        ),
        ("C[N+](C)([O-])CC[N+](C)([S-])C", "2-[dimethyl(sulfanylidene)-\u03bb5-azanyl]-N,N-dimethylethan-1-amine N-oxide"),
        ("[O-][NH2+]CCc1cccc(C[NH2+][O-])c1", "2-{3-[(oxo-\u03bb5-azanyl)methyl]phenyl}ethan-1-amine N-oxide"),
        ("C[N+](C)([O-])c1ccccc1CC[N+](C)(C)[O-]", "2-{2-[dimethyl(oxo)-\u03bb5-azanyl]ethyl}-N,N-dimethylaniline N-oxide"),
        pytest.param(
            "[O-][N+]1(C)CCCCC1CC[N+](C)(C)[O-]",
            "2-{2-[dimethyl(oxo)-\u03bb5-azanyl]ethyl}-1-methylpiperidine 1-oxide",
            id="saturated_ring_nitrogen_oxide_outranks_the_chain_oxide",
        ),
        pytest.param(
            "[O-][n+]1ccccc1CC[N+](C)(C)[O-]",
            "2-{2-[dimethyl(oxo)-\u03bb5-azanyl]ethyl}pyridine 1-oxide",
            id="aromatic_ring_nitrogen_oxide_outranks_the_chain_oxide",
        ),
    ],
)
def test_further_amine_oxides_are_oxo_azanyl_prefixes_of_the_senior_amine_oxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=[N+]([O-])Cl", "N-chloromethanimine N-oxide"),
        ("OI(O)c1ccccc1", "phenyl-λ3-iodanediol"),
        ("ClI(Cl)C", "dichloro(methyl)-λ3-iodane"),
    ],
)
def test_imine_oxide_with_n_halogen_and_lambda_3_iodanes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NNNNC(=O)c1ccccc1", "N-(triazan-1-yl)benzamide"),
        ("CC(=O)N(C)NNN", "N-methyl-N-(triazan-1-yl)acetamide"),
        ("O=C(c1ccccc1)NNNNC(=O)c1ccccc1", "N,N'-(hydrazine-1,2-diyl)dibenzamide"),
        ("CC(=O)NNNC(C)=O", "N,N'-azanediyldiacetamide"),
        ("NNNNC(=O)O", "tetraazane-1-carboxylic acid"),
        ("OC(=O)NNNNC(=O)O", "tetraazane-1,4-dicarboxylic acid"),
        ("CNNNC(=O)O", "3-methyltriazane-1-carboxylic acid"),
        ("NNNN=C=O", "1-isocyanatotriazane"),
        ("CNNN=C=O", "1-isocyanato-2-methylhydrazine"),
        ("O=C=NN[N+](=O)[O-]", "N-isocyanatonitramide"),
    ],
)
def test_nitrogen_chain_broken_to_express_the_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "CCNC(=NCc1ccc(OC)c(F)c1)NCCCOCC1CCOCC1",
            "N-ethyl-N''-[(3-fluoro-4-methoxyphenyl)methyl]-N'-{3-[(oxan-4-yl)methoxy]propyl}guanidine",
            id="ether_in_a_substituent",
        ),
        pytest.param("CN=C(NC)NCCO", "N-(2-hydroxyethyl)-N',N''-dimethylguanidine", id="hydroxy_prefix"),
        pytest.param("CN=C(NC)NCCC(=O)OC", "methyl 3-[(N,N'-dimethylcarbamimidoyl)amino]propanoate", id="ester_outranks_guanidine"),
    ],
)
def test_guanidines_take_junior_groups_as_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "Cc1ccc(O)c(NC(=O)NC2CC(=O)N(Cc3ccccc3)C2)c1",
            "N-(1-benzyl-5-oxopyrrolidin-3-yl)-N'-(2-hydroxy-5-methylphenyl)urea",
            id="lactam_is_a_pseudoketone_junior_to_urea",
        ),
        pytest.param(
            "C=CCN(CC(=O)N1CCCN(c2ccc(-c3ccccc3)nn2)CC1)C(=O)NCC",
            "N'-ethyl-N-{[4-(6-phenylpyridazin-3-yl)-1,4-diazepane-1-carbonyl]methyl}-N-(prop-2-en-1-yl)urea",
            id="amide_of_a_ring_nitrogen_is_a_hidden_amide_junior_to_urea",
        ),
    ],
)
def test_ring_amides_do_not_outrank_a_urea(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "NNC(Cc1csc2ccccc12)CC1CCCO1",
            "[1-(1-benzothiophen-3-yl)-3-(oxolan-2-yl)propan-2-yl]hydrazine",
            id="rings_with_sulfur_and_oxygen_are_junior_to_hydrazine",
        ),
        pytest.param(
            "NNC(c1ccc(F)c(F)c1)C1CCOCC1",
            "[(3,4-difluorophenyl)(oxan-4-yl)methyl]hydrazine",
            id="two_ring_substituents_on_one_carbon",
        ),
        pytest.param("C=C(C)c1ccc(NN)cc1", "[4-(prop-1-en-2-yl)phenyl]hydrazine", id="unsaturated_substituent"),
        pytest.param(
            "CCS(=O)(=O)CCNN(C)C",
            "2-[2-(ethanesulfonyl)ethyl]-1,1-dimethylhydrazine",
            id="sulfonyl_in_a_substituent",
        ),
    ],
)
def test_hydrazine_is_the_parent_beside_rings_and_groups_junior_to_it(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NOCCN1CCCC1=O", "1-[2-(aminooxy)ethyl]pyrrolidin-2-one", id="beside_a_lactam"),
        pytest.param(
            "NOCc1cc(Br)c(C=O)c2c1OCCO2",
            "8-[(aminooxy)methyl]-6-bromo-2,3-dihydro-1,4-benzodioxine-5-carbaldehyde",
            id="beside_an_aldehyde_on_a_fused_ring",
        ),
        pytest.param(
            "Cc1ccc(CN(C)C2CCC(ON)CC2)cc1",
            "4-(aminooxy)-N-methyl-N-[(4-methylphenyl)methyl]cyclohexan-1-amine",
            id="beside_an_amine",
        ),
    ],
)
def test_aminooxy_is_a_prefix_beside_a_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("S=C(Nc1ccccc1)N1CCCC1", "N-phenylpyrrolidine-1-carbothioamide", id="pyrrolidine"),
        pytest.param("S=C(NC)N1CCOCC1", "N-methylmorpholine-4-carbothioamide", id="morpholine"),
        pytest.param("CCNC(=S)N1CCC(N)CC1", "4-amino-N-ethylpiperidine-1-carbothioamide", id="amino_on_the_ring"),
        pytest.param("S=C(N)N1CCCC1", "pyrrolidine-1-carbothioamide", id="unsubstituted_nitrogen"),
    ],
)
def test_thiourea_with_a_ring_nitrogen_is_the_carbothioamide_of_the_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)ON(C)C", "1-[(dimethylamino)oxy]ethan-1-one", id="blue_book_example"),
        pytest.param("CC(=O)ON=C(C)C", "1-{[(propan-2-ylidene)amino]oxy}ethan-1-one", id="ketoxime"),
        pytest.param("CC(=O)ON=C1CCCCC1", "1-[(cyclohexylideneamino)oxy]ethan-1-one", id="ring_ylidene"),
        pytest.param("O=C(ON=Cc1ccccc1)c1cccs1", "[(benzylideneamino)oxy](thiophen-2-yl)methanone", id="aldoxime_of_a_ring_acid"),
        pytest.param("CNC(=O)ON=C(C)C", "N-[(methylcarbamoyl)oxy]propan-2-imine", id="carbamate_of_a_ketoxime_is_a_prefix_of_the_imine"),
        pytest.param("CC(C)=NOC(=O)OC", "N-[(methoxycarbonyl)oxy]propan-2-imine", id="carbonate_of_a_ketoxime_is_a_prefix_of_the_imine"),
        pytest.param("CNC(=O)ON(C)C", "N-methyl-N-[(methylcarbamoyl)oxy]methanamine", id="carbamate_of_an_acyclic_amine_is_a_prefix_of_the_amine"),
        pytest.param("COC(=O)ON(C)C", "N-[(methoxycarbonyl)oxy]-N-methylmethanamine", id="carbonate_of_an_acyclic_amine_is_a_prefix_of_the_amine"),
        pytest.param("NC(=O)ON(C)C", "N-(carbamoyloxy)-N-methylmethanamine", id="unsubstituted_carbamate_of_an_acyclic_amine"),
        pytest.param("CNC(=O)ONC", "N-[(methylcarbamoyl)oxy]methanamine", id="carbamate_of_a_secondary_amine"),
        pytest.param("CNC(=O)ON1C(=O)CCC1=O", "2,5-dioxopyrrolidin-1-yl methylcarbamate", id="cyclic_nitrogen_keeps_the_traditional_ester"),
    ],
)
def test_acyl_derivative_of_an_oxime_or_amine_oxygen_is_a_pseudoketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)ONC(C)=O", "N-(acetyloxy)acetamide", id="acyloxy_on_an_amide_nitrogen"),
        pytest.param("CC(=O)ON(C)C(C)=O", "N-(acetyloxy)-N-methylacetamide", id="acyloxy_beside_a_second_n_substituent"),
        pytest.param("CCC(=O)ONC(=O)c1ccco1", "N-(propanoyloxy)furan-2-carboxamide", id="acyloxy_on_the_nitrogen_of_a_ring_carboxamide"),
        pytest.param("COC(=O)ONC(C)=O", "N-[(methoxycarbonyl)oxy]acetamide", id="alkoxycarbonyloxy_on_an_amide_nitrogen"),
        pytest.param("CS(=O)(=O)ONC(C)=O", "N-[(methanesulfonyl)oxy]acetamide", id="sulfonyloxy_on_an_amide_nitrogen"),
    ],
)
def test_acyloxy_group_on_an_amide_nitrogen_is_a_prefix_of_the_amide(smiles, expected):
    # P-65.6.3.2.3, P-65.6.3.4.1: an acyclic nitrogen gives no traditional ester
    assert smiles_to_iupac(smiles) == expected


# P-21.2.2, P-21.2.3.1, P-44.1.2, P-62.2, P-63.4.2, P-68.3.1: bonds between different heteroatoms
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CCN(C)NO", "1-ethyl-2-hydroxy-1-methylhydrazine", id="hydrazine_with_a_hydroxy_prefix"),
        pytest.param("CCNONC", "N-[(methylamino)oxy]ethanamine", id="amine_with_an_amino_oxy_prefix"),
        pytest.param("CN(C)SSS", "N-methyl-N-trisulfanylmethanamine", id="sulfur_chain_on_an_amine_nitrogen"),
        pytest.param("CNOS", "N-(sulfanyloxy)methanamine", id="chalcogen_pair_on_nitrogen_is_no_peroxol"),
        pytest.param("CCNP", "N-ethylphosphinous amide", id="phosphinous_amide_without_organyl_groups"),
                pytest.param("COSSS", "methoxytrisulfane", id="senior_single_chalcogen_as_an_ether_prefix"),
        pytest.param("CONOC", "N-methoxy-O-methylhydroxylamine", id="hydroxylamine_as_the_functional_parent"),
        pytest.param("CONSC", "O-methyl-N-(methylsulfanyl)hydroxylamine", id="sulfanyl_group_on_the_hydroxylamine_nitrogen"),
        pytest.param("CSNSC", "S-methyl-N-(methylsulfanyl)(thiohydroxylamine)", id="thiohydroxylamine_when_no_oxygen"),
        pytest.param("CON[SiH3]", "N-methoxysilanamine", id="alkoxy_group_on_the_silanamine_nitrogen"),
    ],
)
def test_bonds_between_different_heteroatoms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-44.1.2.2, P-31.1.4.2.4: an N=N bond closed in a ring is part of the ring parent, never a diazene
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC1CCN=N1", "3-methyl-4,5-dihydro-3H-pyrazole", id="five_membered_ring_azo"),
        pytest.param("C1CCCN=N1", "3,4,5,6-tetrahydropyridazine", id="six_membered_ring_azo"),
    ],
)
def test_ring_azo_bond_is_not_a_diazene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-66.3.2.3, P-66.3.5.3: an acylhydrazide group is cited as a hydrazido prefix, also with a substituted acyl nitrogen
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)N(CC)Nc1ccc(cc1)S(=O)(=O)O", "4-(N-ethylacetohydrazido)benzene-1-sulfonic acid", id="substituted_acyl_nitrogen"),
        pytest.param("NNC(=O)NNCC(=O)O", "(hydrazinecarbohydrazido)acetic acid", id="hydrazide_ahead_of_an_n_substituted_glycine"),
        pytest.param("O=C(O)CCCNNC(=O)OC(=O)NN", "4-{[(hydrazinecarbonyl)oxy]formohydrazido}butanoic acid", id="carbonic_residue_as_an_oxy_group"),
    ],
)
def test_hydrazido_prefix_of_a_substituted_or_carbonic_acylhydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-93.4.2.1.3, P-93.5.1.4.2.1, P-68.3.1.2.5, P-66.4.4: the configuration of a C=N bond is cited by the ylidene group, or in
# front of the name when the carbon belongs to a suffix that the ring parent does not number
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC/C=N/NC(N)=O", "2-[(1E)-propylidene]hydrazine-1-carboxamide", id="semicarbazone"),
        pytest.param("C/C(=N\\NC(N)=O)CC", "2-[(2E)-butan-2-ylidene]hydrazine-1-carboxamide", id="semicarbazone_of_a_ketone"),
        pytest.param("CC/C=N/NC(N)=S", "2-[(1E)-propylidene]hydrazine-1-carbothioamide", id="thiosemicarbazone"),
        pytest.param("c1ccccc1/C=N/Nc1ccccc1", "[(E)-benzylidene](phenyl)hydrazine", id="hydrazone_with_a_single_atom_ylidene"),
        pytest.param("N/C(=N\\O)c1ccccc1", "(Z)-N'-hydroxybenzenecarboximidamide", id="amidoxime_on_a_ring_without_a_locant"),
        pytest.param("N/C(=N/O)C1CC1", "(E)-N'-hydroxycyclopropanecarboximidamide", id="amidoxime_other_configuration"),
    ],
)
def test_configuration_of_the_c_n_bond_of_hydrazones_and_amidoximes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_azine_ylidenes_that_differ_in_configuration_are_not_ordered():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C/C=N/N=C\\C")


def test_stereodefined_ring_amine_with_a_chalcogen_in_its_n_substituent():
    assert smiles_to_iupac("CSCCN[C@H]1CCOC1") == "(3S)-N-[2-(methylsulfanyl)ethyl]oxolan-3-amine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCC=NO", "propanal oxime", id="oxime"),
        pytest.param("CCC=NOCC", "propanal O-ethyloxime", id="oxime_o_substituent"),
        pytest.param("O/N=C(\\c1ccc(Cl)cc1)c1ccccc1", "(Z)-(4-chlorophenyl)(phenyl)methanone oxime", id="configuration_of_the_c_n_bond"),
        pytest.param("CC(=NO)C(C)=NO", "butane-2,3-dione dioxime", id="every_carbonyl_converted"),
        pytest.param("CC(=NOc1ccccc1)CC(C)=NOc1ccccc1", "pentane-2,4-dione bis(O-phenyloxime)", id="multiplied_substituted_unit"),
        pytest.param("CC(=NO)C(C)=O", "butane-2,3-dione oxime", id="equivalent_carbonyls_need_no_locant"),
        pytest.param("CC(=NO)C(=O)CC", "pentane-2,3-dione 2-oxime", id="unbranched_chain_locant"),
        pytest.param("ON=C1CC(C)C(=O)CC1", "2-methylcyclohexane-1,4-dione 4-oxime", id="ring_numbered_by_its_substituent"),
        pytest.param("ON=C1C(=O)c2ccccc2C=C1", "naphthalene-1,2-dione 2-oxime", id="fused_ring_locant"),
        pytest.param("CC(=NO)C(=O)C(C)=O", "pentane-2,3,4-trione 2-oxime", id="tied_chain_ends_take_the_lower_locant"),
        pytest.param("O=CC(C)CC=NO", "2-methylbutanedial 4-oxime", id="branched_chain_aldehyde"),
        pytest.param("O=Cc1ccc(C=NO)c(C)c1", "2-methylbenzene-1,4-dicarbaldehyde 1-oxime", id="ring_aldehyde_locant_of_the_ring_atom"),
        pytest.param(
            "CC(=NO)CC(=NOC)CC(=NOCC)C",
            "heptane-2,4,6-trione 2-(O-ethyloxime) 4-(O-methyloxime) 6-oxime",
            id="different_units_tied_ends_lowest_locant_to_the_first_cited",
        ),
        pytest.param("CC/C(=N/O)/C(=N/O)C", "(2E,3Z)-pentane-2,3-dione dioxime", id="configurations_of_several_groups"),
        pytest.param("CC[C@H](C)/C=N/O", "(1E,2S)-2-methylbutanal oxime", id="configuration_beside_a_parent_stereocentre"),
        pytest.param("CC(C)=NN(C)C", "acetone dimethylhydrazone", id="hydrazone_with_retained_parent_name"),
        pytest.param(
            "CC=NNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]",
            "acetaldehyde (2,4-dinitrophenyl)hydrazone",
            id="hydrazone_compound_substituent_with_nitro_charges",
        ),
        pytest.param("CC(C)=NN=C(C)C", "acetone azine", id="symmetrical_azine"),
        pytest.param("CCC(C)=NN=C1CCCCC1", "cyclohexanone butan-2-ylidenehydrazone", id="azine_ring_component_senior"),
        pytest.param("C1CCCCC1=NN=C1CCCC1", "cyclohexanone cyclopentylidenehydrazone", id="azine_larger_ring_component_senior"),
        pytest.param("CCC=NN=CC", "propanal ethylidenehydrazone", id="azine_longer_chain_component_senior"),
        pytest.param("CC(C)=NNC(N)=O", "acetone semicarbazone", id="semicarbazone"),
        pytest.param(
            "CCCC(c1ccccc1)=NNC(=S)N(c1ccccc1)c1ccccc1",
            "1-phenylbutan-1-one 4,4-diphenylthiosemicarbazone",
            id="chalcogen_semicarbazone_with_locants",
        ),
    ],
)
def test_functional_modifier_names_of_imine_derivatives(smiles, expected):
    # P-15.2.2, P-68.3.1.1.2, P-68.3.1.2.2, P-68.3.1.2.3, P-68.3.1.2.5
    assert [a.name for a in smiles_to_alternative_names(smiles) if a.kind == "functional-modifier"] == [expected]
