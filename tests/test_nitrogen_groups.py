import pytest
from smiles_to_iupac import smiles_to_iupac
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
        pytest.param(
            "Nc1ccc(Nc2ccc(Nc3ccccc3)cc2)cc1",
            "N1-(4-aminophenyl)-N4-phenylbenzene-1,4-diamine",
            id="ring_parent_aminophenyl_branch",
        ),
        pytest.param("NCc1ccc(NC)cc1", "4-(aminomethyl)-N-methylaniline", id="ring_with_one_amine_beats_chain"),
    ],
)
def test_polyamine_n_locants_and_ring_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCNCNCC", "N,N'-methylenediethanamine", id="methylene_joins_monoamine_nitrogens"),
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
    ],
)
def test_n_substituted_aniline_and_multi_ring_chain_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC=CC", "prop-1-en-1-amine", id="enamine"),
        pytest.param("NC1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexan-1-amine", id="cyclic_amine_branch_stereocenter"),
    ],
)
def test_enamine_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_primary_amine_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("C/C=C/[C@@H](C)N") == "(2R,3E)-pent-3-en-2-amine"
    assert smiles_to_iupac("C/C=C/[C@H](C)N") == "(2S,3E)-pent-3-en-2-amine"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[N+](C)(C)[C@H](C)/C=C/C", "(2R,3E)-N,N,N-trimethylpent-3-en-2-aminium", id="ammonium_stereocenter_with_ez_double_bond_coexistence"),
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


def test_unsaturated_ring_amine_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC#C1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NCC1CCCCC1", "cyclohexylmethanamine", id="substituent_chain_amine"),
        pytest.param("NC1CCCCC1CN", "2-(aminomethyl)cyclohexan-1-amine", id="with_amine_chain_amine_tie"),
        pytest.param("NC1CCCCC1C(N)CN", "1-(2-aminocyclohexyl)ethane-1,2-diamine", id="with_amine_chain_amine_chain_wins"),
    ],
)
def test_ring_substituent_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C[NH2+][O-]", id="primary_amine_oxide_not_matched"),
    ],
)
def test_primary_amine_oxide_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    ],
)
def test_azines_and_chalcogen_semicarbazones_are_hydrazine_derivatives(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_azine_with_a_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(C)=NN=C(C)[C@H](C)CC")


def test_methyl_n_n_dimethylcarbamate():
    assert smiles_to_iupac("COC(=O)N(C)C") == "methyl dimethylcarbamate"


def test_plain_cyclic_r_group():
    assert smiles_to_iupac("O=C(N)OC1CCCCC1") == "cyclohexyl carbamate"
    assert smiles_to_iupac("NC(=O)OC1CCCC1") == "cyclopentyl carbamate"


def test_n_methyl_cyclic_r_group():
    assert smiles_to_iupac("CNC(=O)OC1CCCCC1") == "cyclohexyl methylcarbamate"


def test_two_cyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#COCCOC#N")


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
    ],
)
def test_methanediazonium_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
    assert smiles_to_iupac("c1ccccc1CC[N+]#N") == "2-phenylethane-1-diazonium"
    assert smiles_to_iupac("c1ccccc1CCC[N+]#N") == "3-phenylpropane-1-diazonium"


def test_dimethyl_benzenediazonium():
    assert smiles_to_iupac("Cc1ccc(C)cc1[N+]#N") == "2,5-dimethylbenzenediazonium"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("Cc1ccccc1CC[N+]#N", id="substituted_benzene_ring_diazonium_raises"),
        pytest.param("C=Cc1ccccc1CC[N+]#N", id="chain_diazonium_unsaturation_raises"),
    ],
)
def test_phenyl_substituted_benzene_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("N1CCCCC1=N", id="ring_fused_guanidine_still_not_supported"),
    ],
)
def test_substituted_phenyl_n_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
        pytest.param("ClC(C(=O)N(N(C)C)C)C", "2-chloro-N,N',N'-trimethylpropanehydrazide", id="n_locants_ordered_with_other_prefixes"),
        pytest.param("CC[C@@H](C)C(=O)N(C)N", "(2R)-N,2-dimethylbutanehydrazide", id="n_substituted_stereocenter_hydrazide"),
    ],
)
def test_formohydrazide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("O=CN(C)N", id="n_substituted_formohydrazide_raises"),
        pytest.param("CCCC(=O)NNC(=O)C", id="diacylhydrazide_raises"),
    ],
)
def test_n_substituted_formohydrazide_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


def test_diacyl_hydrazide_formyl():
    # Mononuclear retained name on both sides.
    assert smiles_to_iupac("O=CNNC=O") == "N'-formylformohydrazide"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(=O)NNC(=O)c1ccccc1", id="diacyl_hydrazide_aromatic_acyl_raises"),
        pytest.param("CC(C)C(=O)NNC(=O)C(C)C", id="diacyl_hydrazide_branched_acyl_raises"),
    ],
)
def test_diacyl_hydrazide_aromatic_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    ],
)
def test_ring_nitrogen_amino_is_amine_not_hydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_on_nitrogen_raises__hydrazine():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClNN")


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
    ],
)
def test_branched_symmetric_imide_and_related(smiles, expected):
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


def test_ring_not_supported__isocyanide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N+]#[C-])CC1")


def test_phenyl_isocyanide_chain():
    assert smiles_to_iupac("c1ccccc1C[N+]#[C-]") == "(isocyanomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CC[N+]#[C-]") == "(2-isocyanoethyl)benzene"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C=Cc1ccccc1C[N+]#[C-]", id="phenyl_isocyanide_unsaturation_not_supported"),
        pytest.param("C=CC[N+]#[C-]", id="unsaturated_chain_not_supported__isocyanide"),
        pytest.param("[C-]#[N+]C[N+]#[C-]", id="two_isocyanide_groups_not_supported"),
        pytest.param("CO[N+](=O)[O-].C", id="multiple_fragments_raises"),
        pytest.param("[O-][N+](=O)OCCCO[N+](=O)[O-]", id="two_nitrate_groups_raises"),
    ],
)
def test_phenyl_isocyanide_unsaturation_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("O=NOCCCON=O", id="two_nitrite_groups_raises"),
    ],
)
def test_multiple_fragments_raises_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    ],
)
def test_n_methyl_n_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_semicarbazide():
    assert smiles_to_iupac("NC(=O)NN") == "hydrazinecarboxamide"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NNC(=O)NN", id="double_amino_substituted_urea_not_supported"),
    ],
)
def test_semicarbazide_with_n_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCCCC1=NN", "cyclohexylidenehydrazine"),
        ("C1CCCCC1=NNC", "1-cyclohexylidene-2-methylhydrazine"),
        ("CC(C)=NN(C)C", "1,1-dimethyl-2-(propan-2-ylidene)hydrazine"),
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
    ],
)
def test_ring_imine_and_oxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)=Nc1ccccc1", "N-phenylpropan-2-imine", id="n_aryl_imine"),
        pytest.param("c1ccccc1C=Nc1ccccc1", "N,1-diphenylmethanimine", id="n_aryl_aldimine_of_an_aryl_aldehyde"),
        pytest.param("CC=NOCc1ccccc1", "N-(benzyloxy)ethanimine", id="oxime_ether_with_an_aromatic_chain"),
        pytest.param("CC=NOC(C)C", "N-[(propan-2-yl)oxy]ethanimine", id="branched_oxime_o_substituent"),
        pytest.param("N=CC=N", "ethane-1,2-diimine", id="polyimine"),
        pytest.param("CC(=NO)N", "N'-hydroxyethanimidamide", id="amidoxime"),
        pytest.param("NC(=NO)c1ccccc1", "N'-hydroxybenzenecarboximidamide", id="aryl_amidoxime"),
        pytest.param("N=C(N)c1ccccc1C(=N)N", "benzene-1,2-dicarboximidamide", id="ring_diamidine"),
        pytest.param("CC(=NC)N(C)C", "N,N,N'-trimethylethanimidamide", id="n_and_n_prime_substituents"),
        pytest.param("NC(=N)CC(=N)N", "propanediimidamide", id="chain_diamidine"),
    ],
)
def test_imine_and_amidine_n_substituents_and_polyfunctional_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC=NNC(N)=O", "2-ethylidenehydrazine-1-carboxamide", id="semicarbazone"),
        pytest.param("CNC(=O)NNC", "N,2-dimethylhydrazine-1-carboxamide", id="substituted_semicarbazide"),
        pytest.param("CC=NNC(=O)c1ccccc1", "N'-ethylidenebenzohydrazide", id="acylhydrazone"),
        pytest.param("O=C(NN=Cc1ccccc1)c1ccncc1", "N'-benzylidenepyridine-4-carbohydrazide", id="ring_carbohydrazide_hydrazone"),
        pytest.param("CC(=O)N(c1ccccc1)Nc1ccccc1", "N,N'-diphenylacetohydrazide", id="n_and_n_prime_substituents"),
        pytest.param("NC(=Nc1ccccn1)N", "N''-(pyridin-2-yl)guanidine", id="guanidine_with_heteroaromatic_ring"),
        pytest.param("CN(C)C(=NC)N(C)C", "N,N,N',N',N''-pentamethylguanidine", id="guanidine_with_substituents_on_every_nitrogen"),
    ],
)
def test_hydrazide_semicarbazone_and_guanidine_n_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)NC(=O)OC", "methyl N-acetylcarbamate", id="n_acyl_carbamate"),
        pytest.param("O=C(NC(=O)c1ccccc1)OC", "methyl N-benzoylcarbamate", id="n_aroyl_carbamate"),
        pytest.param("CC(=O)N(C)C(=O)OC", "methyl N-acetyl-N-methylcarbamate", id="n_acyl_n_alkyl_carbamate"),
    ],
)
def test_n_acyl_carbamates(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acylated_nitrogen_is_not_a_hetero_parent_for_a_carboxylic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC(=O)O")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CNC(=O)CCCC(=O)NC", "N1,N5-dimethylpentanediamide", id="n_locants_of_a_chain_diamide"),
        pytest.param("CNC(=O)CCC(N)=O", "N1-methylbutanediamide", id="one_substituted_group_of_a_diamide"),
        pytest.param("CCNC(=O)c1cccc(C(=O)NC)c1", "N1-ethyl-N3-methylbenzene-1,3-dicarboxamide", id="n_locants_of_a_ring_diamide"),
        pytest.param("CC(=O)NCC(N)=O", "2-acetamidoacetamide", id="second_amide_cited_as_a_prefix"),
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
        pytest.param("N=C(NN)C1CCCCC1", "cyclohexanecarboximidohydrazide", id="carboximidohydrazide_on_a_ring"),
        pytest.param("CN(C)C=NN=C(C)C", "N,N-dimethyl-N'-(propan-2-ylidene)methanehydrazonamide", id="ylidene_on_the_terminal_nitrogen"),
        pytest.param("CCN=C(c1ccccc1)N(C)Nc1ccccc1", "N''-ethyl-N-methyl-N'-phenylbenzenecarboximidohydrazide", id="three_nitrogen_locants"),
        pytest.param("NC(=NN)C(N)=NN", "ethanedihydrazonamide", id="two_hydrazonamide_groups"),
        pytest.param("NNC(=N)C(=N)NN", "ethanediimidohydrazide", id="two_imidohydrazide_groups"),
        pytest.param("NC=N", "methanimidamide", id="amidine_of_a_one_carbon_parent"),
    ],
)
def test_amidrazones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CNO", "N-hydroxymethanamine", id="n_hydroxy_amine"),
        pytest.param("CNCl", "methylhypochlorous amide", id="amide_of_hypochlorous_acid"),
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
        pytest.param("OC(=O)CN(C)C(=S)C(Cl)C", "(2-chloro-N-methylpropanethioamido)acetic acid", id="n_substituted_amido_prefix_merges_with_the_acyl_prefixes"),
        pytest.param("NC(=S)NCC(=O)O", "(carbamothioylamino)acetic acid", id="thiourea_group_stays_an_acylamino_prefix"),
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
        pytest.param("NC(=O)NC(=O)c1ccccc1", "N-carbamoylbenzamide", id="urea_beneath_a_carboxamide"),
        pytest.param("NC(=O)NS(=O)(=O)c1ccccc1", "N-carbamoylbenzenesulfonamide", id="urea_beneath_a_sulfonamide"),
        pytest.param("CC(=O)NCCNC(N)=O", "N-[2-(carbamoylamino)ethyl]acetamide", id="carbamoylamino_prefix_beneath_an_amide"),
        pytest.param("O=CNCCCNC(N)=O", "N-[3-(carbamoylamino)propyl]formamide", id="formamide_outranks_urea"),
        pytest.param("NC(=O)NC(=O)O", "1-amido-2-imidodicarbonic acid", id="carbamoyl_on_the_nitrogen_of_carbamic_acid"),
        pytest.param("CNC(=O)NC(=O)O", "(methylcarbamoyl)carbamic acid", id="substituted_carbamoyl_on_carbamic_acid"),
    ],
)
def test_ureas_with_further_substituents_and_ureas_beneath_senior_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CN(C)N=O", "dimethylnitrous amide", id="nitroso_on_amine_nitrogen"),
        pytest.param("CN(C)[N+](=O)[O-]", "dimethylnitramide", id="nitro_on_amine_nitrogen"),
        pytest.param("NN[N+](=O)[O-]", "nitric hydrazide", id="nitro_on_hydrazine_nitrogen"),
        pytest.param("NNN=O", "nitrous hydrazide", id="nitroso_on_hydrazine_nitrogen"),
        pytest.param("CN([N+](=O)[O-])[N+](=O)[O-]", "methyl(nitro)nitramide", id="second_nitro_group_is_a_prefix_of_nitramide"),
        pytest.param("CCCCN(CC)N=O", "butyl(ethyl)nitrous amide", id="nitrous_amide_prefixes_without_locants"),
        pytest.param("O=C1CCCCC1N=O", "2-nitrosocyclohexan-1-one", id="nitroso_on_a_ring_that_carries_a_ketone"),
        pytest.param("CN(N=O)c1ccc(cc1)C(=O)O", "4-[methyl(nitroso)amino]benzoic acid", id="nitrosoamino_prefix_under_an_acid"),
        pytest.param("O=NN1CCCC1", "1-nitrosopyrrolidine", id="nitroso_on_a_ring_nitrogen"),
        pytest.param("OC(=O)CN(N=O)CC(=O)O", "2,2'-(nitrosoazanediyl)diacetic acid", id="nitroso_on_a_linking_nitrogen"),
    ],
)
def test_nitroso_and_nitro_groups_on_nitrogen_and_on_rings_with_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)NOC", "N-methoxyacetamide", id="n_alkoxy_amide"),
        pytest.param("CC(=O)NOCC(C)C", "N-(2-methylpropoxy)acetamide", id="branched_alkoxy_on_amide_nitrogen"),
        pytest.param("CC(=O)NOc1ccccc1", "N-phenoxyacetamide", id="n_aryloxy_amide"),
        pytest.param("CC(=O)N(C)OC", "N-methoxy-N-methylacetamide", id="alkoxy_and_alkyl_on_one_amide_nitrogen"),
        pytest.param("O=C(NOC)CCC(=O)NOC", "N1,N4-dimethoxybutanediamide", id="n_alkoxy_groups_of_a_diamide"),
        pytest.param("CC(=O)NOC(C)C", "N-[(propan-2-yl)oxy]acetamide", id="compound_organyloxy_prefix_is_enclosed"),
        pytest.param("CS(=O)(=O)NOC", "N-methoxymethanesulfonamide", id="n_alkoxy_sulfonamide"),
    ],
)
def test_n_alkoxy_and_n_aryloxy_amides(smiles, expected):
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
        ("N#CCc1ccccc1", "phenylacetonitrile"),
        ("OCC=O", "hydroxyacetaldehyde"),
        ("C1OC1c1ccccc1", "phenyloxirane"),
        ("CC1CN1", "2-methylaziridine"),
        ("ClCC(N)=O", "2-chloroacetamide"),
    ],
)
def test_locant_is_omitted_only_when_the_parent_has_one_kind_of_substitutable_hydrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
        ("NC(=O)ONC", "(methylamino) carbamate"),
        ("NC(=O)ON", "amino carbamate"),
        ("NC(=O)ON(C)C", "(dimethylamino) carbamate"),
        ("CNC(=O)ONC", "(methylamino) N-methylcarbamate"),
    ],
)
def test_carbamic_acid_esters_of_amino_groups_cite_the_amino_group_as_the_ester_word(smiles, expected):
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
        ("N=CNc1ccc(C(=O)O)cc1", "4-methanimidamidobenzoic acid"),
        ("CC(=N)Nc1ccc(C(=O)O)cc1", "4-ethanimidamidobenzoic acid"),
        ("NS(=N)CCC(=O)O", "3-(S-aminosulfinimidoyl)propanoic acid"),
        ("NC(N)=NCCCC(=O)O", "4-[(diaminomethylidene)amino]butanoic acid"),
        ("C(N)(=N)NC=O", "N-carbamimidoylformamide"),
        ("C(N)(=N)NC(C)=O", "N-carbamimidoylacetamide"),
        ("C(N)(=N)NC(N)=O", "N-carbamimidoylurea"),
        ("C(N)(=N)NC(N)=N", "imidodicarbonimidic diamide"),
        ("C(N)(=N)NC(=N)NC(N)=N", "diimidotricarbonimidic diamide"),
        ("C(C)N=C(N(C1=CC=CC=C1)C1=CC=CC=C1)NC(N)=N", "N'1-ethyl-N1,N1-diphenylimidodicarbonimidic diamide"),
        ("N=C(NC(N)=N)NC(NC(NC(N)=N)=N)=N", "3,5,7-triimino-2,4,6,8-tetraazanonane-1,9-diimidamide"),
        ("NN=C(C)NCCC(=O)O", "3-(ethanehydrazonamido)propanoic acid"),
        ("NN=CNc1ccc(C(=O)O)cc1", "4-(methanehydrazonamido)benzoic acid"),
        ("NNC=Nc1ccc(C(=O)O)cc1", "4-[(hydrazinylmethylidene)amino]benzoic acid"),
        ("N=CNNCCC(=O)O", "3-(methanimidohydrazido)propanoic acid"),
        ("CC(=N)NNCC(=O)OC", "methyl (ethanimidohydrazido)acetate"),
        ("N=C(NNc1ccc(C(=O)O)cc1)c1ccccc1", "4-(benzenecarboximidohydrazido)benzoic acid"),
        ("NNC(=N)c1cccc(C(=O)O)c1", "3-(hydrazinecarboximidoyl)benzoic acid"),
    ],
)
def test_amidine_locants_prefixes_and_condensed_guanidines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
