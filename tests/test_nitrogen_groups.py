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
    "smiles",
    [
        pytest.param("CC(N)N(C)C", id="geminal_diamine_with_substituent_raises"),
        pytest.param("NCCN(C)CC=C", id="two_amines_with_unsaturated_n_substituent_raises"),
    ],
)
def test_phenyl_chain_amidine_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    "smiles",
    [
        pytest.param("Cc1ccccc1CN=[N+]=[N-]", id="substituted_benzene_ring_chain_not_supported"),
        pytest.param("c1ccc2ccccc2c1CN=[N+]=[N-]", id="fused_aromatic_ring_chain_not_supported"),
        pytest.param("C1CCC(N=[N+]=[N-])CC1", id="ring_not_supported"),
        pytest.param("C=CCN=[N+]=[N-]", id="unsaturated_chain_not_supported"),
    ],
)
def test_substituted_benzene_ring_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_acetone_azine_name():
    assert smiles_to_iupac("CC(C)=NN=C(C)C") == "N-(propan-2-ylideneamino)propan-2-imine"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC=NN=C(C)C", id="asymmetric_azine_raises"),
        pytest.param("C1CCCCC1=NN=C1CCCCC1", id="ring_raises"),
        pytest.param("c1ccccc1C=NN=Cc1ccccc1", id="aromatic_carbon_raises"),
        pytest.param("CC[C@@H](C)C(C)=NN=C(C)[C@H](C)CC", id="azine_specified_chain_stereocenter_raises"),
    ],
)
def test_asymmetric_azine_raises_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    "smiles",
    [
        pytest.param("C1(=[N+]=[N-])CCCCC1", id="ring_not_supported__diazo"),
        pytest.param("C=CC=[N+]=[N-]", id="unsaturated_chain_not_supported__diazo"),
        pytest.param("[N-]=[N+]=CC=[N+]=[N-]", id="two_diazo_groups_not_supported"),
    ],
)
def test_ring_not_supported_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
        pytest.param("CC[C@@H](C)C(=O)N(C)N", "(2R)-N-methyl-2-methylbutanehydrazide", id="n_substituted_stereocenter_hydrazide"),
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


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CNO", id="n_substituted_not_supported"),
        pytest.param("CNOC", id="n_o_disubstituted_not_supported"),
        pytest.param("ONC1CCCCC1", id="ring_not_supported__hydroxylamine"),
        pytest.param("C=CON", id="unsaturated_not_supported"),
    ],
)
def test_n_substituted_not_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_phenyl_o_substituent():
    assert smiles_to_iupac("NOc1ccccc1") == "O-phenylhydroxylamine"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NOCc1ccccc1", id="benzyl_o_substituent_still_not_supported"),
        pytest.param("NOc1ccccc1C", id="substituted_phenyl_o_substituent_not_supported"),
    ],
)
def test_benzyl_o_substituent_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)C(=O)NC(=O)C(C)C", "N-(2-methylpropanoyl)-2-methylpropanamide", id="branched_symmetric_imide"),
        pytest.param("ClCC(=O)NC(=O)CCl", "N-(chloroacetyl)-2-chloroethanamide", id="halogen_substituent"),
    ],
)
def test_branched_symmetric_imide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(=O)NC(=O)CC", id="unsymmetric_imide_not_supported"),
        pytest.param("C=CC(=O)NC(=O)C=C", id="unsaturated_imide_not_supported"),
    ],
)
def test_unsymmetric_imide_not_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_n_substituted_cyclic_imide_named_via_ketone_module():
    assert smiles_to_iupac("CN1C(=O)CCC1=O") == "1-methylpyrrolidine-2,5-dione"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("c1ccccc1CC(=O)NC(=O)CCc1ccccc1", id="phenyl_chain_unsymmetric_raises"),
        pytest.param("c1ccccc1C(=O)NC(=O)c1ccccc1", id="phenyl_ring_directly_on_acyl_carbon_raises"),
        pytest.param("c1ccccc1CC(=O)NC(=O)CC", id="one_sided_phenyl_chain_raises"),
        pytest.param("Cc1ccccc1CC(=O)NC(=O)Cc1ccccc1C", id="substituted_benzene_ring_raises"),
    ],
)
def test_imide_phenyl_chain_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    "smiles",
    [
        pytest.param("C1CCC(N=C=O)CC1", id="ring_not_supported__isocyanate"),
        pytest.param("C=CCN=C=O", id="unsaturated_chain_not_supported__isocyanate"),
        pytest.param("O=C=NCN=C=O", id="two_isocyanate_groups_not_supported"),
    ],
)
def test_ring_not_supported_and_related_raise_2(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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


def test_formonitrile_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C#N")


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
    assert smiles_to_iupac("C=Cc1ccccc1CC#N") == "2-(2-ethenylphenyl)acetonitrile"


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
        pytest.param("N#CCC1CCCCC1c1ccccc1", "2-(2-phenylcyclohexyl)acetonitrile", id="two_ring_aromatic_substituent_nitrile_chain_nitrile"),
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
        pytest.param("CNC(=O)N", "N-methylurea", id="n_methylurea"),
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
        pytest.param("CC(=O)NCC(N)=O", "2-(acetylamino)acetamide", id="second_amide_cited_as_a_prefix"),
    ],
)
def test_n_locants_of_polyamides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(OC(=O)c1ccccc1)C(=O)Nc1ccccc1", "1-(phenylcarbamoyl)ethyl benzoate", id="n_phenyl_carbamoyl_prefix"),
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
