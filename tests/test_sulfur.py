import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_ethyl_disulfide():
    assert smiles_to_iupac("CSSCC") == "(methyldisulfanyl)ethane"



@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CSSSSC", "dimethyltetrasulfane", id="tetrasulfide_chain"),
        pytest.param("CSS", "methanedithioperoxol", id="terminal_persulfide_methane"),
        pytest.param("CCSS", "ethanedithioperoxol", id="terminal_persulfide_ethane"),
        pytest.param("c1ccccc1SS", "benzenedithioperoxol", id="terminal_persulfide_benzene"),
        pytest.param("CCCOS", "propane-1-OS-thioperoxol", id="oxygen_sulfur_peroxol_analogue"),
        pytest.param("CSO", "methane-SO-thioperoxol", id="sulfur_oxygen_peroxol_analogue"),
        pytest.param("CC(=O)[Se]N", "1-(aminoselanyl)ethan-1-one", id="acyl_selenium_nitrogen_is_a_pseudoketone"),
        pytest.param("NSc1ccncc1", "4-(aminosulfanyl)pyridine", id="sulfur_nitrogen_link_on_a_heteroaromatic_ring"),
        pytest.param("CCOP(=O)(Cl)CSC", "ethyl [(methylsulfanyl)methyl]phosphonochloridate", id="thioether_in_the_organyl_group_of_a_phosphonic_derivative"),
    ],
)
def test_tetrasulfide_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=CSS", "ethenedithioperoxol", id="unsaturated_chain"),
        pytest.param("SSC1=CCC=CC1", "cyclohexa-1,4-diene-1-dithioperoxol", id="unsaturated_ring"),
        pytest.param("SSC1CCCCC1", "cyclohexanedithioperoxol", id="saturated_ring"),
        pytest.param("[SeH][Se]C1CCCCC1", "cyclohexanediselenoperoxol", id="selenium_analogue_on_a_ring"),
        pytest.param("OSC1CCCCC1", "cyclohexane-SO-thioperoxol", id="mixed_oxygen_sulfur_on_a_ring"),
        pytest.param("CC(=NC)CCSS", "3-(methylimino)butane-1-dithioperoxol", id="imine_is_an_imino_prefix"),
    ],
)
def test_chalcogen_peroxols_on_rings_and_unsaturated_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("SSCCSS", "ethane-1,2-bis(dithioperoxol)", id="two_groups_on_a_chain"),
        pytest.param("SSc1ccc(cc1)SS", "benzene-1,4-bis(dithioperoxol)", id="two_groups_on_a_ring"),
        pytest.param("SSCC(SS)CSS", "propane-1,2,3-tris(dithioperoxol)", id="three_groups"),
        pytest.param("OSCCSO", "ethane-1,2-bis(SO-thioperoxol)", id="mixed_pair_analogue"),
        pytest.param("SSCCOO", "2-disulfanylethane-1-peroxol", id="peroxol_outranks_the_analogue"),
        pytest.param("SOCCSS", "2-disulfanylethane-1-OS-thioperoxol", id="oxygen_first_analogue_outranks_dithioperoxol"),
        pytest.param("NCCSS", "2-aminoethane-1-dithioperoxol", id="analogue_outranks_amine"),
    ],
)
def test_several_chalcogen_peroxol_groups_take_a_bis_suffix(smiles, expected):
    # P-63.4.2.1, P-16.3.6, P-72.2.2.2.2
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)SSCC", "(2S)-2-(ethyldisulfanyl)butane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



def test_phenyl_disulfide_direct_bond():
    assert smiles_to_iupac("c1ccccc1SSC") == "(methyldisulfanyl)benzene"
    assert smiles_to_iupac("c1ccccc1SSCC") == "(ethyldisulfanyl)benzene"


@pytest.mark.parametrize(
    "smiles",
    [
    ],
)
def test_phenyl_disulfide_chain_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_phenyl_isothiocyanate_chain():
    assert smiles_to_iupac("c1ccccc1CN=C=S") == "(isothiocyanatomethyl)benzene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COS(=O)OC", "dimethyl sulfite", id="sulfite_ester_unaffected"),
        pytest.param("COS(=O)(=O)[O-].[Na+]", "sodium methyl sulfate", id="salt_of_partial_ester"),
    ],
)
def test_sulfite_ester_unaffected_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=S1Oc2ccc(C)cc2O1", "5-methyl-2H-1,3,2λ4-benzodioxathiol-2-one", id="sulfite"),
        pytest.param("O=S1(=O)Oc2ccc(C)cc2O1", "5-methyl-2H-1,3,2λ6-benzodioxathiole-2,2-dione", id="sulfate"),
    ],
)
def test_sulfur_in_a_ring_is_a_heterocycle_not_a_sulfite_or_sulfate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OS(=O)(=O)OCCOS(=O)(=O)O", "ethane-1,2-diyl bis(hydrogen sulfate)", id="sulfate"),
        pytest.param("OS(=O)OCCOS(=O)O", "ethane-1,2-diyl bis(hydrogen sulfite)", id="sulfite"),
        pytest.param("OS(=O)(=O)OCC(Cl)COS(=O)(=O)O", "2-chloropropane-1,3-diyl bis(hydrogen sulfate)", id="substituted_diyl"),
        pytest.param("OS(=O)(=O)OCC(O)COS(=O)(=O)O", "2-hydroxypropane-1,3-diyl bis(hydrogen sulfate)", id="hydroxy_diyl"),
        pytest.param("OP(O)(=O)OCCOP(O)(O)=O", "ethane-1,2-diyl bis(dihydrogen phosphate)", id="phosphate"),
        pytest.param("OP(O)(=O)OCC(OP(O)(O)=O)COP(O)(O)=O", "propane-1,2,3-triyl tris(dihydrogen phosphate)", id="triple"),
        pytest.param("OP(O)(=O)OC1CCCCC1OP(O)(O)=O", "cyclohexane-1,2-diyl bis(dihydrogen phosphate)", id="ring_diyl"),
        pytest.param("OP(O)OCCOP(O)O", "ethane-1,2-diyl bis(dihydrogen phosphite)", id="phosphite"),
        pytest.param("COP(=O)(O)OCCCOP(=O)(O)OC", "dimethyl propane-1,3-diyl bis(hydrogen phosphate)", id="partial_esters_cite_hydrogen"),
        pytest.param("CCOP(O)OCCCOP(O)OCC", "diethyl propane-1,3-diyl bis(hydrogen phosphite)", id="partial_phosphite_esters"),
        pytest.param("CS(=O)(=O)SCCSS(C)(=O)=O", "S,S\u2032-ethane-1,2-diyl dimethanesulfonothioate", id="sulfonothioate_cites_the_element_symbols"),
    ],
)
def test_multiplicative_esters_of_sulfur_and_phosphorus_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC[C@H](C)OS(=O)(=O)O", "(2S)-butan-2-yl hydrogen sulfate", id="sulfate_monoester"),
        pytest.param("CC[C@H](C)OS(=O)(=O)OC", "(2S)-butan-2-yl methyl sulfate", id="sulfate_diester"),
        pytest.param("CC[C@H](C)OS(=O)(=O)[O-]", "(2S)-butan-2-yl sulfate", id="sulfate_anion"),
        pytest.param("CC[C@H](C)OS(=O)(=O)[O-].[Na+]", "sodium (2S)-butan-2-yl sulfate", id="sulfate_salt"),
        pytest.param("C[C@H](O)COS(=O)(=O)O", "(2S)-2-hydroxypropyl hydrogen sulfate", id="hydroxy_beside_the_ester"),
        pytest.param("CC[C@H](C)O[N+](=O)[O-]", "(2S)-butan-2-yl nitrate", id="nitrate"),
        pytest.param("CC[C@H](C)ON=O", "(2S)-butan-2-yl nitrite", id="nitrite"),
    ],
)
def test_stereocentre_of_the_organyl_group_of_a_noncarbon_oxoacid_ester(smiles, expected):
    # P-67.1.3.2, P-93.5: the CIP descriptor stands in front of the organyl group of the ester
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(OS(C)(=O)=O)OS(C)(=O)=O", "ethane-1,1-diyl dimethanesulfonate", id="geminal_sulfonate_diester"),
        pytest.param(
            "CC#CC#CC(OS(=O)(=O)c1ccc(C)cc1)OS(=O)(=O)c1ccc(C)cc1",
            "hexa-2,4-diyne-1,1-diyl bis(4-methylbenzene-1-sulfonate)",
            id="acid_name_with_locants_takes_bis",
        ),
        pytest.param("CS(=O)(=O)OCC(O)COS(C)(=O)=O", "2-hydroxypropane-1,3-diyl dimethanesulfonate", id="hydroxy_prefix_on_diyl"),
        pytest.param("OCC(O)COS(=O)(=O)OS(=O)(=O)O", "2,3-dihydroxypropyl hydrogen disulfate", id="hydroxyalkyl_ester_of_disulfuric_acid"),
    ],
)
def test_sulfonate_diesters_and_hydroxyalkyl_disulfates(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OC(=O)CCOS(=O)(=N)O", "3-[(hydroxysulfonimidoyl)oxy]propanoic acid", id="imidic_o_bound"),
        pytest.param("OC(=O)CCS(=O)(=N)O", "3-(hydroxysulfonimidoyl)propanoic acid", id="imidic_c_bound"),
        pytest.param("OC(=O)CCS(=N)(=N)O", "3-(hydroxysulfonodiimidoyl)propanoic acid", id="diimidic"),
        pytest.param(
            "CCN=S(=O)(CC(=O)O)C(C)C(=O)O",
            "2-(1-carboxy-N-ethylmethanesulfonimidoyl)propanoic acid",
            id="n_substituent_of_an_imidoyl_between_two_acid_skeletons",
        ),
        pytest.param("OC(=O)CCOS(=S)(=O)O", "3-[(thiosulfo)oxy]propanoic acid", id="thio_o_bound"),
        pytest.param("OC(=O)CCS(=O)(=S)O", "3-(thiosulfo)propanoic acid", id="thio_c_bound"),
        pytest.param("OC(=O)CCS(=O)(=O)S", "3-(thiosulfo)propanoic acid", id="thio_s_acid"),
        pytest.param("OC(=O)CCOS(=S)(=S)O", "3-[(dithiosulfo)oxy]propanoic acid", id="dithio"),
        pytest.param("OC(=O)CCS(=S)O", "3-(thiosulfino)propanoic acid", id="thiosulfino"),
    ],
)
def test_modified_sulfonic_acid_groups_as_prefixes_beside_a_senior_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_methylsulfanylpropane():
    assert smiles_to_iupac("CSCCC") == "1-(methylsulfanyl)propane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)SC(C)C", "2-[(propan-2-yl)sulfanyl]propane"),
        ("CC(C)CSC(C)CC", "2-[(2-methylpropyl)sulfanyl]butane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)SCC", "(2S)-2-(ethylsulfanyl)butane"),
    ],
)
def test_stereocenter_on_parent_chain__sulfide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_substituent_branch():
    assert smiles_to_iupac("CCCCCS[C@H](C)CC") == "1-{[(2R)-butan-2-yl]sulfanyl}pentane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1SC", "(methylsulfanyl)benzene"),
        ("c1ccccc1SC(C)C", "[(propan-2-yl)sulfanyl]benzene"),
        ("c1ccccc1CSCC", "[(ethylsulfanyl)methyl]benzene"),
        ("c1ccccc1CSC(C)C", "{[(propan-2-yl)sulfanyl]methyl}benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CSc1ccccc1SC", "1,2-bis(methylsulfanyl)benzene", id="multiple_substituents"),
        pytest.param("c1ccccc1S[C@H](C)CC", "{[(2R)-butan-2-yl]sulfanyl}benzene", id="stereocenter"),
    ],
)
def test_benzene_ring_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C=C(S(=O)N)C", "prop-1-ene-2-sulfinamide", id="ene_carbon"),
        pytest.param("NS(=O)CS(=O)N", "methanedisulfinamide", id="two_sulfinamides"),
        pytest.param(
            "CN(CC)S(=O)(=O)C(C)S(=O)(=O)NC", "N1-ethyl-N1,N'1-dimethylethane-1,1-disulfonamide", id="geminal_groups_told_apart_by_primes"
        ),
        pytest.param(
            "CNS(=O)(=O)C1(S(=O)(=O)NC)CCCCC1", "N1,N'1-dimethylcyclohexane-1,1-disulfonamide", id="geminal_groups_on_a_ring_atom"
        ),
    ],
)
def test_ene_carbon_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=S(N)C1CCCCC1", "cyclohexanesulfinamide", id="cyclohexanesulfinamide"),
        pytest.param("O=S(N)C1CCCCC1Cl", "2-chlorocyclohexane-1-sulfinamide", id="2_chlorocyclohexane_1_sulfinamide"),
    ],
)
def test_cyclohexanesulfinamide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCC=C1") == "cyclohex-2-ene-1-sulfinamide"
    assert smiles_to_iupac("O=S(N)C1CC=CCC1") == "cyclohex-3-ene-1-sulfinamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=S(N)C1CCCC=C1C", "2-methylcyclohex-2-ene-1-sulfinamide", id="with_substituent"),
        pytest.param("O=S(NC)C1CCCC=C1", "N-methylcyclohex-2-ene-1-sulfinamide", id="with_n_substituent"),
    ],
)
def test_unsaturated_ring_sulfinamide_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_substituent_chain_sulfinamide():
    assert smiles_to_iupac("NS(=O)CC1CCCCC1") == "cyclohexylmethanesulfinamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NS(=O)CC1CCC(C)CC1", "(4-methylcyclohexyl)methanesulfinamide", id="ring_substituent_chain_sulfinamide_ring_with_substituent"),
        pytest.param("NS(=O)CCO", "2-hydroxyethane-1-sulfinamide", id="sulfinamide_with_alcohol"),
    ],
)
def test_ring_substituent_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1S(=O)N", "benzenesulfinamide", id="benzenesulfinamide"),
        pytest.param("Cc1ccc(cc1)S(=O)NC", "N,4-dimethylbenzene-1-sulfinamide", id="substituted_benzenesulfinamide_n_alkyl"),
    ],
)
def test_benzenesulfinamide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_sulfinamide_n_alkyl():
    assert smiles_to_iupac("c1ccccc1CCCS(=O)NC") == "N-methyl-3-phenylpropane-1-sulfinamide"


def test_phenyl_chain_sulfinamide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CS(=O)N)cc1") == "(4-chlorophenyl)methanesulfinamide"


def test_phenyl_chain_sulfinamide_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CCS(=O)N") == "2-(2-ethenylphenyl)ethane-1-sulfinamide"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CS(=O)NC", "N-methylmethanesulfinamide", id="methylmethanesulfinamide"),
        pytest.param("O=S(NC)C1CCCCC1", "N-methylcyclohexanesulfinamide", id="methylcyclohexanesulfinamide"),
    ],
)
def test_n_methylmethanesulfinamide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CS(=O)NCCCl", "N-(2-chloroethyl)methanesulfinamide", id="halogenated_n_substituent"),
        pytest.param("CS(=O)NCC=C", "N-(prop-2-en-1-yl)methanesulfinamide", id="unsaturated_n_substituted_sulfinamide"),
        pytest.param("CC[C@@H](C)S(=O)N", "(2R)-butane-2-sulfinamide", id="sulfinamide_specified_chain_stereocenter"),
    ],
)
def test_halogenated_n_substituent_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_sulfinamide_specified_sulfur_stereocenter_cited():
    assert smiles_to_iupac("CCC[S@](=O)N") == "(S)-propane-1-sulfinamide"
    assert smiles_to_iupac("CCC[S@@](=O)N") == "(R)-propane-1-sulfinamide"
    assert smiles_to_iupac("CC[C@H](C)[S@](=O)N") == "(2S,S)-butane-2-sulfinamide"
    assert smiles_to_iupac("CC[C@H](C)[S@@](=O)c1ccccc1") == "[(R)-(2S)-butane-2-sulfinyl]benzene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(S(=O)O)C", "propane-2-sulfinic acid"),
    ],
)
def test_saturated_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_2_chlorocyclohexane_1_sulfinic_acid():
    assert smiles_to_iupac("OS(=O)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfinic acid"


def test_unsaturated_ring_sulfinic_acid():
    assert smiles_to_iupac("OS(=O)C1CCCC=C1") == "cyclohex-2-ene-1-sulfinic acid"
    assert smiles_to_iupac("OS(=O)C1CC=CCC1") == "cyclohex-3-ene-1-sulfinic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OS(=O)CC1CCCCC1", "cyclohexylmethanesulfinic acid", id="ring_substituent_chain_sulfinic_acid"),
        pytest.param("c1ccccc1S(=O)O", "benzenesulfinic acid", id="benzenesulfinic_acid"),
        pytest.param("Cc1ccccc1S(=O)O", "2-methylbenzene-1-sulfinic acid", id="substituted_benzenesulfinic_acid"),
        pytest.param("Clc1ccc(CCCS(=O)O)cc1", "3-(4-chlorophenyl)propane-1-sulfinic acid", id="phenyl_chain_sulfinic_acid_ring_halogen"),
        pytest.param("CCc1ccc(cc1)CS(=O)O", "(4-ethylphenyl)methanesulfinic acid", id="phenyl_chain_sulfinic_acid_ring_ethyl"),
    ],
)
def test_sulfinic_acid_chains_with_ring_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COS(=O)(=O)OC", "dimethyl sulfate", id="sulfate_ester_unaffected"),
        pytest.param("COS(=O)[O-].[Na+]", "sodium methyl sulfite", id="salt_of_partial_ester__sulfite"),
    ],
)
def test_sulfate_ester_unaffected_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=S(=O)(N)C1CCCCC1", "cyclohexanesulfonamide", id="cyclohexanesulfonamide"),
        pytest.param("O=S(=O)(NC)C1CCCCC1", "N-methylcyclohexanesulfonamide", id="n_methylcyclohexanesulfonamide"),
    ],
)
def test_cyclohexanesulfonamide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC=C1") == "cyclohex-2-ene-1-sulfonamide"
    assert smiles_to_iupac("O=S(=O)(N)C1CC=CCC1") == "cyclohex-3-ene-1-sulfonamide"


def test_unsaturated_ring_sulfonamide_with_substituent():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-sulfonamide"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NS(=O)(=O)CC1CCCCC1", "cyclohexylmethanesulfonamide", id="ring_substituent_chain_sulfonamide"),
        pytest.param("NS(=O)(=O)CC1CCC(C)CC1", "(4-methylcyclohexyl)methanesulfonamide", id="ring_substituent_chain_sulfonamide_ring_with_substituent"),
        pytest.param("c1ccccc1S(=O)(=O)N", "benzenesulfonamide", id="benzenesulfonamide"),
        pytest.param("Cc1ccc(cc1)S(=O)(=O)NC", "N,4-dimethylbenzene-1-sulfonamide", id="substituted_benzenesulfonamide_n_alkyl"),
        pytest.param("c1ccccc1CCCS(=O)(=O)NC", "N-methyl-3-phenylpropane-1-sulfonamide", id="phenyl_chain_sulfonamide_n_alkyl"),
        pytest.param("Clc1ccc(CCCS(=O)(=O)N)cc1", "3-(4-chlorophenyl)propane-1-sulfonamide", id="phenyl_chain_sulfonamide_ring_halogen"),
        pytest.param("CCc1ccc(cc1)CS(N)(=O)=O", "(4-ethylphenyl)methanesulfonamide", id="phenyl_chain_sulfonamide_ring_ethyl"),
        pytest.param("C=Cc1ccccc1CCS(=O)(=O)N", "2-(2-ethenylphenyl)ethane-1-sulfonamide", id="phenyl_chain_sulfonamide_unsaturation"),
        pytest.param("CS(=O)(=O)N(CC)C(C)(C)C", "N-tert-butyl-N-ethylmethanesulfonamide", id="two_different_n_substituents_alphabetized_ignoring_italic_prefix"),
        pytest.param("CS(=O)(=O)NCCCl", "N-(2-chloroethyl)methanesulfonamide", id="halogenated_n_substituent"),
        pytest.param("CS(=O)(=O)NCC=C", "N-(prop-2-en-1-yl)methanesulfonamide", id="unsaturated_n_substituted_sulfonamide"),
        pytest.param("C[C@@H](Cl)S(=O)(=O)N", "(1R)-1-chloroethane-1-sulfonamide", id="acyclic_sulfonamide_stereocenter_with_coexisting_substituent"),
        pytest.param("N[S](=O)(=O)[C@H]1CCCC[C@@H]1Cl", "(1S,2S)-2-chlorocyclohexane-1-sulfonamide", id="cyclic_sulfonamide_stereocenter"),
        pytest.param("NS(=O)(=O)C1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexane-1-sulfonamide", id="cyclic_sulfonamide_branch_stereocenter"),
    ],
)
def test_ring_substituent_chain_and_related_2(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCS(=O)(=O)[O-]", "chloromethanesulfonate"),
    ],
)
def test_sulfonate_two_carbon_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)S(=O)(=O)[O-]", "(2S)-butane-2-sulfonate"),
    ],
)
def test_sulfonate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-]S(=O)(=O)CCCS(=O)(=O)[O-]", "propane-1,3-disulfonate", id="multiple_sulfonate_groups_is_named"),
        pytest.param("[O-]S(=O)(=O)C1CCCCC1", "cyclohexanesulfonate", id="sulfonate_on_ring_is_named"),
        pytest.param("C=CS(=O)(=O)[O-]", "ethenesulfonate", id="sulfonate_carbon_in_double_bond_is_named"),
        pytest.param("NCCS(=O)(=O)[O-]", "2-aminoethane-1-sulfonate", id="sulfonate_other_heteroatom_is_named"),
    ],
)
def test_multiple_sulfonate_groups_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCS(=O)(=O)CC", "(ethanesulfonyl)ethane"),
        ("CS(=O)(=O)CCCC", "1-(methanesulfonyl)butane"),
    ],
)
def test_sulfone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=S1(=O)CCCCC1", "1λ6-thiane-1,1-dione", id="ring_sulfone_is_a_lambda6_heterone"),
        pytest.param("OS(=O)(=O)C1CCCC=C1C", "2-methylcyclohex-2-ene-1-sulfonic acid", id="unsaturated_ring_sulfonic_acid_with_substituent"),
    ],
)
def test_ring_sulfone_is_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OS(=O)(=O)CC1CCCCC1", "cyclohexylmethanesulfonic acid", id="ring_substituent_chain_sulfonic_acid"),
        pytest.param("OS(=O)(=O)CC1CCC(C)CC1", "(4-methylcyclohexyl)methanesulfonic acid", id="ring_substituent_chain_sulfonic_acid_ring_with_substituent"),
        pytest.param("OS(=O)(=O)CC1CCCC=C1", "(cyclohex-2-en-1-yl)methanesulfonic acid", id="ring_substituent_chain_sulfonic_acid_unsaturated_ring"),
        pytest.param("c1ccccc1S(=O)(=O)O", "benzenesulfonic acid", id="benzenesulfonic_acid"),
        pytest.param("C=Cc1ccccc1CCS(=O)(=O)O", "2-(2-ethenylphenyl)ethane-1-sulfonic acid", id="phenyl_chain_sulfonic_acid_unsaturation"),
        pytest.param("OS(=O)(=O)c1cccnc1", "pyridine-3-sulfonic acid", id="heteroaromatic_direct_attachment_sulfonic_acid"),
        pytest.param("C[C@@H](Cl)S(=O)(=O)O", "(1R)-1-chloroethane-1-sulfonic acid", id="acyclic_sulfonic_acid_stereocenter_with_coexisting_substituent"),
        pytest.param("O=S(=O)(O)[C@H]1CCCC[C@@H]1Cl", "(1S,2S)-2-chlorocyclohexane-1-sulfonic acid", id="cyclic_sulfonic_acid_stereocenter"),
        pytest.param("OS(=O)(=O)C1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexane-1-sulfonic acid", id="cyclic_sulfonic_acid_branch_stereocenter"),
        pytest.param("[SH3+]", "sulfanium", id="sulfanium"),
    ],
)
def test_ring_substituent_chain_and_related_3(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)[SH2+]", "(propan-2-yl)sulfanium", id="branched_substituent"),
        pytest.param("[SH2+]C1CCCCC1", "cyclohexylsulfanium", id="ring_substituent"),
        pytest.param("Cc1ccccc1[SH2+]", "(2-methylphenyl)sulfanium", id="substituted_phenyl"),
        pytest.param("C[S+](C)CC=C", "dimethyl(prop-2-en-1-yl)sulfanium", id="unsaturated_substituent"),
        pytest.param("CCOCC[S+](C)CCOCC", "bis(2-ethoxyethyl)(methyl)sulfanium", id="multiplied_compound_first"),
    ],
)
def test_sulfonium_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)C[S+](C)C", "dimethyl(2-oxopropyl)sulfanium", id="ketone_group_in_a_substituent"),
        pytest.param("OCC[S+](C)C", "(2-hydroxyethyl)di(methyl)sulfanium", id="hydroxy_group_in_a_substituent"),
    ],
)
def test_sulfonium_substituent_with_characteristic_group_is_a_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CS[S+](SC)C", "1,2,3-trimethyltrisulfan-2-ium", id="cationic_centre_can_sit_on_another_sulfur"),
        pytest.param("C[S+](C)SSC", "1,1,3-trimethyltrisulfan-1-ium", id="cationic_centre_at_the_end_of_the_chain"),
        pytest.param("C[S+](C)SC", "trimethyldisulfanium", id="symmetric_chain_with_every_position_filled"),
        pytest.param("F[S+](C)(F)(F)F", "tetrafluoro(methyl)-λ4-sulfanium", id="hypervalent_sulfonium_centre"),
        pytest.param("C[S+](C)(C)=O", "trimethyl(oxo)-λ4-sulfanium", id="double_bonded_group_on_a_hypervalent_centre"),
    ],
)
def test_sulfanium_chain_locants_and_lambda_centre(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        pytest.param("C[S@+](CC)CCC", "(R)-ethyl(methyl)(propyl)sulfanium", id="sulfanium_centre"),
        pytest.param("C[S@@+](CC)c1ccccc1.[Cl-]", "(S)-ethyl(methyl)(phenyl)sulfanium chloride", id="salt"),
        pytest.param("C[Se@+](CC)c1ccccc1.[Cl-]", "(R)-ethyl(methyl)(phenyl)selanium chloride", id="selanium_centre"),
        pytest.param("C[S@+](CC)CC(=O)[O-]", "[(R)-ethyl(methyl)sulfaniumyl]acetate", id="prefix_of_an_anionic_parent"),
        pytest.param("C[N@+](CC)(CCC)CC(=O)[O-]", "[(S)-ethyl(methyl)(propyl)azaniumyl]acetate", id="azaniumyl_prefix"),
        pytest.param("C[S@+]1CCC[C@H]1C", "(1S,2R)-1,2-dimethylthiolan-1-ium", id="ring_centre"),
        pytest.param("[S@+](C)(CC)CCC", "(S)-ethyl(methyl)(propyl)sulfanium", id="centre_without_preceding_atom"),
    ],
)
def test_stereogenic_sulfanium_and_selanium_centres_cite_their_descriptor(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_and_phenyl_substituents():
    assert smiles_to_iupac("C[S+](C)c1ccccc1") == "dimethyl(phenyl)sulfanium"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCCS(=O)CC", "1-(ethanesulfinyl)butane"),
        ("CCS(=O)CC", "(ethanesulfinyl)ethane"),
    ],
)
def test_sulfoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=S1(=O)CCCC1", "1λ6-thiolane-1,1-dione"),
        ("O=S1(=O)CCCCN1", "1λ6,2-thiazinane-1,1-dione"),
        ("O=S1(=O)CCCN1", "1λ6,2-thiazolidine-1,1-dione"),
        ("O=S1(=O)CCCN1C", "2-methyl-1λ6,2-thiazolidine-1,1-dione"),
        ("CC1CCS(=O)(=O)N1", "3-methyl-1λ6,2-thiazolidine-1,1-dione"),
        ("O=S1CCCCN1", "1λ4,2-thiazinan-1-one"),
        ("O=S1(=O)OCCC1", "1,2λ6-oxathiolane-2,2-dione"),
        ("CC1CCCOS1(=O)=O", "3-methyl-1,2λ6-oxathiane-2,2-dione"),
        ("S=S1OCCC1", "1,2λ4-oxathiolane-2-thione"),
        ("O=S1CCOC1", "1,3λ4-oxathiolan-3-one"),
        ("O=C1CCS(=O)(=O)C1", "1λ6-thiolane-1,1,3-trione"),
        ("O=S1C=CC=C1", "1H-1λ4-thiophen-1-one"),
        ("O=S1(=O)C=CC=C1C", "2-methyl-1H-1λ6-thiophene-1,1-dione"),
        ("O=S1N=CC=C1", "1H-1λ4,2-thiazol-1-one"),
        ("O=[Se]1C=CC=CC=C1", "1H-1λ4-selenepin-1-one"),
    ],
)
def test_ring_chalcogen_with_doubly_bonded_chalcogens_is_a_lambda_heterone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)C1CCS(=O)(=O)C1", "1,1-dioxo-1λ6-thiolane-3-carboxylic acid", id="sulfone_ring_beside_an_acid"),
        pytest.param("N#CC1CCS(=O)CC1", "1-oxo-1λ4-thiane-4-carbonitrile", id="sulfoxide_ring_beside_a_nitrile"),
        pytest.param("OC(=O)CC1CCS(=O)(=O)C1", "(1,1-dioxo-1λ6-thiolan-3-yl)acetic acid", id="sulfone_ring_as_a_substituent_group"),
        pytest.param("OC(=O)C1CCS(=O)(=O)N1", "1,1-dioxo-1λ6,2-thiazolidine-3-carboxylic acid", id="hetero_ring_with_a_second_heteroatom"),
        pytest.param("OC(=O)C1CC[Se](=O)(=O)C1", "1,1-dioxo-1λ6-selenolane-3-carboxylic acid", id="selenium_ring"),
        pytest.param("OC(=O)C1CCS(=N)(=O)C1", "1-imino-1-oxo-1λ6-thiolane-3-carboxylic acid", id="imino_and_oxo_on_one_atom"),
        pytest.param("OC(=O)C1CCP(=O)(C)CC1", "1-methyl-1-oxo-1λ5-phosphinane-4-carboxylic acid", id="phosphorus_ring"),
        pytest.param("O=S1(=O)C=CC=C1C(=O)O", "1,1-dioxo-1H-1λ6-thiophene-2-carboxylic acid", id="mancude_ring_cites_the_indicated_hydrogen"),
        pytest.param("OC(=O)c1ccc2S(=O)(=O)C=Cc2c1", "1,1-dioxo-1H-1λ6-benzothiophene-5-carboxylic acid", id="fused_mancude_ring"),
        pytest.param("OC(=O)c1ccc2c(c1)S(=O)(=O)CC2", "1,1-dioxo-2,3-dihydro-1H-1λ6-benzothiophene-6-carboxylic acid", id="fused_ring_with_hydro_prefixes"),
        pytest.param("OC(=O)c1ccc2c(c1)S(=O)(=O)c1ccccc1-2", "5,5-dioxo-5H-5λ6-dibenzo[b,d]thiophene-3-carboxylic acid", id="fused_ring_of_three_rings"),
        pytest.param(
            "OC(=O)c1ccc2c(c1)S(=O)(=O)c1ccccc1S2(=O)=O",
            "5,5,10,10-tetraoxo-5λ6,10λ6-thianthrene-2-carboxylic acid",
            id="two_sulfone_atoms_need_no_indicated_hydrogen",
        ),
    ],
)
def test_ring_heteroatom_oxo_and_imino_prefixes_cite_the_lambda_number(smiles, expected):
    # P-63.6, P-64.4.2, P-14.1.3
    assert smiles_to_iupac(smiles) == expected


# P-68.4.3.6, P-21.1.2, P-22.2.2: the imino groups of a ring sulfoximide are prefixes beside the suffix 'one'
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=S1(=N)CCCC1", "1-imino-1λ6-thiolan-1-one", id="nh_on_a_five_membered_ring"),
        pytest.param("O=S1(=N)CCC1", "1-imino-1λ6-thietan-1-one", id="nh_on_a_four_membered_ring"),
        pytest.param("O=S1(=N)CCCCC1", "1-imino-1λ6-thian-1-one", id="nh_on_a_six_membered_ring"),
        pytest.param("O=S1(=NC)CCC1", "1-(methylimino)-1λ6-thietan-1-one", id="alkyl_on_the_imino_nitrogen"),
        pytest.param("O=S1(=NC(C)C)CCCC1", "1-[(propan-2-yl)imino]-1λ6-thiolan-1-one", id="branched_alkyl_on_the_imino_nitrogen"),
        pytest.param("O=S1(=Nc2ccccc2)CCCC1", "1-(phenylimino)-1λ6-thiolan-1-one", id="aryl_on_the_imino_nitrogen"),
        pytest.param("O=S1(=NC)CC(C)CC1", "3-methyl-1-(methylimino)-1λ6-thiolan-1-one", id="substituent_on_the_ring_carbon"),
        pytest.param("O=S1(=N)CCOCC1", "4-imino-1,4λ6-oxathian-4-one", id="second_heteroatom_in_the_ring"),
        pytest.param("S=S1(=N)CCCC1", "1-imino-1λ6-thiolane-1-thione", id="thione_beside_the_imino_prefix"),
    ],
)
def test_imino_prefixes_of_a_ring_sulfoximide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-41, P-44.1.2.2: a ring with oxygen or sulfur in the imino substituent is senior to the lambda-sulfane chain or ring
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CS(=O)(=NC1COC1)C", "3-{[dimethyl(oxo)-\u03bb6-sulfanylidene]amino}oxetane", id="oxygen_ring_on_an_acyclic_sulfoximide"),
        pytest.param("CS(=O)(=NC1CCSCC1)C", "4-{[dimethyl(oxo)-\u03bb6-sulfanylidene]amino}thiane", id="sulfur_ring_on_an_acyclic_sulfoximide"),
        pytest.param("C1CCS(=O)(=NC2COC2)CC1", "3-[(1-oxo-1\u03bb6-thian-1-ylidene)amino]oxetane", id="oxygen_ring_beside_a_ring_sulfoximide"),
    ],
)
def test_ring_in_the_imino_substituent_is_the_parent_of_a_sulfoximide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-64.4.2, P-21.2.4.1, P-22.2.2, P-68.3.1: ring atoms of nonstandard bonding number and heteroatom groups on ring nitrogen
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("N=S1CCCC1", "1λ4-thiolan-1-imine", id="imine_on_a_ring_sulfur"),
        pytest.param("N=S1C=CSC1", "1λ4,3-dithiol-1(2H)-imine", id="imine_on_a_ring_sulfur_with_added_hydrogen"),
        pytest.param("N=P1CCCC1", "1λ5-phospholan-1-imine", id="imine_on_a_ring_phosphorus"),
        pytest.param("CN(O)N1CC1", "N-hydroxy-N-methylaziridin-1-amine", id="heteroatom_group_on_the_amine_of_a_ring_nitrogen"),
        pytest.param("SON1CCC1", "azetidine-1-OS-thioperoxol", id="chalcogen_pair_on_a_ring_nitrogen"),
        pytest.param("OON1CC1", "aziridine-1-peroxol", id="peroxol_on_a_ring_nitrogen"),
    ],
)
def test_ring_heteroatom_imines_and_heteroatom_groups_on_ring_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,warns",
    [("O=S1C=CC=C1", False), ("O=S1CCCC1", False)],
)
def test_non_pin_retained_name_cases_warn(smiles, warns):
    import warnings

    from smiles_to_iupac import NonPreferredNameWarning

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        smiles_to_iupac(smiles)
    assert any(issubclass(w.category, NonPreferredNameWarning) for w in caught) == warns


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=S1CCCCC1", "1λ4-thian-1-one", id="ring_resolves_via_hetero_ring_oxide"),
        pytest.param("[S-]C(=O)CCC(=O)[S-]", "butanebis(thioate)", id="multiple_thioate_groups_is_named"),
        pytest.param("[S-]C(=O)C1CCCCC1", "cyclohexanecarbothioate", id="thioate_on_ring_is_named"),
        pytest.param("CC(C)C(=O)[S-]", "2-methylpropanethioate", id="thioate_branched_r_group"),
    ],
)
def test_ring_resolves_via_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(=O)[S-]", "(2S)-2-methylbutanethioate"),
    ],
)
def test_thioate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1C(=O)[S-]", "benzenecarbothioate", id="directly_attached_thioate_is_named"),
        pytest.param("Clc1ccc(CC(=O)[S-])cc1", "(4-chlorophenyl)ethanethioate", id="chain_thioate_ring_halogen"),
        pytest.param("C=Cc1ccccc1CC(=O)[S-]", "(2-ethenylphenyl)ethanethioate", id="chain_thioate_unsaturation_is_named"),
    ],
)
def test_phenyl_directly_attached_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("N#CSCCSC#N", "ethane-1,2-diyl dithiocyanate", id="chain"),
        pytest.param("N#C[Se]Cc1ccc(C[Se]C#N)cc1", "1,4-phenylenebis(methylene) diselenocyanate", id="polyol_group_named_as_in_an_acetate"),
        pytest.param("N#COC[C@H](C)OC#N", "(2S)-propane-1,2-diyl dicyanate", id="stereocentre_of_the_organyl_group"),
    ],
)
def test_cyanic_esters_of_a_polyol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_thiocyanate():
    assert smiles_to_iupac("c1ccccc1SC#N") == "phenyl thiocyanate"  # CID 21357


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=S)O", "ethanethioic O-acid"),
        ("O=CS", "methanethioic S-acid"),
        ("O=C(S)C(C)Oc1ccccc1", "2-phenoxypropanethioic S-acid"),
        ("O=C(S)C(C)Sc1ccccc1", "2-(phenylsulfanyl)propanethioic S-acid"),
    ],
)
def test_thioic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C1CCCCC1S", "2-sulfanylcyclohexan-1-one", id="ring"),
        pytest.param("CC(C)Cc1ccc(cc1)C(C)C(=O)S", "2-[4-(2-methylpropyl)phenyl]propanethioic S-acid", id="phenyl_substituent_thioic_acid_branched_chain"),
        pytest.param("SC1CCCC=C1C", "2-methylcyclohex-2-ene-1-thiol", id="unsaturated_ring_thiol_with_substituent"),
    ],
)
def test_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("SCC1CCCC=C1", "(cyclohex-2-en-1-yl)methanethiol", id="ring_substituent_chain_thiol_unsaturated_ring"),
        pytest.param("SC1CCCCC1CCS", "2-(2-sulfanylethyl)cyclohexane-1-thiol", id="ring_with_thiol_chain_thiol_tie"),
        pytest.param("SCCO", "2-sulfanylethan-1-ol", id="thiol_with_alcohol"),
        pytest.param("c1ccccc1S", "benzenethiol", id="benzenethiol"),
        pytest.param("Sc1ccccc1S", "benzene-1,2-dithiol", id="two_direct_ring_thiols"),
        pytest.param("C=Cc1ccccc1CCS", "2-(2-ethenylphenyl)ethane-1-thiol", id="phenyl_chain_thiol_unsaturation"),
    ],
)
def test_ring_substituent_chain_and_related_4(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)S", "(2S)-butane-2-thiol"),
    ],
)
def test_acyclic_thiol_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("SC1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexane-1-thiol", id="cyclic_thiol_branch_stereocenter"),
        pytest.param("S[C@H]1CCCCC1Cl", "(1S)-2-chlorocyclohexane-1-thiol", id="thiol_partially_specified_stereocenters_cites_the_specified_elements"),
    ],
)
def test_cyclic_thiol_branch_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("SCCc1cccnc1", "2-(pyridin-3-yl)ethane-1-thiol"),
        ("SC1CCCCC1c1ccc[nH]1", "2-(1H-pyrrol-2-yl)cyclohexane-1-thiol"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("SC1CCCC(c2ccccc2)C1c3ccccc3", "2¹,2²,2³,2⁴,2⁵,2⁶-hexahydro[1¹,2¹:2²,3¹-terphenyl]-2³-thiol", id="three_rings"),
        pytest.param("SCC1CCCCC1c1ccccc1", "(2-phenylcyclohexyl)methanethiol", id="chain_thiol"),
    ],
)
def test_two_ring_aromatic_substituent_thiol_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCC(=S)C", "1-chloropropane-2-thione"),
    ],
)
def test_thione_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=S)C1CCCCC1", "1-cyclohexylethane-1-thione", id="substituent_chain_thione"),
        pytest.param("S=C1CCCCC1CC(=S)C", "2-(2-sulfanylidenepropyl)cyclohexane-1-thione", id="with_thione_chain_thione_tie"),
    ],
)
def test_ring_substituent_chain_and_related_5(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Cc1ccc(cc1)C(=S)C", "1-(4-methylphenyl)ethane-1-thione", id="aromatic_thione_is_named"),
        pytest.param("CC[C@@H](C)C(C)=S", "(3R)-3-methylpentane-2-thione", id="acyclic_thione_stereocenter"),
        pytest.param("S=C1CCCC[C@H]1C", "(2R)-2-methylcyclohexane-1-thione", id="cyclic_thione_stereocenter"),
        pytest.param("S=C1CCC(CC1)[C@@H](C)CC", "4-[(2S)-butan-2-yl]cyclohexane-1-thione", id="cyclic_thione_branch_stereocenter"),
        pytest.param("c1ccccc1CC(=S)C", "1-phenylpropane-2-thione", id="phenyl_chain_thione"),
        pytest.param("c1ccccc1C(=S)C", "1-phenylethane-1-thione", id="phenyl_directly_attached_thione_is_named"),
        pytest.param("C=Cc1ccccc1CC(=S)C", "1-(2-ethenylphenyl)propane-2-thione", id="phenyl_chain_thione_unsaturation_is_named"),
        pytest.param("NC(=S)N", "thiourea", id="thiourea"),
        pytest.param("CNC(=S)N", "methylthiourea", id="n_methylthiourea"),
        pytest.param("CN(C)C(=S)N", "N,N-dimethylthiourea", id="n_n_dimethylthiourea_same_nitrogen"),
        pytest.param("CCN(C)C(=S)N", "N-ethyl-N-methylthiourea", id="n_ethyl_n_methylthiourea_same_nitrogen"),
    ],
)
def test_aromatic_thione_is_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CNC(=S)Nc1ccccc1", "N-methyl-N'-phenylthiourea", id="methyl_n_prime_phenylthiourea_different_nitrogens"),
        pytest.param("c1ccc(NC(=S)Nc2ccccc2)cc1", "N,N'-diphenylthiourea", id="n_prime_diphenylthiourea"),
        pytest.param("CCN(C)C(=S)NC", "N-ethyl-N,N'-dimethylthiourea", id="more_substituents_take_the_unprimed_locant"),
        pytest.param("S=C(Nc1ccccc1)Nc1ccccn1", "N-phenyl-N'-(pyridin-2-yl)thiourea", id="heteroaromatic_ring_n_substituent"),
    ],
)
def test_n_methyl_n_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1=S", "2-hydroxycyclohexane-1-thione"),
        ("OC1CCCCC1C(=S)C1CCCCC1O", "bis(2-hydroxycyclohexyl)methanethione"),
        ("CC(=S)CC(C)=O", "4-sulfanylidenepentan-2-one"),
    ],
)
def test_thione_outranks_hydroxy_and_yields_to_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CS(=O)(=O)NN", "methanesulfonohydrazide", id="sulfonohydrazide"),
        pytest.param("CS(=O)NN", "methanesulfinohydrazide", id="sulfinohydrazide"),
        pytest.param("CS(=O)(=O)NNC", "N'-methylmethanesulfonohydrazide", id="terminal_nitrogen_substituent"),
        pytest.param("CS(=O)(=O)N(C)NC", "N,N'-dimethylmethanesulfonohydrazide", id="both_nitrogens_substituted"),
        pytest.param("NNS(=O)(=O)CCS(=O)(=O)NN", "ethane-1,2-disulfonohydrazide", id="two_sulfonohydrazide_groups"),
        pytest.param("CNNS(=O)(=O)CS(=O)(=O)NNC", "N',N'''-dimethylmethanedisulfonohydrazide", id="geminal_sulfonohydrazide_groups"),
        pytest.param("c1ccccc1S(=O)(=O)NN", "benzenesulfonohydrazide", id="sulfonohydrazide_on_a_ring"),
        pytest.param("NNS(=O)(=O)c1ccccc1C(=O)O", "2-(hydrazinesulfonyl)benzoic acid", id="senior_acid_cites_the_prefix"),
        pytest.param("NNS(=O)(=S)C", "methanesulfonothiohydrazide", id="replaced_oxygen_in_a_sulfonohydrazide"),
        pytest.param("C[Se](=O)(=S)NNC", "N'-methylmethaneselenonothiohydrazide", id="replaced_oxygen_on_a_selenium_centre"),
    ],
)
def test_sulfonohydrazides_and_sulfinohydrazides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CS(=N)(=O)N", "methanesulfonimidamide", id="sulfonimidamide"),
        pytest.param("CS(=N)(=N)N", "methanesulfonodiimidamide", id="sulfonodiimidamide"),
        pytest.param("CS(=N)N", "methanesulfinimidamide", id="sulfinimidamide"),
        pytest.param("CS(=O)(=NC)NC", "N,N'-dimethylmethanesulfonimidamide", id="amino_and_imido_substituents"),
        pytest.param("CS(=N)(=NC)N(C)C", "N,N,N'-trimethylmethanesulfonodiimidamide", id="substituted_nitrogens_take_the_lower_locants"),
        pytest.param("c1ccccc1S(=N)(=O)N", "benzenesulfonimidamide", id="sulfonimidamide_on_a_ring"),
        pytest.param("NS(=O)(=N)CCS(=O)(=N)N", "ethane-1,2-disulfonimidamide", id="two_sulfonimidamide_groups"),
        pytest.param("CNS(=N)(=O)CCS(=N)(=O)NC", "N1,N2-dimethylethane-1,2-disulfonimidamide", id="n_substituted_sulfonimidamide_groups"),
        pytest.param("CS(=S)(=N)N", "methanesulfonimidothioamide", id="replaced_oxygen_in_a_sulfonimidamide"),
        pytest.param("CS(=N)(=O)N1CCOCC1", "4-(methanesulfonimidoyl)morpholine", id="sulfonimidoyl_on_a_ring_nitrogen"),
        pytest.param("CS(=NC)(=O)n1cccc1", "1-(N-methylmethanesulfonimidoyl)-1H-pyrrole", id="n_substituted_sulfonimidoyl_on_an_aromatic_nitrogen"),
        pytest.param("c1ccccc1S(=S)(=N)NC", "N-methylbenzenesulfonimidothioamide", id="replaced_oxygen_on_a_ring"),
        pytest.param("CS(=O)(=N)NN", "methanesulfonimidohydrazide", id="sulfonimidohydrazide"),
        pytest.param("CS(=[Se])(=NC)NNC", "N',N''-dimethylmethanesulfonimidoselenohydrazide", id="primes_continue_over_the_hydrazide_and_imido_nitrogens"),
    ],
)
def test_sulfonimidamides_and_sulfinimidamides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NS(=O)(=O)O", "sulfamic acid", id="sulfamic_acid"),
        pytest.param("CCNS(=O)(=O)O", "N-ethylsulfamic acid", id="substituted_sulfamic_acid"),
        pytest.param("NS(=O)(=O)[O-]", "sulfamate", id="sulfamate"),
        pytest.param("[Na+].NS(=O)(=O)[O-]", "sodium sulfamate", id="sulfamate_in_a_salt"),
        pytest.param("NS(N)(=O)=O", "sulfuric diamide", id="sulfuric_diamide"),
        pytest.param("C[C@H](CC)NS(=O)(=O)N", "N-[(2R)-butan-2-yl]sulfuric diamide", id="stereocentre_in_an_n_substituent"),
        pytest.param("C/C=C/CNS(=O)O", "N-[(2E)-but-2-en-1-yl]sulfuramidous acid", id="double_bond_in_an_n_substituent"),
        pytest.param("CN(C)S(=O)(=O)N", "N,N-dimethylsulfuric diamide", id="substituents_on_one_nitrogen"),
        pytest.param("CNS(=O)(=O)N(C)C", "N,N,N'-trimethylsulfuric diamide", id="substituents_on_both_nitrogens"),
        pytest.param("CNS(N)=O", "N-methylsulfurous diamide", id="sulfurous_diamide"),
        pytest.param("CN(Cl)S(=O)(=O)O", "N-chloro-N-methylsulfamic acid", id="halogen_on_the_nitrogen"),
        pytest.param("CN(N=O)S(=O)(=O)NC", "N,N'-dimethyl-N-nitrososulfuric diamide", id="nitroso_on_one_nitrogen_of_a_diamide"),
        pytest.param("CN(N=O)S(=O)(=O)Cl", "N-methyl-N-nitrososulfamoyl chloride", id="nitroso_on_the_nitrogen_of_an_acid_chloride"),
        pytest.param("CN([N+](=O)[O-])S(=O)(=O)O", "N-methyl-N-nitrosulfamic acid", id="nitro_on_the_nitrogen"),
    ],
)
def test_sulfamic_acid_and_the_amides_of_sulfuric_and_sulfurous_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-72.2.2.2.1.1, P-67.1.2.6: the anion of the amide acid of sulfurous acid ends in 'ite'; the acid, its ester and its chloride keep 'ous acid', 'ite' and 'ous chloride'
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NS(=O)O", "sulfuramidous acid", id="acid"),
        pytest.param("CNS(=O)O", "N-methylsulfuramidous acid", id="n_substituted_acid"),
        pytest.param("NS(=O)OC", "methyl sulfuramidite", id="ester"),
        pytest.param("NS(=O)Cl", "sulfuramidous chloride", id="chloride"),
        pytest.param("NS(=O)[O-]", "sulfuramidite", id="anion"),
        pytest.param("CNS(=O)[O-]", "N-methylsulfuramidite", id="n_methyl_anion"),
        pytest.param("CN(C)S(=O)[O-]", "N,N-dimethylsulfuramidite", id="two_groups_on_the_nitrogen"),
        pytest.param("CCNS(=O)[O-]", "N-ethylsulfuramidite", id="n_ethyl_anion"),
        pytest.param("c1ccccc1NS(=O)[O-]", "N-phenylsulfuramidite", id="n_phenyl_anion"),
        pytest.param("CNS(=O)[O-].[Na+]", "sodium N-methylsulfuramidite", id="anion_in_a_salt"),
    ],
)
def test_sulfuramidous_acid_and_its_anion(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "CS(=O)(=O)S(=O)(=O)c1ccccc1", "1-methyl-2-phenyl-1λ6,2λ6-disulfane-1,1,2,2-tetrone", id="disulfane_tetrone"
        ),
        pytest.param("CSCl", "methyl thiohypochlorite", id="thiohypochlorite_ester"),
        pytest.param("ClSCCSCl", "ethane-1,2-diyl bis(thiohypochlorite)", id="thiohypochlorite_multiplicative"),
        pytest.param("ClNCCSCl", "[2-(chlorosulfanyl)ethyl]hypochlorous amide", id="halogen_amide_outranks_thiohypohalite"),
        pytest.param(
            "OS(I)(I)SC#N", "1-hydroxy-1,1-diiodo-1λ4-disulfane-2-carbonitrile", id="lambda_chain_nitrile"
        ),
        pytest.param("CS(C)(C)SC", "1,1,1,2-tetramethyl-1λ4-disulfane", id="lambda_chain_without_suffix"),
        pytest.param("S=S", "disulfene", id="chain_double_bond_standard_valence"),
        pytest.param("CS(C)(C)(C)=S", "1,1,1,1-tetramethyl-1λ6-disulfene", id="chain_double_bond_lambda_label"),
        pytest.param("C[Se](C)=[Se]", "1,1-dimethyl-1λ4-diselene", id="chain_double_bond_selenium"),
        pytest.param("CSS(C)(C)=S", "2,2,3-trimethyl-2λ6-trisulf-1-ene", id="chain_double_bond_locant_before_prefixes"),
        pytest.param("CS(=S)=S", "2-methyl-2λ6-trisulfa-1,2-diene", id="chain_cumulated_double_bonds"),
        pytest.param("CS(C)(=O)=S", "1,1-dimethyl-1λ6-disulfen-1-one", id="chain_double_bond_beside_heterone_oxygen"),
        # P-64.7.2: the heterone is the parent when it holds more doubly bonded oxygen atoms than the ketone units
        pytest.param(
            "O=C1CCC(CC1)S(=O)(=O)S(=O)(=O)C1CCC(=O)CC1",
            "1,2-bis(4-oxocyclohexyl)-1λ6,2λ6-disulfane-1,1,2,2-tetrone",
            id="heterone_outranks_ring_ketone_units",
        ),
        pytest.param(
            "C1CCCCC1S(=O)(=O)S(=O)(=O)C1CCCCC1",
            "1,2-dicyclohexyl-1λ6,2λ6-disulfane-1,1,2,2-tetrone",
            id="heterone_outranks_ring_units_without_a_group",
        ),
        pytest.param(
            "OC1CCC(CC1)S(=O)(=O)S(=O)(=O)C1CCC(O)CC1",
            "1,2-bis(4-hydroxycyclohexyl)-1λ6,2λ6-disulfane-1,1,2,2-tetrone",
            id="heterone_outranks_hydroxy_units",
        ),
        pytest.param(
            "O=C1CC(=O)CC(C1)S(=O)S(=O)C1CC(=O)CC(=O)C1",
            "5,5'-(disulfinyl)di(cyclohexane-1,3-dione)",
            id="equal_oxygen_count_leaves_the_ring_units_as_parent",
        ),
        pytest.param(
            "CC(=O)CS(=O)(=O)S(=O)(=O)CC(C)=O",
            "1,2-bis(2-oxopropyl)-1λ6,2λ6-disulfane-1,1,2,2-tetrone",
            id="heterone_outranks_chain_ketone_arms",
        ),
        pytest.param(
            "CCS(=O)(=O)S(=O)(=O)CC",
            "1,2-diethyl-1λ6,2λ6-disulfane-1,1,2,2-tetrone",
            id="heterone_outranks_chain_arms_without_a_group",
        ),
        pytest.param(
            "OC(=O)CS(=O)(=O)S(=O)(=O)CC(O)=O",
            "2,2'-(disulfonyl)diacetic acid",
            id="carboxylic_acid_arms_outrank_the_heterone",
        ),
    ],
)
def test_chalcogen_chain_heterones_and_halogen_acid_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("SO", "sulfanol", id="hydroxy_outranks_thiol"),
        pytest.param("CSSO", "methyldisulfanol", id="identical_run_kept_whole_under_hydroxy"),
        pytest.param("c1ccccc1SSSO", "phenyltrisulfanol", id="organic_end_on_a_homogeneous_run"),
        pytest.param("COOS", "methyldioxidanethiol", id="thiol_on_dioxidane"),
        pytest.param("COOSS", "methyldioxidanedithioperoxol", id="two_atom_suffix_of_identical_atoms"),
        pytest.param("CSSOO", "methyldisulfaneperoxol", id="peroxol_suffix_after_an_identical_run"),
        pytest.param("COS[Se][SeH]", "methoxysulfanediselenoperoxol", id="alkoxy_prefix_beside_a_peroxol_suffix"),
        pytest.param("C[Se][Se][Se]OS", "methyltriselane-OS-thioperoxol", id="mixed_pair_suffix_with_locants"),
        pytest.param("CSOSO", "[(methylsulfanyl)oxy]sulfanol", id="hydroxy_outranks_thioperoxol"),
        pytest.param("[SeH][Te][SeH]", "tellanediselenol", id="same_suffix_at_both_ends"),
        pytest.param("OOSO", "(hydroperoxy)sulfanol", id="junior_peroxol_end_as_a_prefix"),
        pytest.param("OOSS", "disulfaneperoxol", id="peroxol_outranks_dithioperoxol"),
        pytest.param("[TeH][Te][Te][Se][SeH]", "tritellane(diselenoperoxol)", id="suffix_in_parentheses_without_prefixes"),
        pytest.param("OSOSO", "dithioxanediol", id="alternating_parent_with_two_suffixes"),
        pytest.param("SSSCSSS", "1,1'-methylenebis(trisulfane)", id="multiplicative_chalcogen_chain_parents"),
        pytest.param("SSSCCSSS", "1,1'-(ethane-1,2-diyl)bis(trisulfane)", id="multiplicative_chain_linker_with_locants"),
        pytest.param("C(CSSSC)SSSC", "1,1'-(ethane-1,2-diyl)bis(3-methyltrisulfane)", id="trisulfane_is_no_heterounit_of_skeletal_replacement"),
        pytest.param("CSSSCSSSC", "1,1'-methylenebis(3-methyltrisulfane)", id="substituent_on_the_far_end_of_each_polysulfane"),
        pytest.param("CSCSSSCSC", "bis[(methylsulfanyl)methyl]trisulfane", id="other_sulfur_units_stay_in_the_prefixes_of_the_polysulfane"),
        pytest.param("CSSSCCSSSCC", "ethyl[2-(methyltrisulfanyl)ethyl]trisulfane", id="different_end_groups_leave_a_substitutive_name"),
        pytest.param("ClCSSSCCSSSCCl", "1,1'-(ethane-1,2-diyl)bis[3-(chloromethyl)trisulfane]", id="compound_end_group_nests_the_marks"),
        pytest.param("ClSSSC", "chloro(methyl)trisulfane", id="halogen_end_of_a_homogeneous_run"),
        pytest.param("ClSSBr", "bromo(chloro)disulfane", id="pair_with_halogen_ends_only_is_the_parent"),
        pytest.param("OCCSSSCl", "2-(chlorotrisulfanyl)ethan-1-ol", id="hydroxy_on_an_arm_outranks_the_run_as_parent"),
    ],
)
def test_hydrogen_terminated_chalcogen_chains_are_parent_hydrides_with_suffixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[AsH5]", "λ5-arsane", id="pnictogen_hypervalent"),
        pytest.param("[IH3]", "λ3-iodane", id="halogen_hypervalent"),
        pytest.param("C[SH3]", "methyl-λ4-sulfane", id="substituted_hypervalent"),
        pytest.param("C[Po]CC", "ethyl(methyl)polane", id="polonium_and_astatine_standard_bonding_number"),
        pytest.param("c1ccccc1[IH2]", "phenyl-λ3-iodane", id="terminal_halogen_with_added_hydrogens"),
        pytest.param("SSS", "trisulfane", id="homogeneous_chalcogen_chain"),
        pytest.param("S[SH2]S", "2λ4-trisulfane", id="chain_lambda_locant"),
        pytest.param("OO", "dioxidane", id="dioxidane"),
        pytest.param("CC(=O)OI(OC(C)=O)c1ccccc1", "bis(acetyloxy)(phenyl)-λ3-iodane", id="acyloxy_groups_on_centre"),
        pytest.param(
            "OI(OS(=O)(=O)c1ccc(C)cc1)c1ccccc1",
            "[(4-methylbenzene-1-sulfonyl)oxy](phenyl)-λ3-iodanol",
            id="sulfonyloxy_beside_centre_hydroxy",
        ),
        pytest.param("c1ccc(I(c2ccccc2)OS(=O)(=O)C)cc1", "[(methanesulfonyl)oxy]di(phenyl)-λ3-iodane", id="sulfonyloxy_on_centre"),
        pytest.param("CC(=O)O[IH2]", "(acetyloxy)-λ3-iodane", id="acyloxy_on_halogen_with_added_hydrogens"),
        pytest.param("COI(C)C", "methoxydi(methyl)-λ3-iodane", id="alkoxy_on_centre"),
        pytest.param("C=S(C)C", "dimethyl(methylidene)-λ4-sulfane", id="ylidene_on_centre"),
        pytest.param("CC(=O)[IH2]", "acetyl-λ3-iodane", id="acyl_on_halogen_with_added_hydrogens_is_no_acyl_halide"),
        pytest.param("CC(=O)I(C)C", "acetyldi(methyl)-λ3-iodane", id="acyl_on_a_substituted_hypervalent_halogen"),
        pytest.param("C/C=C/S(C)(C)C", "trimethyl[(1E)-prop-1-en-1-yl]-λ4-sulfane", id="alkenyl_on_centre_cites_e"),
        pytest.param("C/C=C\\S(C)(C)C", "trimethyl[(1Z)-prop-1-en-1-yl]-λ4-sulfane", id="alkenyl_on_centre_cites_z"),
        pytest.param("C=S1CCCC1", "1-methylidene-1λ4-thiolane", id="ylidene_on_ring_centre"),
        pytest.param("OC(C[SH2]C)CSC", "1-(methyl-λ4-sulfanyl)-3-(methylsulfanyl)propan-2-ol", id="nonstandard_prefix_cited_first"),
        pytest.param("OC(C[SH2]C)C[SH4]C", "1-(methyl-λ6-sulfanyl)-3-(methyl-λ4-sulfanyl)propan-2-ol", id="higher_bonding_number_cited_first"),
        pytest.param("C=[Se]1C=CC=C1", "1-methylidene-1H-1λ4-selenophene", id="ylidene_on_mancude_ring_centre"),
        pytest.param("C=S1C=CC=C1C(=O)O", "1-methylidene-1H-1λ4-thiophene-2-carboxylic acid", id="ylidene_on_mancude_ring_centre_with_acid"),
    ],
)
def test_nonstandard_bonding_number_hydrides_and_chalcogen_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OC(=O)CI(O)O", "(dihydroxy-λ3-iodanyl)acetic acid", id="centre_hydroxy_under_acid"),
        pytest.param("OC(=O)C[IH2]", "(λ3-iodanyl)acetic acid", id="halogen_with_added_hydrogens"),
        pytest.param("OC(=O)C[PH4]", "(λ5-phosphanyl)acetic acid", id="pnictogen_hydrogens_counted"),
        pytest.param("OC(=O)CC[SH4]C", "3-(methyl-λ6-sulfanyl)propanoic acid", id="chalcogen_lambda6"),
        pytest.param("C=S(C)CC(=O)O", "[methyl(methylidene)-λ4-sulfanyl]acetic acid", id="ylidene_on_centre_under_acid"),
        pytest.param("C=[Se](C)(=O)CC(=O)O", "[methyl(methylidene)(oxo)-λ6-selanyl]acetic acid", id="oxo_and_ylidene_on_centre_under_acid"),
        pytest.param(
            "CC(=O)OI(OC(C)=O)c1ccc(cc1)C(O)=O", "4-[bis(acetyloxy)-λ3-iodanyl]benzoic acid", id="acyloxy_on_centre"
        ),
        pytest.param(
            "OCCI(OC(C)=O)OC(C)=O", "2-[bis(acetyloxy)-λ3-iodanyl]ethan-1-ol", id="acyloxy_on_centre_under_alcohol"
        ),
        pytest.param("NCC[SH2]C", "2-(methyl-λ4-sulfanyl)ethan-1-amine", id="centre_without_hydroxy_under_amine"),
        pytest.param("OCC[SH2]C", "2-(methyl-λ4-sulfanyl)ethan-1-ol", id="centre_without_hydroxy_under_alcohol"),
        pytest.param("OCI(O)O", "(hydroxymethyl)-λ3-iodanediol", id="centre_with_more_hydroxy_groups_is_parent"),
        pytest.param("OC(=O)C[SH2]CC(=O)O", "2,2'-(λ4-sulfanediyl)diacetic acid", id="lambda_linker"),
        pytest.param("OC(=O)C[SH](C)CC(=O)O", "2,2'-(methyl-λ4-sulfanediyl)diacetic acid", id="substituted_lambda_linker"),
    ],
)
def test_lambda_centres_are_prefixes_or_linkers_beside_senior_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CP(O)(O)(O)O", "methyl-λ5-phosphanetetrol", id="hydroxy_groups_beyond_three"),
        pytest.param("O[Te](O)(O)(O)(O)O", "λ6-tellanehexol", id="final_a_of_the_prefix_elided_before_ol"),
        pytest.param("NS(N)(N)(N)N", "λ6-sulfanepentamine", id="amino_groups_beyond_three"),
    ],
)
def test_lambda_hydrides_with_many_principal_groups_multiply_the_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OS(=O)(=O)OCCO", "2-hydroxyethyl hydrogen sulfate", id="sulfate_beside_alcohol"),
        pytest.param("COS(=O)(=O)OCCO", "2-hydroxyethyl methyl sulfate", id="mixed_sulfate_beside_alcohol"),
        pytest.param("OS(=O)OCC(=O)C", "2-oxopropyl hydrogen sulfite", id="sulfite_beside_ketone"),
        pytest.param("OS(=O)(=O)Oc1ccc(N)cc1", "4-aminophenyl hydrogen sulfate", id="sulfate_beside_amine"),
        pytest.param("NS(=O)(=O)OCC=O", "2-oxoethyl sulfamate", id="sulfamate_beside_aldehyde"),
        pytest.param("NS(=O)(=O)OCC(=O)Cl", "2-chloro-2-oxoethyl sulfamate", id="sulfamate_beside_acid_halide"),
        pytest.param("OS(=O)(=O)OCC(=O)OC", "methyl (sulfooxy)acetate", id="sulfate_yields_to_carboxylic_ester"),
        pytest.param("c1ccncc1COS(N)(=O)=O", "(pyridin-3-yl)methyl sulfamate", id="sulfamate_of_a_ring_nitrogen_alcohol_part"),
        pytest.param("Cc1ncsc1COS(N)(=O)=O", "(4-methyl-1,3-thiazol-5-yl)methyl sulfamate", id="sulfamate_of_a_ring_chalcogen_alcohol_part"),
    ],
)
def test_noncarbon_acid_esters_outrank_junior_groups_of_the_alcohol_component(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OS(=O)(=O)OCCC(O)=O", "3-(sulfooxy)propanoic acid", id="sulfo"),
        pytest.param("COS(=O)OCCC(O)=O", "3-[(methoxysulfinyl)oxy]propanoic acid", id="alkoxysulfinyl"),
        pytest.param("ClS(=O)(=O)OCCC(O)=O", "3-[(chlorosulfonyl)oxy]propanoic acid", id="halosulfonyl"),
        pytest.param("NS(=O)(=O)OCCC(O)=O", "3-(sulfamoyloxy)propanoic acid", id="sulfamoyl"),
        pytest.param("OS(=O)(=O)SCCC(O)=O", "3-(sulfosulfanyl)propanoic acid", id="attached_through_sulfur"),
        pytest.param("COS(=O)(=O)NCCC(O)=O", "3-[(methoxysulfonyl)amino]propanoic acid", id="attached_through_nitrogen"),
        pytest.param("OC(=O)CCOS(=O)(=O)C#N", "3-[(cyanosulfonyl)oxy]propanoic acid", id="cyanide_on_sulfur"),
        pytest.param(
            "OC(=O)CCOS(=O)(=O)N=C=S", "3-[(isothiocyanatosulfonyl)oxy]propanoic acid", id="pseudohalide_nitrogen_not_sulfamoyl"
        ),
        pytest.param("OS(=O)(=O)NCCC(O)=O", "3-(sulfoamino)propanoic acid", id="sulfamic_acid_yields_to_carboxylic_acid"),
    ],
)
def test_sulfur_acid_groups_attached_through_oxygen_under_a_carboxylic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[SH2]CCSCCSCCSC", "2λ4,5,8,11-tetrathiadodecane"),
        ("CSCC[SH4]CC[SH2]CCSC", "2,5λ6,8λ4,11-tetrathiadodecane"),
        ("C[SH2]CCOCCOCCOC", "2,5,8-trioxa-11λ4-thiadodecane"),
        ("COCCOCC[SH2]CCOC", "2,5,11-trioxa-8λ4-thiadodecane"),
        ("O1CC[SH2]CCCCCCCCCC1", "1-oxa-4λ4-thiacyclotetradecane"),
        ("O1CC[SH4]CCCCCCC[SH2]CC1", "1-oxa-4λ6,12λ4-dithiacyclotetradecane"),
        ("O1C[SH2]CC1", "1,3λ4-oxathiolane"),
        ("[SH2]1C=CC=CC=C1", "1H-1λ4-thiepine"),
    ],
)
def test_skeletal_chalcogen_with_nonstandard_bonding_number_keeps_lambda(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se][Se]CCSSC", "1-(methyldiselanyl)-2-(methyldisulfanyl)ethane", id="mixed_diselanyl_and_disulfanyl_on_ethane"),
        pytest.param("CC(C)[Se][Se]CCC", "1-[(propan-2-yl)diselanyl]propane", id="branched_diselanyl_prefix"),
        pytest.param("CSSCCSSC", "1,2-bis(methyldisulfanyl)ethane", id="two_disulfanyl_prefixes"),
        pytest.param("BrSSC", "(bromodisulfanyl)methane", id="halogen_on_the_far_end_of_a_disulfanyl_prefix"),
        pytest.param("I[Te][Te]CC", "(iododitellanyl)ethane", id="halogen_on_a_ditellanyl_prefix"),
        pytest.param("CSOCC", "[(methylsulfanyl)oxy]ethane", id="sulfur_oxygen_pair_prefix"),
        pytest.param("COSC1CCCCC1", "(methoxysulfanyl)cyclohexane", id="alkoxysulfanyl_prefix_on_a_ring"),
        pytest.param("C1(=CC=CC=C1)[Se][Te]C1=CC=CC=C1", "[(phenylselanyl)tellanyl]benzene", id="selenium_tellurium_pair_prefix"),
        pytest.param("CCCCCSS[C@H](C)CC", "1-{[(2R)-butan-2-yl]disulfanyl}pentane", id="stereocentre_in_a_disulfanyl_prefix"),
        pytest.param("CSCSSCCSCCSC", "2,4,5,8,11-pentathiadodecane", id="four_units_take_skeletal_replacement"),
        pytest.param("CSCSSCCSCC[Se]C", "2,4,5,8-tetrathia-11-selenadodecane", id="four_units_with_selenium"),
    ],
)
def test_chains_of_two_chalcogens_are_prefixes_and_four_units_make_a_replacement_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OB(O)C(C)CSS", "(1-disulfanylpropan-2-yl)boronic acid", id="boronic_acid"),
        pytest.param("C[N+](C)(C)CCSS", "2-disulfanyl-N,N,N-trimethylethan-1-aminium", id="quaternary_ammonium"),
        pytest.param("[NH3+]CC[Se][SeH]", "2-diselanylethan-1-aminium", id="primary_ammonium"),
        pytest.param("C[P+](C)(C)CCSS", "(2-disulfanylethyl)tri(methyl)phosphanium", id="phosphonium"),
        pytest.param("SC(=O)CCSS", "3-disulfanylpropanethioic S-acid", id="thioic_acid"),
        pytest.param("NCCSS", "2-aminoethane-1-dithioperoxol", id="amine_ranks_below_the_dithioperoxol"),
    ],
)
def test_a_dithioperoxol_group_is_a_prefix_only_beside_a_senior_class(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereogenic_sulfinyl_group_cites_its_descriptor():
    assert smiles_to_iupac("C[S@](=O)c1ccccc1") == "(S)-(methanesulfinyl)benzene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C/C=C/S(=O)(=O)c1ccccc1", "[(1E)-prop-1-ene-1-sulfonyl]benzene", id="alkenesulfonyl_on_a_ring"),
        pytest.param("C/C=C\\S(=O)c1ccccc1", "[(1Z)-prop-1-ene-1-sulfinyl]benzene", id="alkenesulfinyl_on_a_ring"),
        pytest.param("C/C=C/S(=O)(=O)CC(=O)O", "[(1E)-prop-1-ene-1-sulfonyl]acetic acid", id="alkenesulfonyl_on_a_chain"),
    ],
)
def test_double_bond_of_a_sulfonyl_or_sulfinyl_prefix_cites_its_descriptor(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[S@](=O)CC", "(S)-(methanesulfinyl)ethane"),
        ("[N+](=O)([O-])C1=CC=C(C=C1)[S@](=O)OCC", "ethyl (R)-4-nitrobenzene-1-sulfinate"),
        ("CC[C@H](C)[S@](=O)OC", "methyl (2S,R)-butane-2-sulfinate"),
    ],
)
def test_stereogenic_sulfoxide_and_sulfinate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OS1(=NCCC1)=O", "1-hydroxy-4,5-dihydro-3H-1λ6,2-thiazol-1-one"),
        ("OS1(=NCCCCCCCCCC1)=O", "1-hydroxy-1λ6-thia-2-azacyclododec-1-en-1-one"),
        ("S1(=NCCCC1)O", "3,4,5,6-tetrahydro-1λ4,2-thiazin-1-ol"),
        ("O=S1(=O)NC=CC=C1", "1λ6,2-thiazine-1,1(2H)-dione"),
    ],
)
def test_ring_chalcogen_with_oxo_or_hydroxy_in_an_unsaturated_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=S1NCCCCCCCCC1", "1λ4-thia-2-azacycloundecan-1-one", id="sultim_of_eleven_ring_members"),
    ],
)
def test_sultim_above_ten_ring_members_is_named_by_skeletal_replacement(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("S1=CN=CC=C1", "1λ4,3-thiazine"),
        ("[SiH]1=CC=CC=C1", "siline"),
    ],
)
def test_lambda_and_silicon_mancude_monocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OS(=O)(=O)Cl", "sulfurochloridic acid"),
        ("NS(=O)(=O)S", "sulfamothioic S-acid"),
        ("OS(=O)(=O)[N+]#[C-]", "sulfurisocyanidic acid"),
        ("FS(=N)N=C=O", "sulfurimidisocyanatidous fluoride"),
        ("CNS(=O)(Br)=NC", "N,N'-dimethylsulfuramidimidic bromide"),
        ("CNS(=O)(=O)Cl", "N-methylsulfamoyl chloride"),
        ("O=S(=O)(Cl)N1CCOCC1", "morpholine-4-sulfonyl chloride"),
        ("COS(=O)(=O)N1CCCCC1", "methyl piperidine-1-sulfonate"),
        ("CN=S(Cl)Cl", "N-methylsulfurimidous dichloride"),
        ("CN(C)S(=O)(=NC1=CC=CC=C1)N(C)C", "N,N,N',N'-tetramethyl-N''-phenylsulfurimidic diamide"),
        ("CN=S(O)S(=N)O", "N1-methyl-1,2-diimidodithionous acid"),
    ],
)
def test_sulfur_oxoacids_modified_by_infixes_and_their_halides_and_amides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)C[Cl+][S-]", "(thiochlorosyl)acetic acid"),
        ("[S-][Cl+]c1ccccc1", "(thiochlorosyl)benzene"),
    ],
)
def test_thiohalogen_oxo_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("SSSSC(=O)c1ccccc1", "phenyl(tetrasulfanyl)methanone"),
        ("SS(=O)(=O)SSSS(=O)(=O)S", "trisulfanedisulfonothioic S-acid"),
        ("SSSS(=O)(=O)O", "trisulfanesulfonic acid"),
    ],
)
def test_sulfur_chain_broken_to_express_the_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(Oc2ccccc2[SH5])cc1", "4-[2-(λ6-sulfanyl)phenoxy]benzoic acid"),
        ("OC(=O)C(C[SH3])C[SH5]", "3-(λ6-sulfanyl)-2-[(λ4-sulfanyl)methyl]propanoic acid"),
        (
            "OC(=O)C(CCCP[PH4])CCC[PH3]P",
            "5-(1λ5-diphosphan-1-yl)-2-[3-(2λ5-diphosphan-1-yl)propyl]pentanoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)C[SH2]COC(C[SH4]CSC)C[SH2]C[SH2]C",
            "4-[5-({[(methylsulfanyl)methyl]-λ6-sulfanyl}methyl)-4-oxa-2λ4,7λ4,9λ4-trithiadecan-1-yl]benzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)C[SiH2]C[SiH2]C[SiH2]C(CCC)[SiH2]C",
            "4-(7-propyl-2,4,6,8-tetrasilanonan-1-yl)benzoic acid",
        ),
        ("OC(=O)c1ccc(cc1)[SiH2]C[SiH2]C[SiH2]C[SiH2]C", "4-(1,3,5,7-tetrasilaoctan-1-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)[SiH](S[SiH3])O[SiH3]", "4-[1-(silylsulfanyl)disiloxanyl]benzoic acid"),
        ("OC(=O)c1ccc(cc1)[SiH](O[SiH3])O[SiH3]", "4-[1-(silyloxy)disiloxanyl]benzoic acid"),
    ],
)
def test_nonstandard_bonding_in_substituents_and_substituent_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[SH]1=CCC=C1", "3H-1λ4-thiophene"),
        ("[SH3]1=CC=CC=C1", "1λ6-thiopyran"),
        ("[SH2]1C=Cc2ccccc12", "1H-1λ4-benzothiophene"),
        ("[IH]1C=Cc2ccccc12", "1H-1λ3-benziodole"),
        ("[IH]1CCCCC1", "1λ3-iodinane"),
        ("[IH3]1C=CC=C1", "1H-1λ5-iodole"),
        ("C1=CC=C[I]=C1", "1λ3-iodinine"),
        ("COI1OC(=O)c2ccccc12", "1-methoxy-1λ3,2-benziodoxol-3(1H)-one"),
        ("O1C[SH4]CC1C1C[SH2]CO1", "5-(1,3λ4-oxathiolan-5-yl)-1,3λ6-oxathiolane"),
        ("O1CC=CC2=[PH2]C=CC=C12", "2H-5λ5-phosphinino[3,2-b]pyran"),
        ("O1C=CC2=[PH2]C=CC=C12", "4λ5-phosphinino[3,2-b]furan"),
    ],
)
def test_indicated_hydrogen_and_halogen_atoms_of_nonstandard_bonding_in_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "CS(=O)(=O)N(CCOc1cc(Cl)cc(Cl)c1)CCc1cccc(CC(=O)O)c1",
            "[3-(2-{[2-(3,5-dichlorophenoxy)ethyl](methanesulfonyl)amino}ethyl)phenyl]acetic acid",
            id="sulfonyl_and_alkyl_on_the_nitrogen_of_an_amino_prefix",
        ),
        pytest.param("CC(=O)N(S(C)(=O)=O)CCC(=O)O", "3-[N-(methanesulfonyl)acetamido]propanoic acid", id="acyl_and_sulfonyl"),
    ],
)
def test_amino_prefix_with_a_sulfonyl_group_on_the_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CCS(CC)=Nc1ccccc1", "S,S-diethyl-N-phenyl-\u03bb4-sulfanimine", id="sulfimide_cites_s_and_n_locants"),
        pytest.param("C[Se](C)=Nc1ccccc1", "Se,Se-dimethyl-N-phenyl-\u03bb4-selanimine", id="selenium_analogue_cites_its_symbol"),
        pytest.param("O=S=Nc1ccccc1", "(phenylimino)-\u03bb4-sulfanone", id="sulfinylamine"),
        pytest.param("CN=S(=O)=O", "(methylimino)-\u03bb6-sulfanedione", id="sulfonylamine"),
        pytest.param("CS(C)(=O)=Nc1ccccc1", "dimethyl(phenylimino)-\u03bb6-sulfanone", id="sulfoximide"),
        pytest.param("CS(=O)(=Nc1ncccc1)C", "dimethyl[(pyridin-2-yl)imino]-\u03bb6-sulfanone", id="heteroaryl_on_the_imino_nitrogen"),
        pytest.param("CS(=O)(=Nc1nccs1)C", "dimethyl[(1,3-thiazol-2-yl)imino]-\u03bb6-sulfanone", id="heteroaryl_with_two_ring_heteroatoms"),
        pytest.param("COCCN=S(=O)(C)C", "[(2-methoxyethyl)imino]di(methyl)-\u03bb6-sulfanone", id="ether_in_the_imino_substituent"),
        pytest.param("CSCCN=S(=O)(C)C", "dimethyl{[2-(methylsulfanyl)ethyl]imino}-\u03bb6-sulfanone", id="sulfide_in_the_imino_substituent"),
        pytest.param("CS(=O)(=NO)C", "(hydroxyimino)di(methyl)-\u03bb6-sulfanone", id="hydroxy_on_the_imino_nitrogen"),
        pytest.param("CS(=O)(=NN)C", "(aminoimino)di(methyl)-\u03bb6-sulfanone", id="amino_on_the_imino_nitrogen"),
        pytest.param("N=S(=N)(c1ccccc1)c1ccccc1", "diphenyl-\u03bb6-sulfanediimine", id="sulfonediimine_with_sulfur_groups_only"),
        pytest.param("CN=S=NCC", "ethyl(methyl)-\u03bb4-sulfanediimine", id="sulfur_diimide_without_locants"),
        pytest.param("N=S(=N)=N", "\u03bb6-sulfanetriimine", id="sulfur_triimide"),
        pytest.param("C[S@](=O)(=N)c1ccccc1", "(R)-imino(methyl)(phenyl)-\u03bb6-sulfanone", id="stereogenic_sulfur_of_a_sulfoximide"),
        pytest.param("C[C@@H](Cl)[S@](C)(=N)=O", "(R)-[(1R)-1-chloroethyl](imino)(methyl)-\u03bb6-sulfanone", id="stereocentre_in_a_group_of_a_stereogenic_sulfur"),
    ],
)
def test_sulfimides_sulfoximides_and_sulfur_imides_on_the_lambda_sulfane_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "c1ccccc1S(=O)(=O)N=S(c1ccccc1)c1ccccc1",
            "N-(diphenyl-\u03bb4-sulfanylidene)benzenesulfonamide",
            id="sulfimide_under_a_sulfonamide",
        ),
        pytest.param(
            "CS(C)(=O)=NS(=O)(=O)c1ccccc1",
            "N-[dimethyl(oxo)-\u03bb6-sulfanylidene]benzenesulfonamide",
            id="sulfoximide_under_a_sulfonamide",
        ),
        pytest.param("CS(=O)(=NC#N)C", "[dimethyl(oxo)-\u03bb6-sulfanylidene]cyanamide", id="sulfoximide_under_a_cyanamide"),
        pytest.param("N#CN=S(C)C", "(dimethyl-\u03bb4-sulfanylidene)cyanamide", id="sulfimide_under_a_cyanamide"),
        pytest.param(
            "CS(=O)(=N[N+](=O)[O-])C", "[dimethyl(oxo)-\u03bb6-sulfanylidene]nitramide", id="sulfoximide_under_a_nitramide"
        ),
        pytest.param(
            "CS(C)(=O)=NC(=O)NC", "N-[dimethyl(oxo)-\u03bb6-sulfanylidene]-N'-methylurea", id="sulfoximide_on_a_urea_nitrogen"
        ),
        pytest.param(
            "CS(C)(=O)=NC(=O)OC", "methyl [dimethyl(oxo)-\u03bb6-sulfanylidene]carbamate", id="sulfoximide_on_a_carbamate_nitrogen"
        ),
    ],
)
def test_chalcogen_imide_cited_as_a_ylidene_prefix_beside_a_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-16.3.5: a substituted ester group of sulfurous or sulfuric acid is multiplied with bis
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("ClCOS(=O)OCCl", "bis(chloromethyl) sulfite", id="sulfite"),
        pytest.param("ClC(Cl)OS(=O)(=O)OC(Cl)Cl", "bis(dichloromethyl) sulfate", id="sulfate_with_a_prefix_of_its_own_multiplier"),
    ],
)
def test_multiplier_of_identical_substituted_sulfur_ester_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-65.2.1.2, P-65.6.3.3.7.2.1, P-77.3.1: the anion of an organyl ester of a thiocarbonic acid is the ester name without 'hydrogen'
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CCOC(=S)[S-].[K+]", "potassium O-ethyl carbonodithioate", id="xanthate_salt"),
        pytest.param("CCOC(=S)[S-]", "O-ethyl carbonodithioate", id="lone_xanthate_anion"),
        pytest.param("CCSC(=S)[S-].[K+]", "potassium ethyl carbonotrithioate", id="trithiocarbonate_monoester"),
        pytest.param("CCSC(=O)[O-].[K+]", "potassium S-ethyl carbonothioate", id="s_ester_with_the_anion_on_oxygen"),
        pytest.param("CCOC(=O)[S-].[Na+]", "sodium O-ethyl carbonothioate", id="tautomer_of_the_o_ester"),
    ],
)
def test_anions_and_salts_of_thiocarbonic_acid_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
