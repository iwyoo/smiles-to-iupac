import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_ethyl_disulfide():
    assert smiles_to_iupac("CSSCC") == "(methyldisulfanyl)ethane"


def test_branched_disulfanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)SSC(C)C")


def test_tetrasulfide_chain():
    assert smiles_to_iupac("CSSSSC") == "(methyltetrasulfanyl)methane"


def test_terminal_persulfide_methane():
    assert smiles_to_iupac("CSS") == "disulfanylmethane"


def test_terminal_persulfide_ethane():
    assert smiles_to_iupac("CCSS") == "disulfanylethane"


def test_both_terminal_disulfane_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SS")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)SSCC", "(2S)-2-(ethyldisulfanyl)butane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCSS[C@H](C)CC")


def test_phenyl_disulfide_direct_bond():
    assert smiles_to_iupac("c1ccccc1SSC") == "(methyldisulfanyl)benzene"
    assert smiles_to_iupac("c1ccccc1SSCC") == "(ethyldisulfanyl)benzene"


def test_phenyl_disulfide_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CSSC")


def test_phenyl_disulfide_perthiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1SS")


def test_phenyl_disulfide_branched_other_side_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1SSC(C)C")


def test_phenyl_disulfide_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1SSC")


def test_phenyl_disulfide_unsaturation_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CSSC")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(N=C=S)CC1")


def test_phenyl_isothiocyanate_chain():
    assert smiles_to_iupac("c1ccccc1CN=C=S") == "(isothiocyanatomethyl)benzene"


def test_phenyl_isothiocyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1N=C=S")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCN=C=S")


def test_two_isothiocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C=NCN=C=S")


def test_sulfuric_acid_itself_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)O")


def test_sulfite_ester_unaffected():
    # One fewer double-bonded O routes to `_sulfite.py` instead, unchanged.
    assert smiles_to_iupac("COS(=O)OC") == "dimethyl sulfite"


def test_salt_of_partial_ester():
    assert smiles_to_iupac("COS(=O)(=O)[O-].[Na+]") == "sodium methyl sulfate"


def test_salt_of_partial_ester_multivalent_cation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COS(=O)(=O)[O-].[Ca+2]")


def test_methylsulfanylpropane():
    assert smiles_to_iupac("CSCCC") == "1-methylsulfanylpropane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)SC(C)C", "2-(propan-2-yl)sulfanylpropane"),
        ("CC(C)CSC(C)CC", "2-(2-methylpropyl)sulfanylbutane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)SCC", "(2S)-2-ethylsulfanylbutane"),
    ],
)
def test_stereocenter_on_parent_chain__sulfide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_substituent_branch():
    assert smiles_to_iupac("CCCCCS[C@H](C)CC") == "1-{[(2R)-butan-2-yl]sulfanyl}pentane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1SC", "methylsulfanylbenzene"),
        ("c1ccccc1SC(C)C", "(propan-2-yl)sulfanylbenzene"),
        ("c1ccccc1CSCC", "(ethylsulfanylmethyl)benzene"),
        ("c1ccccc1CSC(C)C", "[(propan-2-yl)sulfanylmethyl]benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents():
    assert smiles_to_iupac("CSc1ccccc1SC") == "1,2-bis(methylsulfanyl)benzene"


def test_benzene_ring_stereocenter():
    assert smiles_to_iupac("c1ccccc1S[C@H](C)CC") == "{[(2R)-butan-2-yl]sulfanyl}benzene"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C(S(=O)N)C")


def test_two_sulfinamides_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CS(=O)N")


def test_cyclohexanesulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCCC1") == "cyclohexanesulfinamide"


def test_2_chlorocyclohexane_1_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfinamide"


def test_unsaturated_ring_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCC=C1") == "cyclohex-2-ene-1-sulfinamide"
    assert smiles_to_iupac("O=S(N)C1CC=CCC1") == "cyclohex-3-ene-1-sulfinamide"


def test_unsaturated_ring_sulfinamide_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(N)C1CCCC=C1C")


def test_unsaturated_ring_sulfinamide_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(N)C1CCCC#C1")


def test_unsaturated_ring_sulfinamide_with_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(NC)C1CCCC=C1")


def test_ring_substituent_chain_sulfinamide():
    assert smiles_to_iupac("NS(=O)CC1CCCCC1") == "cyclohexylmethanesulfinamide"


def test_ring_substituent_chain_sulfinamide_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CC1CCC(C)CC1")


def test_sulfinamide_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CCO")


def test_benzenesulfinamide():
    assert smiles_to_iupac("c1ccccc1S(=O)N") == "benzenesulfinamide"


def test_substituted_benzenesulfinamide_n_alkyl():
    assert smiles_to_iupac("Cc1ccc(cc1)S(=O)NC") == "N,4-dimethylbenzenesulfinamide"


def test_phenyl_chain_sulfinamide_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCCS(=O)NC")


def test_phenyl_chain_sulfinamide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CS(=O)N)cc1") == "(4-chlorophenyl)methanesulfinamide"


def test_phenyl_chain_sulfinamide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCS(=O)N")


def test_n_methylmethanesulfinamide():
    assert smiles_to_iupac("CS(=O)NC") == "N-methylmethanesulfinamide"


def test_n_methylcyclohexanesulfinamide():
    assert smiles_to_iupac("O=S(NC)C1CCCCC1") == "N-methylcyclohexanesulfinamide"


def test_halogenated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)NCCCl")


def test_unsaturated_n_substituted_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)NCC=C")


def test_sulfinamide_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)S(=O)N")


def test_sulfinamide_specified_sulfur_stereocenter_cited():
    assert smiles_to_iupac("CCC[S@](=O)N") == "(S)-propane-1-sulfinamide"
    assert smiles_to_iupac("CCC[S@@](=O)N") == "(R)-propane-1-sulfinamide"


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


def test_unsaturated_ring_sulfinic_acid_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CCCC#C1")


def test_ring_substituent_chain_sulfinic_acid():
    assert smiles_to_iupac("OS(=O)CC1CCCCC1") == "cyclohexylmethanesulfinic acid"


def test_benzenesulfinic_acid():
    assert smiles_to_iupac("c1ccccc1S(=O)O") == "benzenesulfinic acid"  # CID 12057


def test_substituted_benzenesulfinic_acid():
    assert smiles_to_iupac("Cc1ccccc1S(=O)O") == "2-methylbenzene-1-sulfinic acid"  # CID 12661295


def test_phenyl_chain_sulfinic_acid_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCCS(=O)O)cc1") == "3-(4-chlorophenyl)propane-1-sulfinic acid"


def test_phenyl_chain_sulfinic_acid_ring_ethyl():
    assert smiles_to_iupac("CCc1ccc(cc1)CS(=O)O") == "(4-ethylphenyl)methanesulfinic acid"


def test_sulfurous_acid_itself_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)O")


def test_sulfate_ester_unaffected():
    # A second S=O bond routes to `_sulfate.py` instead, unchanged.
    assert smiles_to_iupac("COS(=O)(=O)OC") == "dimethyl sulfate"


def test_salt_of_partial_ester__sulfite():
    assert smiles_to_iupac("COS(=O)[O-].[Na+]") == "sodium methyl sulfite"


def test_salt_of_partial_ester_multivalent_cation_raises__sulfite():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COS(=O)[O-].[Ca+2]")


def test_cyclohexanesulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1") == "cyclohexanesulfonamide"


def test_n_methylcyclohexanesulfonamide():
    assert smiles_to_iupac("O=S(=O)(NC)C1CCCCC1") == "N-methylcyclohexanesulfonamide"


def test_unsaturated_ring_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC=C1") == "cyclohex-2-ene-1-sulfonamide"
    assert smiles_to_iupac("O=S(=O)(N)C1CC=CCC1") == "cyclohex-3-ene-1-sulfonamide"


def test_unsaturated_ring_sulfonamide_with_substituent():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-sulfonamide"


def test_unsaturated_ring_sulfonamide_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CCCC#C1")


def test_unsaturated_ring_sulfonamide_with_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(NC)C1CCCC=C1")


def test_ring_substituent_chain_sulfonamide():
    assert smiles_to_iupac("NS(=O)(=O)CC1CCCCC1") == "cyclohexylmethanesulfonamide"


def test_ring_substituent_chain_sulfonamide_ring_with_substituent():
    assert smiles_to_iupac("NS(=O)(=O)CC1CCC(C)CC1") == "(4-methylcyclohexyl)methanesulfonamide"


def test_benzenesulfonamide():
    assert smiles_to_iupac("c1ccccc1S(=O)(=O)N") == "benzenesulfonamide"


def test_substituted_benzenesulfonamide_n_alkyl():
    assert smiles_to_iupac("Cc1ccc(cc1)S(=O)(=O)NC") == "N,4-dimethylbenzenesulfonamide"


def test_phenyl_chain_sulfonamide_n_alkyl():
    assert smiles_to_iupac("c1ccccc1CCCS(=O)(=O)NC") == "N-methyl-3-phenylpropane-1-sulfonamide"


def test_phenyl_chain_sulfonamide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCCS(=O)(=O)N)cc1") == "3-(4-chlorophenyl)propane-1-sulfonamide"


def test_phenyl_chain_sulfonamide_ring_ethyl():
    assert smiles_to_iupac("CCc1ccc(cc1)CS(N)(=O)=O") == "(4-ethylphenyl)methanesulfonamide"


def test_phenyl_chain_sulfonamide_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CCS(=O)(=O)N") == "2-(2-ethenylphenyl)ethanesulfonamide"


def test_two_different_n_substituents_alphabetized_ignoring_italic_prefix():
    assert smiles_to_iupac("CS(=O)(=O)N(CC)C(C)(C)C") == "N-tert-butyl-N-ethylmethanesulfonamide"


def test_halogenated_n_substituent():
    assert smiles_to_iupac("CS(=O)(=O)NCCCl") == "N-(2-chloroethyl)methanesulfonamide"


def test_unsaturated_n_substituted_sulfonamide():
    assert smiles_to_iupac("CS(=O)(=O)NCC=C") == "N-(prop-2-en-1-yl)methanesulfonamide"


def test_acyclic_sulfonamide_stereocenter_with_coexisting_substituent():
    assert smiles_to_iupac("C[C@@H](Cl)S(=O)(=O)N") == "(1R)-1-chloroethanesulfonamide"


def test_cyclic_sulfonamide_stereocenter():
    assert (
        smiles_to_iupac("N[S](=O)(=O)[C@H]1CCCC[C@@H]1Cl")
        == "(1S,2S)-2-chlorocyclohexane-1-sulfonamide"
    )


def test_cyclic_sulfonamide_branch_stereocenter():
    assert (
        smiles_to_iupac("NS(=O)(=O)C1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-sulfonamide"
    )


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
    "smiles,expected",
    [
        ("CCS(=O)(=O)CC", "(ethanesulfonyl)ethane"),
        ("CS(=O)(=O)CCCC", "1-(methanesulfonyl)butane"),
    ],
)
def test_sulfone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported__sulfone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S1(=O)CCCCC1")


def test_unsaturated_ring_sulfonic_acid_with_substituent():
    assert smiles_to_iupac("OS(=O)(=O)C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-sulfonic acid"


def test_unsaturated_ring_sulfonic_acid_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC#C1")


def test_ring_substituent_chain_sulfonic_acid():
    assert smiles_to_iupac("OS(=O)(=O)CC1CCCCC1") == "cyclohexylmethanesulfonic acid"


def test_ring_substituent_chain_sulfonic_acid_ring_with_substituent():
    assert smiles_to_iupac("OS(=O)(=O)CC1CCC(C)CC1") == "(4-methylcyclohexyl)methanesulfonic acid"


def test_ring_substituent_chain_sulfonic_acid_unsaturated_ring():
    assert smiles_to_iupac("OS(=O)(=O)CC1CCCC=C1") == "(cyclohex-2-en-1-yl)methanesulfonic acid"


def test_benzenesulfonic_acid():
    assert smiles_to_iupac("c1ccccc1S(=O)(=O)O") == "benzenesulfonic acid"  # CID 7371


def test_phenyl_chain_sulfonic_acid_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CCS(=O)(=O)O") == "2-(2-ethenylphenyl)ethanesulfonic acid"


def test_heteroaromatic_direct_attachment_sulfonic_acid():
    assert smiles_to_iupac("OS(=O)(=O)c1cccnc1") == "pyridine-3-sulfonic acid"


def test_acyclic_sulfonic_acid_stereocenter_with_coexisting_substituent():
    assert smiles_to_iupac("C[C@@H](Cl)S(=O)(=O)O") == "(1R)-1-chloroethanesulfonic acid"


def test_cyclic_sulfonic_acid_stereocenter():
    assert (
        smiles_to_iupac("O=S(=O)(O)[C@H]1CCCC[C@@H]1Cl")
        == "(1S,2S)-2-chlorocyclohexane-1-sulfonic acid"
    )


def test_cyclic_sulfonic_acid_branch_stereocenter():
    assert (
        smiles_to_iupac("OS(=O)(=O)C1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-sulfonic acid"
    )


def test_sulfanium():
    assert smiles_to_iupac("[SH3+]") == "sulfanium"


def test_branched_substituent_not_supported__sulfonium():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[SH2+]")


def test_ring_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SH2+]C1CCCCC1")


def test_sulfonium_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[S@+](CC)CCC")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[S@@+](CC)CCC")


def test_mixed_alkyl_and_phenyl_substituents():
    assert smiles_to_iupac("C[S+](C)c1ccccc1") == "dimethyl(phenyl)sulfanium"


def test_substituted_phenyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[SH2+]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCCS(=O)CC", "1-(ethanesulfinyl)butane"),
        ("CCS(=O)CC", "(ethanesulfinyl)ethane"),
    ],
)
def test_sulfoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_resolves_via_hetero_ring_oxide():
    assert smiles_to_iupac("O=S1CCCCC1") == "thiane 1-oxide"


def test_thioate_branched_r_group():
    assert smiles_to_iupac("CC(C)C(=O)[S-]") == "2-methylpropanethioate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(=O)[S-]", "(2S)-2-methylbutanethioate"),
    ],
)
def test_thioate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_thioate_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CC(=O)[S-])cc1") == "2-(4-chlorophenyl)ethanethioate"


def test_two_thiocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#CSCCSC#N")


def test_phenyl_thiocyanate():
    assert smiles_to_iupac("c1ccccc1SC#N") == "phenyl thiocyanate"  # CID 21357


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=S)O", "ethanethioic O-acid"),
        ("O=CS", "methanethioic S-acid"),
    ],
)
def test_thioic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring():
    assert smiles_to_iupac("O=C1CCCCC1S") == "2-sulfanylcyclohexan-1-one"


def test_phenyl_substituent_thioic_acid_branched_chain():
    assert (
        smiles_to_iupac("CC(C)Cc1ccc(cc1)C(C)C(=O)S")
        == "2-[4-(2-methylpropyl)phenyl]propanethioic S-acid"
    )


def test_unsaturated_ring_thiol_with_substituent():
    assert smiles_to_iupac("SC1CCCC=C1C") == "2-methylcyclohex-2-ene-1-thiol"


def test_unsaturated_ring_thiol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC#C1")


def test_ring_substituent_chain_thiol_unsaturated_ring():
    assert smiles_to_iupac("SCC1CCCC=C1") == "(cyclohex-2-en-1-yl)methanethiol"


def test_ring_with_thiol_chain_thiol_tie():
    assert smiles_to_iupac("SC1CCCCC1CCS") == "2-(2-sulfanylethyl)cyclohexane-1-thiol"


def test_thiol_with_alcohol():
    assert smiles_to_iupac("SCCO") == "2-sulfanylethanol"


def test_benzenethiol():
    assert smiles_to_iupac("c1ccccc1S") == "benzenethiol"  # CID 7969


def test_two_direct_ring_thiols():
    assert smiles_to_iupac("Sc1ccccc1S") == "benzene-1,2-dithiol"


def test_phenyl_chain_thiol_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CCS") == "2-(2-ethenylphenyl)ethanethiol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)S", "(2S)-butane-2-thiol"),
    ],
)
def test_acyclic_thiol_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_thiol_branch_stereocenter():
    assert smiles_to_iupac("SC1(CCCCC1)[C@@H](C)CC") == "1-[(2S)-butan-2-yl]cyclohexane-1-thiol"


def test_thiol_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("S[C@H]1CCCCC1Cl") == "(1S)-2-chlorocyclohexane-1-thiol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("SCCc1cccnc1", "2-(pyridin-3-yl)ethanethiol"),
        ("SC1CCCCC1c1ccc[nH]1", "2-(1H-pyrrol-2-yl)cyclohexane-1-thiol"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_aromatic_substituent_thiol_three_rings():
    assert smiles_to_iupac("SC1CCCC(c2ccccc2)C1c3ccccc3") == "2,3-diphenylcyclohexane-1-thiol"


def test_two_ring_aromatic_substituent_thiol_chain_thiol():
    assert smiles_to_iupac("SCC1CCCCC1c1ccccc1") == "(2-phenylcyclohexyl)methanethiol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCC(=S)C", "1-chloropropane-2-thione"),
    ],
)
def test_thione_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_thione_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CCCC#C1")


def test_ring_substituent_chain_thione():
    assert smiles_to_iupac("CC(=S)C1CCCCC1") == "1-cyclohexylethanethione"


def test_ring_with_thione_chain_thione_tie():
    assert smiles_to_iupac("S=C1CCCCC1CC(=S)C") == "2-(2-sulfanylidenepropyl)cyclohexane-1-thione"


def test_thial_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=S")


def test_aromatic_thione_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc(cc1)C(=S)C")


def test_acyclic_thione_stereocenter():
    assert smiles_to_iupac("CC[C@@H](C)C(C)=S") == "(3R)-3-methylpentane-2-thione"


def test_cyclic_thione_stereocenter():
    assert smiles_to_iupac("S=C1CCCC[C@H]1C") == "(2R)-2-methylcyclohexane-1-thione"


def test_cyclic_thione_branch_stereocenter():
    assert smiles_to_iupac("S=C1CCC(CC1)[C@@H](C)CC") == "4-[(2S)-butan-2-yl]cyclohexane-1-thione"


def test_phenyl_chain_thione():
    assert smiles_to_iupac("c1ccccc1CC(=S)C") == "1-phenylpropane-2-thione"


def test_phenyl_directly_attached_thione_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=S)C")


def test_phenyl_chain_thione_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=S)C")


def test_thiourea():
    assert smiles_to_iupac("NC(=S)N") == "thiourea"


def test_n_methylthiourea():
    assert smiles_to_iupac("CNC(=S)N") == "N-methylthiourea"


def test_n_n_dimethylthiourea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=S)N") == "N,N-dimethylthiourea"


def test_n_ethyl_n_methylthiourea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=S)N") == "N-ethyl-N-methylthiourea"


def test_different_substituent_counts_on_different_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCN(C)C(=S)NC")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=S)N")


def test_n_methyl_n_prime_phenylthiourea_different_nitrogens():
    assert smiles_to_iupac("CNC(=S)Nc1ccccc1") == "N-methyl-N'-phenylthiourea"


def test_n_n_prime_diphenylthiourea():
    assert smiles_to_iupac("c1ccc(NC(=S)Nc2ccccc2)cc1") == "N,N'-diphenylthiourea"


def test_substituted_phenyl_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=S)Nc1ccc(C)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=S)N")
