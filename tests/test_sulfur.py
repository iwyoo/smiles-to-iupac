import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methyl_ethyl_disulfide():
    assert smiles_to_iupac("CSSCC") == "(methyldisulfanyl)ethane"



@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CSSSSC", "(methyltetrasulfanyl)methane", id="tetrasulfide_chain"),
        pytest.param("CSS", "methanedithioperoxol", id="terminal_persulfide_methane"),
        pytest.param("CCSS", "ethanedithioperoxol", id="terminal_persulfide_ethane"),
        pytest.param("c1ccccc1SS", "benzenedithioperoxol", id="terminal_persulfide_benzene"),
        pytest.param("CCCOS", "propane-1-OS-thioperoxol", id="oxygen_sulfur_peroxol_analogue"),
        pytest.param("CSO", "methane-SO-thioperoxol", id="sulfur_oxygen_peroxol_analogue"),
    ],
)
def test_tetrasulfide_chain_and_related(smiles, expected):
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


def test_salt_of_partial_ester_multivalent_cation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COS(=O)(=O)[O-].[Ca+2]")


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
    "smiles",
    [
        pytest.param("C=C(S(=O)N)C", id="ene_carbon_not_supported"),
        pytest.param("NS(=O)CS(=O)N", id="two_sulfinamides_not_supported"),
    ],
)
def test_ene_carbon_not_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
    "smiles",
    [
        pytest.param("O=S(N)C1CCCC=C1C", id="with_substituent_raises"),
        pytest.param("O=S(N)C1CCCC#C1", id="triple_bond_raises"),
        pytest.param("O=S(NC)C1CCCC=C1", id="with_n_substituent_raises"),
    ],
)
def test_unsaturated_ring_sulfinamide_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_ring_substituent_chain_sulfinamide():
    assert smiles_to_iupac("NS(=O)CC1CCCCC1") == "cyclohexylmethanesulfinamide"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NS(=O)CC1CCC(C)CC1", id="ring_substituent_chain_sulfinamide_ring_with_substituent_raises"),
        pytest.param("NS(=O)CCO", id="sulfinamide_with_alcohol_not_supported"),
    ],
)
def test_ring_substituent_chain_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1S(=O)N", "benzenesulfinamide", id="benzenesulfinamide"),
        pytest.param("Cc1ccc(cc1)S(=O)NC", "N,4-dimethylbenzenesulfinamide", id="substituted_benzenesulfinamide_n_alkyl"),
    ],
)
def test_benzenesulfinamide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_sulfinamide_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCCS(=O)NC")


def test_phenyl_chain_sulfinamide_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CS(=O)N)cc1") == "(4-chlorophenyl)methanesulfinamide"


def test_phenyl_chain_sulfinamide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCS(=O)N")


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
    "smiles",
    [
        pytest.param("CS(=O)NCCCl", id="halogenated_n_substituent_not_supported"),
        pytest.param("CS(=O)NCC=C", id="unsaturated_n_substituted_sulfinamide_not_supported"),
        pytest.param("CC[C@@H](C)S(=O)N", id="sulfinamide_specified_chain_stereocenter_raises"),
    ],
)
def test_halogenated_n_substituent_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


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
def test_ring_substituent_chain_and_related(smiles, expected):
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


def test_salt_of_partial_ester_multivalent_cation_raises__sulfite():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COS(=O)[O-].[Ca+2]")


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


def test_unsaturated_ring_sulfonamide_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CCCC#C1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NS(=O)(=O)CC1CCCCC1", "cyclohexylmethanesulfonamide", id="ring_substituent_chain_sulfonamide"),
        pytest.param("NS(=O)(=O)CC1CCC(C)CC1", "(4-methylcyclohexyl)methanesulfonamide", id="ring_substituent_chain_sulfonamide_ring_with_substituent"),
        pytest.param("c1ccccc1S(=O)(=O)N", "benzenesulfonamide", id="benzenesulfonamide"),
        pytest.param("Cc1ccc(cc1)S(=O)(=O)NC", "N,4-dimethylbenzenesulfonamide", id="substituted_benzenesulfonamide_n_alkyl"),
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
        pytest.param("C=CS(=O)(=O)[O-]", "eth-1-ene-1-sulfonate", id="sulfonate_carbon_in_double_bond_is_named"),
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


def test_unsaturated_ring_sulfonic_acid_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC#C1")


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
        pytest.param("CC(C)[SH2+]", "propan-2-ylsulfanium", id="branched_substituent"),
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


def test_sulfonium_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[S@+](CC)CCC")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[S@@+](CC)CCC")


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
        pytest.param("Clc1ccc(CC(=O)[S-])cc1", "2-(4-chlorophenyl)ethanethioate", id="chain_thioate_ring_halogen"),
        pytest.param("C=Cc1ccccc1CC(=O)[S-]", "(2-ethenylphenyl)ethanethioate", id="chain_thioate_unsaturation_is_named"),
    ],
)
def test_phenyl_directly_attached_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


def test_unsaturated_ring_thiol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC#C1")


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
        pytest.param("SC1CCCC(c2ccccc2)C1c3ccccc3", "2,3-diphenylcyclohexane-1-thiol", id="three_rings"),
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


def test_unsaturated_ring_thione_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CCCC#C1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=S)C1CCCCC1", "1-cyclohexylethane-1-thione", id="substituent_chain_thione"),
        pytest.param("S=C1CCCCC1CC(=S)C", "2-(2-sulfanylidenepropyl)cyclohexane-1-thione", id="with_thione_chain_thione_tie"),
    ],
)
def test_ring_substituent_chain_and_related_5(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_thial_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=S")


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
        pytest.param("c1ccccc1S(=O)(=O)NN", "benzenesulfonohydrazide", id="sulfonohydrazide_on_a_ring"),
        pytest.param("NNS(=O)(=O)c1ccccc1C(=O)O", "2-(hydrazinesulfonyl)benzoic acid", id="senior_acid_cites_the_prefix"),
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
        pytest.param("CN(C)S(=O)(=O)N", "N,N-dimethylsulfuric diamide", id="substituents_on_one_nitrogen"),
        pytest.param("CNS(=O)(=O)N(C)C", "N,N,N'-trimethylsulfuric diamide", id="substituents_on_both_nitrogens"),
        pytest.param("CNS(N)=O", "N-methylsulfurous diamide", id="sulfurous_diamide"),
    ],
)
def test_sulfamic_acid_and_the_amides_of_sulfuric_and_sulfurous_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "CS(=O)(=O)S(=O)(=O)c1ccccc1", "1-methyl-2-phenyl-1λ6,2λ6-disulfane-1,1,2,2-tetrone", id="disulfane_tetrone"
        ),
        pytest.param("CSCl", "methyl thiohypochlorite", id="thiohypochlorite_ester"),
    ],
)
def test_chalcogen_chain_heterones_and_halogen_acid_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[AsH5]", "λ5-arsane", id="pnictogen_hypervalent"),
        pytest.param("[IH3]", "λ3-iodane", id="halogen_hypervalent"),
                pytest.param("C[SH3]", "methyl-λ4-sulfane", id="substituted_hypervalent"),
        pytest.param("SSS", "trisulfane", id="homogeneous_chalcogen_chain"),
        pytest.param("S[SH2]S", "2λ4-trisulfane", id="chain_lambda_locant"),
        pytest.param("OO", "dioxidane", id="dioxidane"),
    ],
)
def test_nonstandard_bonding_number_hydrides_and_chalcogen_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OS(=O)(=O)OCCC(O)=O", "3-(sulfooxy)propanoic acid", id="sulfo"),
        pytest.param("COS(=O)OCCC(O)=O", "3-[(methoxysulfinyl)oxy]propanoic acid", id="alkoxysulfinyl"),
        pytest.param("ClS(=O)(=O)OCCC(O)=O", "3-[(chlorosulfonyl)oxy]propanoic acid", id="halosulfonyl"),
        pytest.param("NS(=O)(=O)OCCC(O)=O", "3-(sulfamoyloxy)propanoic acid", id="sulfamoyl"),
    ],
)
def test_sulfur_acid_groups_attached_through_oxygen_under_a_carboxylic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[SH2]CCSCCSCCSC", "2λ4,5,8,11-tetrathiadodecane"),
        ("CSCC[SH4]CC[SH2]CCSC", "2,5λ6,8λ4,11-tetrathiadodecane"),
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
