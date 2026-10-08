import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se][Se]C", "(methyldiselanyl)methane", id="dimethyl_diselenide"),
        pytest.param("CC[Se][Se]CC", "(ethyldiselanyl)ethane", id="diethyl_diselenide"),
        pytest.param("C[Se][Se]CC", "(methyldiselanyl)ethane", id="methyl_ethyl_diselenide"),
    ],
)
def test_dimethyl_diselenide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se][SeH]", "methanediselenoperoxol", id="methane"),
        pytest.param("c1ccccc1[Se][SeH]", "benzenediselenoperoxol", id="benzene"),
        pytest.param("CC[Se][SeH]", "ethanediselenoperoxol", id="ethane"),
    ],
)
def test_terminal_perselenol_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_parent_chain():
    assert smiles_to_iupac("CC[C@H](C)[Se][Se]CC") == "(2S)-2-(ethyldiselanyl)butane"



def test_phenyl_diselenide_direct_bond():
    assert smiles_to_iupac("c1ccccc1[Se][Se]C") == "(methyldiselanyl)benzene"
    assert smiles_to_iupac("c1ccccc1[Se][Se]CC") == "(ethyldiselanyl)benzene"


@pytest.mark.parametrize(
    "smiles",
    [
    ],
)
def test_phenyl_diselenide_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Te][Te]C", "(methylditellanyl)methane", id="dimethyl_ditelluride"),
        pytest.param("CC[Te][Te]CC", "(ethylditellanyl)ethane", id="diethyl_ditelluride"),
        pytest.param("C[Te][Te]CC", "(methylditellanyl)ethane", id="methyl_ethyl_ditelluride"),
    ],
)
def test_dimethyl_ditelluride_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC[Te][TeH]", "ethaneditelluroperoxol", id="ethane"),
        pytest.param("C[Te][TeH]", "methaneditelluroperoxol", id="methane"),
    ],
)
def test_terminal_pertellurol_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_parent_chain__ditelluride():
    assert smiles_to_iupac("CC[C@H](C)[Te][Te]CC") == "(2S)-2-(ethylditellanyl)butane"



def test_phenyl_ditelluride_direct_bond():
    assert smiles_to_iupac("c1ccccc1[Te][Te]C") == "(methylditellanyl)benzene"


@pytest.mark.parametrize(
    "smiles",
    [
    ],
)
def test_phenyl_ditelluride_chain_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_phenyl_isoselenocyanate_chain():
    assert smiles_to_iupac("c1ccccc1C[N]=C=[Se]") == "(isoselenocyanatomethyl)benzene"


def test_isothiocyanate_still_works():
    # Sanity check: the sulfur analogue must not be misrouted here.
    assert smiles_to_iupac("CCN=C=S") == "isothiocyanatoethane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[N]=C=[Te]", "isotellurocyanatoethane"),
    ],
)
def test_isotellurocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_isotellurocyanate_chain():
    assert smiles_to_iupac("c1ccccc1C[N]=C=[Te]") == "(isotellurocyanatomethyl)benzene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC[N]=C=[Se]", "isoselenocyanatoethane", id="isoselenocyanate_still_works"),
        pytest.param("C[Se]CC", "(methylselanyl)ethane", id="ethyl_methyl_selenide"),
        pytest.param("CCC[Se]C", "1-(methylselanyl)propane", id="methyl_propyl_selenide"),
    ],
)
def test_isoselenocyanate_still_works_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[Se]C(C)C", "2-[(propan-2-yl)selanyl]propane"),
        ("CC(C)C[Se]C(C)CC", "2-[(2-methylpropyl)selanyl]butane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1[Se]CC", "(ethylselanyl)benzene"),  # PubChem CID 140285
        ("c1ccccc1[Se]C(C)C", "[(propan-2-yl)selanyl]benzene"),
        ("c1ccccc1C[Se]CC", "[(ethylselanyl)methyl]benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se]c1ccccc1[Se]C", "1,2-bis(methylselanyl)benzene", id="benzene_ring_multiple_substituents_named_with_prefix"),
        pytest.param("c1ccccc1[Se][C@H](C)CC", "{[(2R)-butan-2-yl]selanyl}benzene", id="benzene_ring_stereocenter_named_with_prefix"),
        pytest.param("CC(Cl)[Se](=O)O", "1-chloroethane-1-seleninic acid", id="halogen_substituent"),
    ],
)
def test_benzene_ring_multiple_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=C([Se](=O)O)C", "prop-1-ene-2-seleninic acid", id="ene_carbon"),
        pytest.param("O[Se](=O)C[Se](=O)O", "methanediseleninic acid", id="two_seleninic_acids"),
        pytest.param("O[Se](=O)CCO", "2-hydroxyethane-1-seleninic acid", id="seleninic_acid_with_alcohol"),
    ],
)
def test_seleninic_acid_with_unsaturation_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_seleninic_acid():
    assert smiles_to_iupac("O=[Se](O)Cc1ccccc1") == "phenylmethaneseleninic acid"
    assert smiles_to_iupac("c1ccccc1CCC[Se](=O)O") == "3-phenylpropane-1-seleninic acid"


def test_benzeneseleninic_acid():
    assert smiles_to_iupac("c1ccccc1[Se](=O)O") == "benzeneseleninic acid"


def test_substituted_benzeneseleninic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Se](=O)O") == "2-methylbenzene-1-seleninic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Se](=O)O") == "4-methylbenzene-1-seleninic acid"
    )  # PubChem PUG REST


def test_phenyl_chain_seleninic_acid_with_an_unsaturated_ring_substituent():
    assert smiles_to_iupac("C=Cc1ccccc1CC[Se](=O)O") == "2-(2-ethenylphenyl)ethane-1-seleninic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("ClCC(=O)[Se-]", "2-chloroethaneselenoate", id="selenoate_with_halogen_substituent"),
        pytest.param("[Se-]C(=O)CCC(=O)[Se-]", "butanebis(selenoate)", id="multiple_selenoate_groups_is_named"),
        pytest.param("[Se-]C(=O)C1CCCCC1", "cyclohexanecarboselenoate", id="selenoate_on_ring_is_named"),
        pytest.param("CC(C)C(=O)[Se-]", "2-methylpropaneselenoate", id="selenoate_branched_r_group"),
    ],
)
def test_selenoate_with_halogen_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(=O)[Se-]", "(2S)-2-methylbutaneselenoate"),
    ],
)
def test_selenoate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_selenoate():
    assert smiles_to_iupac("c1ccccc1CC(=O)[Se-]") == "2-phenylethaneselenoate"
    assert smiles_to_iupac("c1ccccc1CCC(=O)[Se-]") == "3-phenylpropaneselenoate"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1C(=O)[Se-]", "benzenecarboselenoate", id="directly_attached_selenoate_is_named"),
        pytest.param("Cc1ccccc1CC(=O)[Se-]", "2-(2-methylphenyl)ethaneselenoate", id="substituted_benzene_ring_selenoate_is_named"),
        pytest.param("C=Cc1ccccc1CC(=O)[Se-]", "(2-ethenylphenyl)ethaneselenoate", id="chain_selenoate_unsaturation_is_named"),
    ],
)
def test_phenyl_directly_attached_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_selenocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#C[Se]CC[Se]C#N")


def test_phenyl_selenocyanate():
    assert smiles_to_iupac("c1ccccc1[Se]C#N") == "phenyl selenocyanate"  # CID 555340


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[Se])O", "ethaneselenoic O-acid"),
        ("O=C[SeH]", "methaneselenoic Se-acid"),
    ],
)
def test_selenoic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_substituent_selenoic_acid_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CC(=O)[SeH])cc1") == "2-(4-chlorophenyl)ethaneselenoic Se-acid"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("[SeH]C1CCCC#C1", id="triple_bond_raises"),
    ],
)
def test_unsaturated_ring_selenol_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C12(CCC(CC1)CC2)[SeH]", "bicyclo[2.2.2]octane-1-selenol", id="polycyclic_selenol_on_bridgehead"),
        pytest.param("C1CCCCC1C[SeH]", "cyclohexylmethaneselenol", id="ring_substituent_chain_selenol"),
        pytest.param("[SeH]C1CCCCC1CC[SeH]", "2-(2-selanylethyl)cyclohexane-1-selenol", id="ring_with_selenol_chain_selenol_tie"),
        pytest.param("C[C@@H](CC)C[SeH]", "(2S)-2-methylbutane-1-selenol", id="acyclic_selenol_stereocenter"),
        pytest.param("[SeH]C1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexane-1-selenol", id="cyclic_selenol_branch_stereocenter"),
        pytest.param("[SeH][C@H]1CCCCC1Cl", "(1S)-2-chlorocyclohexane-1-selenol", id="selenol_partially_specified_stereocenters_cites_the_specified_elements"),
        pytest.param("c1ccccc1[SeH]", "benzeneselenol", id="benzeneselenol"),
        pytest.param("Cc1ccccc1[SeH]", "2-methylbenzeneselenol", id="substituted_benzeneselenol"),
    ],
)
def test_polycyclic_selenol_on_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_oxygen_or_sulfur_heteroaromatic_chain_selenol():
    assert smiles_to_iupac("[SeH]Cc1cccs1") == "(thiophen-2-yl)methaneselenol"
    assert smiles_to_iupac("[SeH]Cc1ccco1") == "(furan-2-yl)methaneselenol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[Se](=O)(=O)CC", "(ethaneselenonyl)ethane"),
        ("C[Se](=O)(=O)CCCC", "1-(methaneselenonyl)butane"),
    ],
)
def test_selenone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported__selenone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se](=O)(=O)C")


def test_ring_selenone_is_a_lambda6_heterone():
    assert smiles_to_iupac("O=[Se]1(=O)CCCCC1") == "1λ6-selenane-1,1-dione"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C=C[Se](=O)(=O)C", id="unsaturated_chain_not_supported__selenone"),
        pytest.param("C[Se](=O)(=O)C[Se](=O)(=O)C", id="two_selenone_groups_not_supported"),
    ],
)
def test_unsaturated_chain_not_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se](=O)(=O)O", "methaneselenonic acid", id="methaneselenonic_acid"),
        pytest.param("C=CCCC[Se](=O)(=O)O", "pent-4-ene-1-selenonic acid", id="pent_4_ene_1_selenonic_acid"),
    ],
)
def test_methaneselenonic_acid_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O[Se](=O)(=O)C[Se](=O)(=O)O", "methanediselenonic acid", id="two_selenonic_acids"),
        pytest.param("O[Se](=O)(=O)CCO", "2-hydroxyethane-1-selenonic acid", id="selenonic_acid_with_alcohol"),
    ],
)
def test_selenonic_acid_with_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzeneselenonic_acid():
    assert smiles_to_iupac("c1ccccc1[Se](=O)(=O)O") == "benzeneselenonic acid"


def test_substituted_benzeneselenonic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Se](=O)(=O)O") == "2-methylbenzene-1-selenonic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Se](=O)(=O)O") == "4-methylbenzene-1-selenonic acid"
    )  # PubChem PUG REST


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C@@H](Cl)[Se](=O)(=O)O", "(1S)-1-chloroethane-1-selenonic acid", id="acyclic_selenonic_acid_stereocenter_with_coexisting_substituent"),
        pytest.param("NC(=[Se])N", "selenourea", id="selenourea"),
        pytest.param("CN(C)C(=[Se])N", "N,N-dimethylselenourea", id="n_n_dimethylselenourea_same_nitrogen"),
        pytest.param("CCN(C)C(=[Se])N", "N-ethyl-N-methylselenourea", id="n_ethyl_n_methylselenourea_same_nitrogen"),
    ],
)
def test_acyclic_selenonic_acid_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC(=[Se])Nc1ccccc1", "phenylselenourea", id="phenylselenourea"),
        pytest.param("c1ccc(NC(=[Se])Nc2ccccc2)cc1", "N,N'-diphenylselenourea", id="n_prime_diphenylselenourea"),
    ],
)
def test_n_phenylselenourea_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCC[Se](=O)CC", "1-(ethaneseleninyl)butane"),
        ("C[Se](C)=O", "(methaneseleninyl)methane"),
    ],
)
def test_selenoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported__selenoxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se](=O)C")


def test_ring_selenoxide_is_a_lambda4_heterone():
    assert smiles_to_iupac("O=[Se]1CCCCC1") == "1λ4-selenan-1-one"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C=C[Se](=O)C", id="unsaturated_chain_not_supported__selenoxide"),
        pytest.param("C[Se](=O)C[Se](=O)C", id="two_selenoxide_groups_not_supported"),
    ],
)
def test_unsaturated_chain_not_and_related_raise_2(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC(=[Se])CC1", "cyclohexaneselone"),
        ("ClCC(=[Se])C", "1-chloropropane-2-selone"),
    ],
)
def test_selone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_selone_with_substituent_is_named():
    assert smiles_to_iupac("[Se]=C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-selone"


def test_unsaturated_ring_selone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCCC#C1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=[Se])C1CCC(C)CC1", "1-(4-methylcyclohexyl)ethane-1-selone", id="substituent_chain_selone_ring_with_substituent"),
        pytest.param("[Se]=C1CCCCC1C(=[Se])C", "2-(ethaneselenoyl)cyclohexane-1-selone", id="with_selone_chain_selone_tie"),
    ],
)
def test_ring_substituent_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OCC(=[Se])C", "1-hydroxypropane-2-selone", id="selone_with_hydroxyl_is_named"),
        pytest.param("CC[C@@H](C)C(C)=[Se]", "(3R)-3-methylpentane-2-selone", id="acyclic_selone_stereocenter"),
        pytest.param("[Se]=C1CCC(CC1)[C@@H](C)CC", "4-[(2S)-butan-2-yl]cyclohexane-1-selone", id="cyclic_selone_branch_stereocenter"),
        pytest.param("c1ccccc1CC(=[Se])C", "1-phenylpropane-2-selone", id="phenyl_chain_selone"),
        pytest.param("c1ccccc1C(=[Se])C", "1-phenylethane-1-selone", id="phenyl_directly_attached_selone_is_named"),
        pytest.param("Cc1ccccc1CC(=[Se])C", "1-(2-methylphenyl)propane-2-selone", id="phenyl_substituted_benzene_ring_selone_is_named"),
        pytest.param("CC(=[Te])C", "propane-2-tellone", id="propane_2_tellone"),
    ],
)
def test_selone_with_hydroxyl_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_tellone():
    assert smiles_to_iupac("[Te]=C1CCCC=C1") == "cyclohex-2-ene-1-tellone"
    assert smiles_to_iupac("[Te]=C1CC=CCC1") == "cyclohex-3-ene-1-tellone"


def test_unsaturated_ring_tellone_with_substituent_is_named():
    assert smiles_to_iupac("[Te]=C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-tellone"


def test_unsaturated_ring_tellone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC#C1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=[Te])C1CCCCC1", "1-cyclohexylethane-1-tellone", id="substituent_chain_tellone"),
        pytest.param("[Te]=C1CCCCC1C(=[Te])C", "2-(ethanetelluroyl)cyclohexane-1-tellone", id="with_tellone_chain_tellone_tie"),
    ],
)
def test_ring_substituent_chain_and_related_2(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[Te]=C1CCC2(CCCCC2)CC1", "spiro[5.5]undecane-3-tellone", id="von_baeyer_spiro_tellone_now_supported"),
        pytest.param("CC[C@@H](C)C(C)=[Te]", "(3R)-3-methylpentane-2-tellone", id="acyclic_tellone_stereocenter"),
        pytest.param("[Te]=C1CCC(CC1)[C@@H](C)CC", "4-[(2S)-butan-2-yl]cyclohexane-1-tellone", id="cyclic_tellone_branch_stereocenter"),
        pytest.param("c1ccccc1CC(=[Te])C", "1-phenylpropane-2-tellone", id="phenyl_chain_tellone"),
        pytest.param("c1ccccc1C(=[Te])C", "1-phenylethane-1-tellone", id="phenyl_directly_attached_tellone_is_named"),
        pytest.param("Cc1ccccc1CC(=[Te])C", "1-(2-methylphenyl)propane-2-tellone", id="phenyl_substituted_benzene_ring_tellone_is_named"),
        pytest.param("C[Te]CC", "(methyltellanyl)ethane", id="methyl_ethyl_telluride"),
        pytest.param("CCC[Te]C", "1-(methyltellanyl)propane", id="methyl_propyl_telluride"),
    ],
)
def test_von_baeyer_spiro_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[Te]C(C)C", "2-[(propan-2-yl)tellanyl]propane"),
        ("CC(C)C[Te]C(C)CC", "2-[(2-methylpropyl)tellanyl]butane"),
    ],
)
def test_both_sides_branched_and_tied__telluride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1[Te]CC", "(ethyltellanyl)benzene"),  # PubChem CID 5325650
        ("c1ccccc1[Te]C(C)C", "[(propan-2-yl)tellanyl]benzene"),
        ("c1ccccc1C[Te]CC", "[(ethyltellanyl)methyl]benzene"),
    ],
)
def test_benzene_ring_parent__telluride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Te]c1ccccc1[Te]C", "1,2-bis(methyltellanyl)benzene", id="multiple_substituents_named_with_prefix__telluride"),
        pytest.param("c1ccccc1[Te][C@H](C)CC", "{[(2R)-butan-2-yl]tellanyl}benzene", id="stereocenter_named_with_prefix__telluride"),
    ],
)
def test_benzene_ring_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Te](=O)O", "methanetellurinic acid"),
    ],
)
def test_saturated_tellurinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituent__tellurinic_acid():
    assert smiles_to_iupac("CC(Cl)[Te](=O)O") == "1-chloroethane-1-tellurinic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=C([Te](=O)O)C", "prop-1-ene-2-tellurinic acid", id="ene_carbon"),
        pytest.param("O[Te](=O)C[Te](=O)O", "methaneditellurinic acid", id="two_tellurinic_acids"),
        pytest.param("O[Te](=O)CCO", "2-hydroxyethane-1-tellurinic acid", id="tellurinic_acid_with_alcohol"),
    ],
)
def test_tellurinic_acid_with_unsaturation_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzenetellurinic_acid():
    assert smiles_to_iupac("c1ccccc1[Te](=O)O") == "benzenetellurinic acid"


def test_substituted_benzenetellurinic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Te](=O)O") == "2-methylbenzene-1-tellurinic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Te](=O)O") == "4-methylbenzene-1-tellurinic acid"
    )  # PubChem PUG REST


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)[Te](=O)O", "(2S)-butane-2-tellurinic acid"),
    ],
)
def test_tellurinic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_tellurocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#C[Te]CC[Te]C#N")


def test_phenyl_tellurocyanate():
    assert smiles_to_iupac("c1ccccc1[Te]C#N") == "phenyl tellurocyanate"  # CID 12553982


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C[TeH]", "methanetelluroic Te-acid"),
        ("CC(=[Te])O", "ethanetelluroic O-acid"),
    ],
)
def test_telluroic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_substituent_telluroic_acid_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CC(=O)[TeH])cc1") == "2-(4-chlorophenyl)ethanetelluroic Te-acid"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("[TeH]C1CCCC#C1", id="triple_bond_raises"),
    ],
)
def test_unsaturated_ring_tellurol_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1=CCCCC1C[TeH]", "(cyclohex-2-en-1-yl)methanetellurol", id="ring_substituent_chain_tellurol_unsaturated_ring"),
        pytest.param("[TeH]C1CCCCC1CC[TeH]", "2-(2-tellanylethyl)cyclohexane-1-tellurol", id="ring_with_tellurol_chain_tellurol_tie"),
        pytest.param("OCC[TeH]", "2-tellanylethan-1-ol", id="tellurol_with_alcohol_named_with_prefix"),
        pytest.param("C[C@@H](CC)C[TeH]", "(2S)-2-methylbutane-1-tellurol", id="acyclic_tellurol_stereocenter"),
        pytest.param("[TeH]C1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexane-1-tellurol", id="cyclic_tellurol_branch_stereocenter"),
        pytest.param("[TeH][C@H]1CCCCC1Cl", "(1S)-2-chlorocyclohexane-1-tellurol", id="tellurol_partially_specified_stereocenters_cites_the_specified_elements"),
        pytest.param("c1ccccc1[TeH]", "benzenetellurol", id="benzenetellurol"),
        pytest.param("Cc1ccccc1[TeH]", "2-methylbenzenetellurol", id="substituted_benzenetellurol"),
    ],
)
def test_ring_substituent_chain_and_related_3(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[TeH]Cc1cccnc1", "(pyridin-3-yl)methanetellurol"),
    ],
)
def test_heteroaromatic_chain_tellurol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[Te](=O)(=O)CC", "(ethanetelluronyl)ethane"),
        ("C[Te](=O)(=O)CCCC", "1-(methanetelluronyl)butane"),
    ],
)
def test_tellurone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported__tellurone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te](=O)(=O)C")


def test_ring_tellurone_is_a_lambda6_heterone():
    assert smiles_to_iupac("O=[Te]1(=O)CCCCC1") == "1λ6-tellurane-1,1-dione"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C=C[Te](=O)(=O)C", id="unsaturated_chain_not_supported__tellurone"),
        pytest.param("C[Te](=O)(=O)C[Te](=O)(=O)C", id="two_tellurone_groups_not_supported"),
    ],
)
def test_unsaturated_chain_not_and_related_raise_3(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Te](=O)(=O)O", "methanetelluronic acid", id="methanetelluronic_acid"),
        pytest.param("C=CCCC[Te](=O)(=O)O", "pent-4-ene-1-telluronic acid", id="pent_4_ene_1_telluronic_acid"),
    ],
)
def test_methanetelluronic_acid_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O[Te](=O)(=O)C[Te](=O)(=O)O", "methaneditelluronic acid", id="two_telluronic_acids"),
        pytest.param("O[Te](=O)(=O)CCO", "2-hydroxyethane-1-telluronic acid", id="telluronic_acid_with_alcohol"),
    ],
)
def test_telluronic_acid_with_other_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzenetelluronic_acid():
    assert smiles_to_iupac("c1ccccc1[Te](=O)(=O)O") == "benzenetelluronic acid"


def test_substituted_benzenetelluronic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Te](=O)(=O)O") == "2-methylbenzene-1-telluronic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Te](=O)(=O)O") == "4-methylbenzene-1-telluronic acid"
    )  # PubChem PUG REST


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C@@H](Cl)[Te](=O)(=O)O", "(1S)-1-chloroethane-1-telluronic acid", id="acyclic_telluronic_acid_stereocenter_with_coexisting_substituent"),
        pytest.param("NC(=[Te])N", "tellurourea", id="tellurourea"),
        pytest.param("CN(C)C(=[Te])N", "N,N-dimethyltellurourea", id="n_n_dimethyltellurourea_same_nitrogen"),
        pytest.param("CCN(C)C(=[Te])N", "N-ethyl-N-methyltellurourea", id="n_ethyl_n_methyltellurourea_same_nitrogen"),
    ],
)
def test_acyclic_telluronic_acid_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC(=[Te])Nc1ccccc1", "phenyltellurourea", id="phenyltellurourea"),
        pytest.param("c1ccc(NC(=[Te])Nc2ccccc2)cc1", "N,N'-diphenyltellurourea", id="n_prime_diphenyltellurourea"),
    ],
)
def test_n_phenyltellurourea_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Te](C)=O", "(methanetellurinyl)methane"),
        ("CCCC[Te](=O)CC", "1-(ethanetellurinyl)butane"),
    ],
)
def test_telluroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported__telluroxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te](=O)C")


def test_ring_telluroxide_is_a_lambda4_heterone():
    assert smiles_to_iupac("O=[Te]1CCCCC1") == "1λ4-telluran-1-one"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C=C[Te](=O)C", id="unsaturated_chain_not_supported__telluroxide"),
        pytest.param("C[Te](=O)C[Te](=O)C", id="two_telluroxide_groups_not_supported"),
    ],
)
def test_unsaturated_chain_not_and_related_raise_4(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC[C@H](C)[Se]C", "(2S)-2-(methylselanyl)butane"),
        ("CC[C@@H](C)[Se]C", "(2R)-2-(methylselanyl)butane"),
        ("CC[C@H](C)[Te]C", "(2S)-2-(methyltellanyl)butane"),
        ("C[C@H]([Se]C)[C@H](C)CC", "(2S,3R)-3-methyl-2-(methylselanyl)pentane"),
        ("C[C@H]([Te]C)C(C)(C)C", "(3S)-2,2-dimethyl-3-(methyltellanyl)butane"),
    ],
)
def test_selenide_telluride_stereodescriptors_like_sulfide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CCNC(=[Se])NC", "N-ethyl-N'-methylselenourea"),
        ("c1ccccc1NC(=[Se])NC", "N-methyl-N'-phenylselenourea"),
        ("CCNC(=[Te])NC", "N-ethyl-N'-methyltellurourea"),
        ("CCNC(=[Te])NCCC", "N-ethyl-N'-propyltellurourea"),
    ],
)
def test_selenourea_tellurourea_different_n_substituents_like_thiourea(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("C[Se][Se][Se]C", "(methyltriselanyl)methane"),
        ("CC[Se][Se][Se]CC", "(ethyltriselanyl)ethane"),
        ("C[Se][Se][Se][Se]C", "(methyltetraselanyl)methane"),
        ("CC[Se][Se][Se][SeH]", "tetraselanylethane"),
        ("c1ccccc1[Se][Se][Se]C", "(methyltriselanyl)benzene"),
        ("CC[C@H](C)[Se][Se][Se]C", "(2S)-2-(methyltriselanyl)butane"),
        ("CC[Te][Te][Te]CC", "(ethyltritellanyl)ethane"),
        ("C[Te][Te][Te][Te]C", "(methyltetratellanyl)methane"),
    ],
)
def test_polyselenides_and_polytellurides_like_polysulfides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [("c1ccc[se]1", "selenophene"), ("Cc1ccc(C)[se]1", "2,5-dimethylselenophene")],
)
def test_aromatic_selenium_ring_is_not_a_selenide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("O[Se](=O)(=O)C1CCCCC1", "cyclohexaneselenonic acid"),
        ("O[Te](=O)(=O)C1CCCCC1", "cyclohexanetelluronic acid"),
        ("C1CC2CC1CC2[Se](=O)(=O)O", "bicyclo[2.2.1]heptane-2-selenonic acid"),
        ("C[C@H]1CCCC[C@@H]1[Se](=O)(=O)O", "(1S,2S)-2-methylcyclohexane-1-selenonic acid"),
        ("O[Se](=O)C1CCCCC1", "cyclohexaneseleninic acid"),
        ("O[Te](=O)C1CCCCC1", "cyclohexanetellurinic acid"),
        ("C1CCCCC1[Se@](=O)O", "(R)-cyclohexaneseleninic acid"),
        ("CC[C@@H](C)[Se](=O)O", "(2R)-butane-2-seleninic acid"),
        ("Clc1ccccc1C[Se](=O)(=O)O", "(2-chlorophenyl)methaneselenonic acid"),
    ],
)
def test_selenium_tellurium_oxo_acids_like_sulfur_analogues(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("Cc1ccccc1CC[SeH]", "2-(2-methylphenyl)ethane-1-selenol"),
        ("Cc1ccccc1CC[TeH]", "2-(2-methylphenyl)ethane-1-tellurol"),
        ("Clc1ccc(cc1)CC[SeH]", "2-(4-chlorophenyl)ethane-1-selenol"),
        ("c1ccccc1C([SeH])CC[SeH]", "1-phenylpropane-1,3-diselenol"),
        ("c1ccccc1C([TeH])CC[TeH]", "1-phenylpropane-1,3-ditellurol"),
        ("c1ccccc1CC([SeH])Cc1ccccc1", "1,3-diphenylpropane-2-selenol"),
        ("Cc1ccccc1CC[Se](=O)O", "2-(2-methylphenyl)ethane-1-seleninic acid"),
        ("Cc1ccccc1CC(=O)[Se-]", "2-(2-methylphenyl)ethaneselenoate"),
    ],
)
def test_selenium_tellurium_substituted_benzene_chains_like_sulfur_analogues(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("[SeH]CCN", "2-aminoethane-1-selenol"),
        ("[TeH]CCN", "2-aminoethane-1-tellurol"),
        ("[SeH]c1cccnc1", "pyridine-3-selenol"),
        ("[SeH]C1CC2CCC1C([SeH])C2", "bicyclo[2.2.2]octane-2,6-diselenol"),
        ("[SeH]C1CCCC=C1C", "2-methylcyclohex-2-ene-1-selenol"),
    ],
)
def test_selenols_and_tellurols_beside_other_groups_and_on_any_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se](=O)(=O)NN", "methaneselenonohydrazide", id="selenonohydrazide"),
        pytest.param("C[Se](=O)NN", "methaneseleninohydrazide", id="seleninohydrazide"),
        pytest.param("C[Te](=O)(=O)NN", "methanetelluronohydrazide", id="telluronohydrazide"),
        pytest.param("c1ccccc1[Te](=O)NN", "benzenetellurinohydrazide", id="tellurinohydrazide_on_a_ring"),
    ],
)
def test_selenium_and_tellurium_hydrazides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[Se](=N)(=O)N", "methaneselenonimidamide", id="selenonimidamide"),
        pytest.param("C[Se](=N)N", "methaneseleninimidamide", id="seleninimidamide"),
        pytest.param("C[Te](=N)(=N)N", "methanetelluronodiimidamide", id="telluronodiimidamide"),
    ],
)
def test_selenium_and_tellurium_imidamides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=[Se](N)c1ccco1", "furan-2-seleninamide", id="seleninamide_on_a_ring"),
        pytest.param("O=[Te](=O)(NC)c1ccccc1", "N-methylbenzenetelluronamide", id="telluronamide_n_substituted"),
        pytest.param("NC(=O)C[Se](=O)(=O)N", "2-(aminoselenonyl)acetamide", id="selenonamide_demoted_to_a_prefix"),
        pytest.param("NC(=O)CS(=O)NC", "2-[(methylamino)sulfinyl]acetamide", id="sulfinamide_demoted_to_a_prefix"),
    ],
)
def test_amides_of_sulfinic_selenium_and_tellurium_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
