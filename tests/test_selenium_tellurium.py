import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_diselenide():
    assert smiles_to_iupac("C[Se][Se]C") == "(methyldiselanyl)methane"


def test_diethyl_diselenide():
    assert smiles_to_iupac("CC[Se][Se]CC") == "(ethyldiselanyl)ethane"


def test_methyl_ethyl_diselenide():
    assert smiles_to_iupac("C[Se][Se]CC") == "(methyldiselanyl)ethane"


def test_branched_diselanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se][Se]C(C)C")


def test_terminal_perselenol_methane():
    assert smiles_to_iupac("C[Se][SeH]") == "diselanylmethane"


def test_terminal_perselenol_ethane():
    assert smiles_to_iupac("CC[Se][SeH]") == "diselanylethane"


def test_both_terminal_diselane_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH][SeH]")


def test_stereocenter_on_parent_chain():
    assert smiles_to_iupac("CC[C@H](C)[Se][Se]CC") == "(2S)-2-(ethyldiselanyl)butane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC[Se][Se][C@H](C)CC")


def test_phenyl_diselenide_direct_bond():
    assert smiles_to_iupac("c1ccccc1[Se][Se]C") == "(methyldiselanyl)benzene"
    assert smiles_to_iupac("c1ccccc1[Se][Se]CC") == "(ethyldiselanyl)benzene"


def test_phenyl_diselenide_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Se][Se]C")


def test_phenyl_diselenide_seh_terminal_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Se][SeH]")


def test_phenyl_diselenide_branched_other_side_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Se][Se]C(C)C")


def test_phenyl_diselenide_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[Se][Se]C")


def test_dimethyl_ditelluride():
    assert smiles_to_iupac("C[Te][Te]C") == "(methylditellanyl)methane"


def test_diethyl_ditelluride():
    assert smiles_to_iupac("CC[Te][Te]CC") == "(ethylditellanyl)ethane"


def test_methyl_ethyl_ditelluride():
    assert smiles_to_iupac("C[Te][Te]CC") == "(methylditellanyl)ethane"


def test_branched_ditellanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te][Te]C(C)C")


def test_terminal_pertellurol_ethane():
    assert smiles_to_iupac("CC[Te][TeH]") == "ditellanylethane"


def test_terminal_pertellurol_methane():
    assert smiles_to_iupac("C[Te][TeH]") == "ditellanylmethane"


def test_both_terminal_ditellane_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH][TeH]")


def test_stereocenter_on_parent_chain__ditelluride():
    assert smiles_to_iupac("CC[C@H](C)[Te][Te]CC") == "(2S)-2-(ethylditellanyl)butane"


def test_stereocenter_on_substituent_branch_not_supported__ditelluride():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC[Te][Te][C@H](C)CC")


def test_phenyl_ditelluride_direct_bond():
    assert smiles_to_iupac("c1ccccc1[Te][Te]C") == "(methylditellanyl)benzene"


def test_phenyl_ditelluride_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Te][Te]C")


def test_phenyl_ditelluride_teh_terminal_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Te][TeH]")


def test_phenyl_ditelluride_branched_other_side_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Te][Te]C(C)C")


def test_phenyl_ditelluride_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[Te][Te]C")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N]=C=[Se])CC1")


def test_phenyl_isoselenocyanate_chain():
    assert smiles_to_iupac("c1ccccc1C[N]=C=[Se]") == "(isoselenocyanatomethyl)benzene"


def test_phenyl_isoselenocyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[N]=C=[Se]")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N]=C=[Se]")


def test_two_isoselenocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C=NCN=C=[Se]")


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


def test_ring_not_supported__isotellurocyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N]=C=[Te])CC1")


def test_phenyl_isotellurocyanate_chain():
    assert smiles_to_iupac("c1ccccc1C[N]=C=[Te]") == "(isotellurocyanatomethyl)benzene"


def test_phenyl_isotellurocyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[N]=C=[Te]")


def test_unsaturated_chain_not_supported__isotellurocyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N]=C=[Te]")


def test_two_isotellurocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C=NCN=C=[Te]")


def test_isoselenocyanate_still_works():
    # Sanity check: the selenium analogue must not be misrouted here.
    assert smiles_to_iupac("CC[N]=C=[Se]") == "isoselenocyanatoethane"


def test_ethyl_methyl_selenide():
    assert smiles_to_iupac("C[Se]CC") == "methylselanylethane"


def test_methyl_propyl_selenide():
    assert smiles_to_iupac("CCC[Se]C") == "1-methylselanylpropane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[Se]C(C)C", "2-(propan-2-yl)selanylpropane"),
        ("CC(C)C[Se]C(C)CC", "2-(2-methylpropyl)selanylbutane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1[Se]CC", "ethylselanylbenzene"),  # PubChem CID 140285
        ("c1ccccc1[Se]C(C)C", "(propan-2-yl)selanylbenzene"),
        ("c1ccccc1C[Se]CC", "(ethylselanylmethyl)benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents_named_with_prefix():
    assert smiles_to_iupac("C[Se]c1ccccc1[Se]C") == "1,2-bis(methylselanyl)benzene"


def test_benzene_ring_stereocenter_named_with_prefix():
    assert smiles_to_iupac("c1ccccc1[Se][C@H](C)CC") == "{[(2R)-butan-2-yl]selanyl}benzene"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)[Se](=O)O") == "1-chloroethaneseleninic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C([Se](=O)O)C")


def test_two_seleninic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)C[Se](=O)O")


def test_ring_seleninic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)C1CCCCC1")


def test_seleninic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)CCO")


def test_seleninic_acid_specified_chain_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)[Se](=O)O")


def test_phenyl_chain_seleninic_acid():
    assert smiles_to_iupac("O=[Se](O)Cc1ccccc1") == "phenylmethaneseleninic acid"
    assert smiles_to_iupac("c1ccccc1CCC[Se](=O)O") == "3-phenylpropane-1-seleninic acid"


def test_benzeneseleninic_acid():
    assert smiles_to_iupac("c1ccccc1[Se](=O)O") == "benzeneseleninic acid"


def test_substituted_benzeneseleninic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Se](=O)O") == "2-methylbenzeneseleninic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Se](=O)O") == "4-methylbenzeneseleninic acid"
    )  # PubChem PUG REST


def test_phenyl_substituted_benzene_ring_seleninic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[Se](=O)O")


def test_phenyl_chain_seleninic_acid_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[Se](=O)O")


def test_selenoate_with_halogen_substituent():
    assert smiles_to_iupac("ClCC(=O)[Se-]") == "2-chloroethaneselenoate"


def test_multiple_selenoate_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se-]C(=O)CCC(=O)[Se-]")


def test_selenoate_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se-]C(=O)C1CCCCC1")


def test_selenoate_branched_r_group():
    assert smiles_to_iupac("CC(C)C(=O)[Se-]") == "2-methylpropaneselenoate"


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


def test_phenyl_directly_attached_selenoate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[Se-]")


def test_phenyl_substituted_benzene_ring_selenoate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)[Se-]")


def test_phenyl_chain_selenoate_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)[Se-]")


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Se]C#N")


def test_cyclic_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1[Se]C#N")


def test_two_selenocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#C[Se]CC[Se]C#N")


def test_phenyl_selenocyanate():
    assert smiles_to_iupac("c1ccccc1[Se]C#N") == "phenyl selenocyanate"  # CID 555340


def test_phenyl_selenocyanate_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Se]C#N")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[Se])O", "ethaneselenoic O-acid"),
        ("O=C[SeH]", "methaneselenoic Se-acid"),
    ],
)
def test_selenoic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)[SeH]")


def test_unsaturated_chain_not_supported__selenoic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)[SeH]")


def test_two_selenoic_acid_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C(=O)CC(=O)[SeH]")


def test_phenyl_directly_attached_selenoic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[SeH]")


def test_phenyl_substituent_selenoic_acid_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CC(=O)[SeH])cc1") == "2-(4-chlorophenyl)ethaneselenoic Se-acid"


def test_unsaturated_ring_selenol_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CCCC=C1C")


def test_unsaturated_ring_selenol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CCCC#C1")


def test_polycyclic_selenol_on_bridgehead():
    assert smiles_to_iupac("C12(CCC(CC1)CC2)[SeH]") == "bicyclo[2.2.2]octane-1-selenol"


def test_ring_substituent_chain_selenol():
    assert smiles_to_iupac("C1CCCCC1C[SeH]") == "cyclohexylmethaneselenol"


def test_ring_with_selenol_chain_selenol_tie():
    assert smiles_to_iupac("[SeH]C1CCCCC1CC[SeH]") == "2-(2-selanylethyl)cyclohexane-1-selenol"


def test_acyclic_selenol_stereocenter():
    assert smiles_to_iupac("C[C@@H](CC)C[SeH]") == "(2S)-2-methylbutane-1-selenol"


def test_cyclic_selenol_branch_stereocenter():
    assert (
        smiles_to_iupac("[SeH]C1(CCCCC1)[C@@H](C)CC") == "1-[(2S)-butan-2-yl]cyclohexane-1-selenol"
    )


def test_selenol_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("[SeH][C@H]1CCCCC1Cl") == "(1S)-2-chlorocyclohexane-1-selenol"


def test_benzeneselenol():
    assert smiles_to_iupac("c1ccccc1[SeH]") == "benzeneselenol"  # CID 69530


def test_substituted_benzeneselenol():
    assert smiles_to_iupac("Cc1ccccc1[SeH]") == "2-methylbenzeneselenol"


def test_phenyl_substituted_benzene_ring_selenol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[SeH]")


def test_phenyl_chain_selenol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[SeH]")


def test_phenyl_chain_diselenol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C([SeH])CC[SeH]")


def test_heteroaromatic_direct_attachment_selenol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]c1cccnc1")


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


def test_ring_not_supported__selenone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Se]1(=O)CCCCC1")


def test_unsaturated_chain_not_supported__selenone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Se](=O)(=O)C")


def test_two_selenone_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se](=O)(=O)C[Se](=O)(=O)C")


def test_methaneselenonic_acid():
    assert smiles_to_iupac("C[Se](=O)(=O)O") == "methaneselenonic acid"


def test_pent_4_ene_1_selenonic_acid():
    assert smiles_to_iupac("C=CCCC[Se](=O)(=O)O") == "pent-4-ene-1-selenonic acid"


def test_diselenonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)(=O)C[Se](=O)(=O)O")


def test_ring_selenonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)(=O)C1CCCCC1")


def test_selenonic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)(=O)CCO")


def test_benzeneselenonic_acid():
    assert smiles_to_iupac("c1ccccc1[Se](=O)(=O)O") == "benzeneselenonic acid"


def test_substituted_benzeneselenonic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Se](=O)(=O)O") == "2-methylbenzeneselenonic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Se](=O)(=O)O") == "4-methylbenzeneselenonic acid"
    )  # PubChem PUG REST


def test_acyclic_selenonic_acid_stereocenter_with_coexisting_substituent():
    assert smiles_to_iupac("C[C@@H](Cl)[Se](=O)(=O)O") == "(1S)-1-chloroethaneselenonic acid"


def test_selenourea():
    assert smiles_to_iupac("NC(=[Se])N") == "selenourea"


def test_n_n_dimethylselenourea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=[Se])N") == "N,N-dimethylselenourea"


def test_n_ethyl_n_methylselenourea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=[Se])N") == "N-ethyl-N-methylselenourea"


def test_different_substituents_on_different_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNC(=[Se])NC")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=[Se])N")


def test_n_phenylselenourea():
    assert smiles_to_iupac("NC(=[Se])Nc1ccccc1") == "N-phenylselenourea"


def test_n_n_prime_diphenylselenourea():
    assert smiles_to_iupac("c1ccc(NC(=[Se])Nc2ccccc2)cc1") == "N,N'-diphenylselenourea"


def test_substituted_phenyl_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=[Se])Nc1ccc(C)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=[Se])N")


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


def test_ring_not_supported__selenoxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Se]1CCCCC1")


def test_unsaturated_chain_not_supported__selenoxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Se](=O)C")


def test_two_selenoxide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se](=O)C[Se](=O)C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC(=[Se])CC1", "cyclohexaneselone"),
        ("ClCC(=[Se])C", "1-chloropropane-2-selone"),
    ],
)
def test_selone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_selone_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCCC=C1C")


def test_unsaturated_ring_selone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCCC#C1")


def test_ring_substituent_chain_selone_ring_with_substituent():
    assert smiles_to_iupac("CC(=[Se])C1CCC(C)CC1") == "1-(4-methylcyclohexyl)ethaneselone"


def test_ring_with_selone_chain_selone_tie():
    assert (
        smiles_to_iupac("[Se]=C1CCCCC1C(=[Se])C") == "2-(1-selanylideneethyl)cyclohexane-1-selone"
    )


def test_selenoaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Se]")


def test_selone_with_hydroxyl_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=[Se])C")


def test_acyclic_selone_stereocenter():
    assert smiles_to_iupac("CC[C@@H](C)C(C)=[Se]") == "(3R)-3-methylpentane-2-selone"


def test_cyclic_selone_branch_stereocenter():
    assert (
        smiles_to_iupac("[Se]=C1CCC(CC1)[C@@H](C)CC") == "4-[(2S)-butan-2-yl]cyclohexane-1-selone"
    )


def test_phenyl_chain_selone():
    assert smiles_to_iupac("c1ccccc1CC(=[Se])C") == "1-phenylpropane-2-selone"


def test_phenyl_directly_attached_selone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=[Se])C")


def test_phenyl_substituted_benzene_ring_selone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=[Se])C")


def test_propane_2_tellone():
    assert smiles_to_iupac("CC(=[Te])C") == "propane-2-tellone"


def test_unsaturated_ring_tellone():
    assert smiles_to_iupac("[Te]=C1CCCC=C1") == "cyclohex-2-ene-1-tellone"
    assert smiles_to_iupac("[Te]=C1CC=CCC1") == "cyclohex-3-ene-1-tellone"


def test_unsaturated_ring_tellone_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC=C1C")


def test_unsaturated_ring_tellone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC#C1")


def test_ring_substituent_chain_tellone():
    assert smiles_to_iupac("CC(=[Te])C1CCCCC1") == "1-cyclohexylethanetellone"


def test_ring_with_tellone_chain_tellone_tie():
    assert (
        smiles_to_iupac("[Te]=C1CCCCC1C(=[Te])C") == "2-(1-tellanylideneethyl)cyclohexane-1-tellone"
    )


def test_telluroaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Te]")


def test_von_baeyer_spiro_tellone_now_supported():
    assert smiles_to_iupac("[Te]=C1CCC2(CCCCC2)CC1") == "spiro[5.5]undecane-3-tellone"


def test_acyclic_tellone_stereocenter():
    assert smiles_to_iupac("CC[C@@H](C)C(C)=[Te]") == "(3R)-3-methylpentane-2-tellone"


def test_cyclic_tellone_branch_stereocenter():
    assert (
        smiles_to_iupac("[Te]=C1CCC(CC1)[C@@H](C)CC") == "4-[(2S)-butan-2-yl]cyclohexane-1-tellone"
    )


def test_phenyl_chain_tellone():
    assert smiles_to_iupac("c1ccccc1CC(=[Te])C") == "1-phenylpropane-2-tellone"


def test_phenyl_directly_attached_tellone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=[Te])C")


def test_phenyl_substituted_benzene_ring_tellone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=[Te])C")


def test_methyl_ethyl_telluride():
    assert smiles_to_iupac("C[Te]CC") == "methyltellanylethane"


def test_methyl_propyl_telluride():
    assert smiles_to_iupac("CCC[Te]C") == "1-methyltellanylpropane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)[Te]C(C)C", "2-(propan-2-yl)tellanylpropane"),
        ("CC(C)C[Te]C(C)CC", "2-(2-methylpropyl)tellanylbutane"),
    ],
)
def test_both_sides_branched_and_tied__telluride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1[Te]CC", "ethyltellanylbenzene"),  # PubChem CID 5325650
        ("c1ccccc1[Te]C(C)C", "(propan-2-yl)tellanylbenzene"),
        ("c1ccccc1C[Te]CC", "(ethyltellanylmethyl)benzene"),
    ],
)
def test_benzene_ring_parent__telluride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents_named_with_prefix__telluride():
    assert smiles_to_iupac("C[Te]c1ccccc1[Te]C") == "1,2-bis(methyltellanyl)benzene"


def test_benzene_ring_stereocenter_named_with_prefix__telluride():
    assert smiles_to_iupac("c1ccccc1[Te][C@H](C)CC") == "{[(2R)-butan-2-yl]tellanyl}benzene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Te](=O)O", "methanetellurinic acid"),
    ],
)
def test_saturated_tellurinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituent__tellurinic_acid():
    assert smiles_to_iupac("CC(Cl)[Te](=O)O") == "1-chloroethanetellurinic acid"


def test_ene_carbon_not_supported__tellurinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C([Te](=O)O)C")


def test_two_tellurinic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)C[Te](=O)O")


def test_ring_tellurinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)C1CCCCC1")


def test_tellurinic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)CCO")


def test_benzenetellurinic_acid():
    assert smiles_to_iupac("c1ccccc1[Te](=O)O") == "benzenetellurinic acid"


def test_substituted_benzenetellurinic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Te](=O)O") == "2-methylbenzenetellurinic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Te](=O)O") == "4-methylbenzenetellurinic acid"
    )  # PubChem PUG REST


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)[Te](=O)O", "(2S)-butane-2-tellurinic acid"),
    ],
)
def test_tellurinic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_r_not_supported__tellurocyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Te]C#N")


def test_cyclic_r_not_supported__tellurocyanate():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1[Te]C#N")


def test_two_tellurocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#C[Te]CC[Te]C#N")


def test_phenyl_tellurocyanate():
    assert smiles_to_iupac("c1ccccc1[Te]C#N") == "phenyl tellurocyanate"  # CID 12553982


def test_phenyl_tellurocyanate_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Te]C#N")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C[TeH]", "methanetelluroic Te-acid"),
        ("CC(=[Te])O", "ethanetelluroic O-acid"),
    ],
)
def test_telluroic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_substituent_not_supported__telluroic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)[TeH]")


def test_unsaturated_chain_not_supported__telluroic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)[TeH]")


def test_two_telluroic_acid_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C(=O)CC(=O)[TeH]")


def test_phenyl_directly_attached_telluroic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[TeH]")


def test_phenyl_substituent_telluroic_acid_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CC(=O)[TeH])cc1") == "2-(4-chlorophenyl)ethanetelluroic Te-acid"


def test_unsaturated_ring_tellurol_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CCCC=C1C")


def test_unsaturated_ring_tellurol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CCCC#C1")


def test_ring_substituent_chain_tellurol_unsaturated_ring():
    assert smiles_to_iupac("C1=CCCCC1C[TeH]") == "(cyclohex-2-en-1-yl)methanetellurol"


def test_ring_with_tellurol_chain_tellurol_tie():
    assert smiles_to_iupac("[TeH]C1CCCCC1CC[TeH]") == "2-(2-tellanylethyl)cyclohexane-1-tellurol"


def test_tellurol_with_alcohol_named_with_prefix():
    assert smiles_to_iupac("OCC[TeH]") == "2-tellanylethanol"


def test_acyclic_tellurol_stereocenter():
    assert smiles_to_iupac("C[C@@H](CC)C[TeH]") == "(2S)-2-methylbutane-1-tellurol"


def test_cyclic_tellurol_branch_stereocenter():
    assert (
        smiles_to_iupac("[TeH]C1(CCCCC1)[C@@H](C)CC") == "1-[(2S)-butan-2-yl]cyclohexane-1-tellurol"
    )


def test_tellurol_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("[TeH][C@H]1CCCCC1Cl") == "(1S)-2-chlorocyclohexane-1-tellurol"


def test_benzenetellurol():
    assert smiles_to_iupac("c1ccccc1[TeH]") == "benzenetellurol"  # CID 5246059


def test_substituted_benzenetellurol():
    assert smiles_to_iupac("Cc1ccccc1[TeH]") == "2-methylbenzenetellurol"


def test_phenyl_substituted_benzene_ring_tellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[TeH]")


def test_phenyl_chain_tellurol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[TeH]")


def test_phenyl_chain_ditellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C([TeH])CC[TeH]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[TeH]Cc1cccnc1", "(pyridin-3-yl)methanetellurol"),
    ],
)
def test_heteroaromatic_chain_tellurol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_heteroaromatic_direct_attachment_tellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]c1cccnc1")


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


def test_ring_not_supported__tellurone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Te]1(=O)CCCCC1")


def test_unsaturated_chain_not_supported__tellurone():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Te](=O)(=O)C")


def test_two_tellurone_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Te](=O)(=O)C[Te](=O)(=O)C")


def test_methanetelluronic_acid():
    assert smiles_to_iupac("C[Te](=O)(=O)O") == "methanetelluronic acid"


def test_pent_4_ene_1_telluronic_acid():
    assert smiles_to_iupac("C=CCCC[Te](=O)(=O)O") == "pent-4-ene-1-telluronic acid"


def test_ditelluronic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)(=O)C[Te](=O)(=O)O")


def test_ring_telluronic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)(=O)C1CCCCC1")


def test_telluronic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)(=O)CCO")


def test_benzenetelluronic_acid():
    assert smiles_to_iupac("c1ccccc1[Te](=O)(=O)O") == "benzenetelluronic acid"


def test_substituted_benzenetelluronic_acid():
    assert smiles_to_iupac("Cc1ccccc1[Te](=O)(=O)O") == "2-methylbenzenetelluronic acid"
    assert (
        smiles_to_iupac("Cc1ccc(cc1)[Te](=O)(=O)O") == "4-methylbenzenetelluronic acid"
    )  # PubChem PUG REST


def test_acyclic_telluronic_acid_stereocenter_with_coexisting_substituent():
    assert smiles_to_iupac("C[C@@H](Cl)[Te](=O)(=O)O") == "(1S)-1-chloroethanetelluronic acid"


def test_tellurourea():
    assert smiles_to_iupac("NC(=[Te])N") == "tellurourea"


def test_n_n_dimethyltellurourea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=[Te])N") == "N,N-dimethyltellurourea"


def test_n_ethyl_n_methyltellurourea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=[Te])N") == "N-ethyl-N-methyltellurourea"


def test_different_substituents_on_different_nitrogens_not_supported__tellurourea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNC(=[Te])NC")


def test_unsaturated_n_substituent_not_supported__tellurourea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=[Te])N")


def test_n_phenyltellurourea():
    assert smiles_to_iupac("NC(=[Te])Nc1ccccc1") == "N-phenyltellurourea"


def test_n_n_prime_diphenyltellurourea():
    assert smiles_to_iupac("c1ccc(NC(=[Te])Nc2ccccc2)cc1") == "N,N'-diphenyltellurourea"


def test_substituted_phenyl_n_substituent_not_supported__tellurourea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=[Te])Nc1ccc(C)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported__tellurourea():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=[Te])N")


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


def test_ring_not_supported__telluroxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=[Te]1CCCCC1")


def test_unsaturated_chain_not_supported__telluroxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[Te](=O)C")


def test_two_telluroxide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Te](=O)C[Te](=O)C")
