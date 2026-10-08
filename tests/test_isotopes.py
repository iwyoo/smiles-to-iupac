import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[2H]CC", "(2H1)ethane", id="one_deuterium_on_ethane"),
        pytest.param("[3H]C", "(3H1)methane", id="tritiated_methane"),
    ],
)
def test_ethane_deuterium_and_tritium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[2H]C([2H])(Cl)Cl", "dichloro(2H2)methane", id="dichlorodideuteromethane_name"),
        pytest.param("[14CH4]", "(14C)methane", id="carbon_14_methane_name"),
        pytest.param("[2H][13CH3]", "(13C,2H1)methane", id="deuterium_and_carbon_isotope_together_methane_name"),
        pytest.param("C[14CH2]C([2H])C", "(2-14C,3-2H1)butane", id="deuterium_and_carbon_isotope_together_chain_name"),
    ],
)
def test_dichlorodideuteromethane_name_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotopically_labeled_halogen_name():
    assert smiles_to_iupac("[2H]C([37Cl])") == "(37Cl)chloro(2H1)methane"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[14CH2]CC", "(2-14C)butane", id="carbon_14_butane_name"),
        pytest.param("FC(F)(F)C[2H]", "1,1,1-trifluoro(2-2H1)ethane", id="trifluoro_deuterio_ethane_name"),
    ],
)
def test_carbon_14_butane_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_single_halogen_ethane_cites_every_locant():
    assert smiles_to_iupac("FCC[2H]") == "1-fluoro(2-2H1)ethane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[13CH2]O", "(1-13C)ethan-1-ol"),
        ("[13CH3]O", "(13C)methanol"),
    ],
)
def test_isotope_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OCCC[18OH]", id="multiple_hydroxyls_raises"),
    ],
)
def test_isotope_alcohol_cases_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[13CH2][18OH]", "(1-13C)ethan-1-(18O)ol", id="combined_carbon_and_oxygen_isotope"),
        pytest.param("[18OH]CCC", "propan-1-(18O)ol", id="longer_chain"),
        pytest.param("[15OH]C", "methan(15O)ol", id="short_lived_oxygen_nuclide"),
        pytest.param("CC(O[2H])CC", "butan-2-(2H)ol", id="deuterated_hydroxyl_on_a_chain"),
    ],
)
def test_isotope_alcohol_oxygen_inserted_before_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotope_alcohol_17o_name():
    assert smiles_to_iupac("[17OH]C") == "methan(17O)ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[18O])O", "(18O)acetic acid"),
    ],
)
def test_isotope_carboxylic_acid_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_chain_unsaturation_alongside_isotope_carboxylic_acid():
    assert smiles_to_iupac("[13CH3]C(=O)OC=C") == "ethenyl (2-13C)acetate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[18O])C", "propan-2-(18O)one"),
        ("CC(=O)[14CH3]", "(1-14C)propan-2-one"),
    ],
)
def test_isotope_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=[13C]1CCCCC1", "(1-13C)cyclohexan-1-one"),
        ("OC1CCCC[13CH2]1", "(2-13C)cyclohexan-1-ol"),
        ("O=C(O)C1CCC[13CH2]C1", "(3-13C)cyclohexane-1-carboxylic acid"),
    ],
)
def test_locants_are_all_cited_once_a_nuclide_needs_one(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[2H]C1CC1", "(2H1)cyclopropane"),
        ("[13CH2]1CCCCC1", "(13C)cyclohexane"),
        ("[13CH2]1CC([2H])CCC1", "(1-13C,3-2H1)cyclohexane"),
        ("CCN(C)C[13CH3]", "N-ethyl-N-methyl(2-13C)ethan-1-amine"),
        ("CCN(CC)CC[13CH3]", "N,N-diethyl(3-13C)propan-1-amine"),
        ("[2H]N(C)CC", "N-methyl(N-2H)ethanamine"),
        ("CC[NH+](C)[2H]", "N-methyl(N-2H)ethanaminium"),
        ("C[NH2+][2H]", "(N-2H1)methanaminium"),
        ("C[15NH3+]", "(15N)methanaminium"),
        ("OCc1ccc([2H])c(Cl)c1", "[3-chloro(4-2H)phenyl]methanol"),
        ("O=C(C)c1cc([2H])c(C)cc1", "1-[4-methyl(3-2H)phenyl]ethan-1-one"),
        ("CC(=O)Oc1ccc(C)cc1[2H]", "4-methyl(2-2H)phenyl acetate"),
        ("CC(=O)Oc1ccc([2H])cc1", "(4-2H)phenyl acetate"),
        ("[2H]c1ccc2ccccc2c1", "(2-2H)naphthalene"),
        ("[2H]c1cccc2ccccc12", "(1-2H)naphthalene"),
        ("[2H]c1c([2H])cc2ccccc2c1", "(2,3-2H2)naphthalene"),
        ("[2H]c1ccc2ccccc2n1", "(2-2H)quinoline"),
        ("Cc1ccc2ccccc2c1[2H]", "2-methyl(1-2H)naphthalene"),
        ("Clc1ccc2cc([2H])ccc2c1", "2-chloro(6-2H)naphthalene"),
        ("[2H]c1ccc2ccccc2[n+]1C", "1-methyl(2-2H)quinolin-1-ium"),
        ("[2H]c1cc[n+](C)cc1", "1-methyl(4-2H)pyridin-1-ium"),
        ("C[n+]1ccccc1[2H]", "1-methyl(2-2H)pyridin-1-ium"),
        ("[2H]c1cccc[n+]1[O-]", "(2-2H)pyridine 1-oxide"),
        ("[2H]c1ccc[nH+]c1", "(3-2H)pyridin-1-ium"),
        ("[13CH2]CCO", "3-hydroxy(1-13C)propyl"),
        ("C[13CH2][CH2]", "(2-13C)propyl"),
        ("[CH2]C(=O)C[13CH3]", "2-oxo(4-13C)butyl"),
        ("[2H]C([2H])([2H])[CH2]", "(2,2,2-2H3)ethyl"),
    ],
)
def test_isotope_ring_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[3H]C1CCCCC1", "(3H1)cyclohexane", id="tritium_on_a_ring"),
        pytest.param("[2H]C1CCCCC1C", "1-methyl(2-2H1)cyclohexane", id="ring_label_beside_a_substituent"),
        pytest.param("[2H]c1ccccc1", "(2H)benzene", id="one_deuterium_on_benzene"),
        pytest.param("[2H]c1ccccc1[2H]", "(1,2-2H2)benzene", id="two_deuteriums_need_locants"),
        pytest.param("[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]", "(2H6)benzene", id="every_position_modified"),
        pytest.param("Cc1ccccc1[2H]", "1-methyl(2-2H)benzene", id="benzene_locants_cited_for_a_ring_label"),
        pytest.param("[2H]c1ccc(cc1)C(=O)O", "(4-2H)benzoic acid", id="ring_label_with_a_principal_group"),
    ],
)
def test_isotope_on_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC=[18O]", "(18O)formic acid"),
    ],
)
def test_isotope_descriptor_before_retained_acid_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[13CH3][14CH2]O", "(2-13C,1-14C)ethan-1-ol", id="two_carbon_nuclides_in_one_series"),
        pytest.param("CC(=O)O[13CH3]", "(13C)methyl acetate", id="isotopic_methyl_group_of_an_ester"),
        pytest.param("O=C(O[13CH2]C)c1ccccc1", "(1-13C)ethyl benzoate", id="isotopic_alkyl_chain_of_an_ester"),
        pytest.param("O=C(c1ccccc1)[13CH3]", "1-phenyl(2-13C)ethan-1-one", id="isotope_in_parent_after_substituent_prefix"),
        pytest.param("Cc1cccnc1[13CH3]", "2-(13C)methyl-3-methylpyridine", id="isotopic_substituent_cited_before_unmodified_one"),
        pytest.param("CC(=O)Nc1ccc([131I])cc1", "N-[4-(131I)iodophenyl]acetamide", id="isotopic_halogen_inside_a_ring_substituent"),
    ],
)
def test_isotopic_descriptor_follows_the_unmodified_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles", ["CC(=O)Nc1ccccc1C(=O)[13CH3]"])
def test_unplaceable_isotope_label_is_never_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)N[2H]", "(N-2H1)acetamide", id="amide_nitrogen_deuterium"),
        pytest.param("CC(=O)[15NH2]", "(15N)acetamide", id="amide_nitrogen_nuclide"),
        pytest.param("[2H]N([2H])c1ccccc1", "(N-2H2)aniline", id="amine_nitrogen_two_deuteriums"),
        pytest.param("CCN([2H])[2H]", "(N-2H2)ethanamine", id="amine_on_a_chain"),
        pytest.param("CC(=O)O[2H]", "(O-2H)acetic acid", id="acid_hydroxyl_deuterium"),
        pytest.param("CCC(=O)O[2H]", "(O-2H)propanoic acid", id="acid_hydroxyl_deuterium_on_a_chain"),
        pytest.param("Clc1ccc(N([2H])[2H])cc1", "4-chloro(N-2H2)aniline", id="group_label_after_prefixes"),
        pytest.param("CC(=[18O])C(C)C", "3-methylbutan-2-(18O)one", id="ketone_oxygen_nuclide"),
        pytest.param("CCC#[15N]", "(15N)propanenitrile", id="nitrile_nitrogen_nuclide"),
        pytest.param("ClCC#[15N]", "chloro(15N)acetonitrile", id="nitrile_nitrogen_after_prefixes"),
        pytest.param("CC(=[18O])[13CH3]", "(1-13C)propan-2-(18O)one", id="skeleton_and_group_nuclides"),
    ],
)
def test_isotope_on_characteristic_group_atom(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)Nc1ccccc1C[13CH3]", "N-[2-(2-13C)ethylphenyl]acetamide", id="labelled_alkyl_on_a_ring_substituent"),
        pytest.param("CC(=O)c1ccc([18F])cc1", "1-[4-(18F)fluorophenyl]ethan-1-one", id="labelled_halogen_on_a_ring_substituent"),
        pytest.param("O=C(O)c1ccc(cc1)C([2H])([2H])[2H]", "4-(2H3)methylbenzoic acid", id="labelled_methyl_beside_a_ring_principal_group"),
        pytest.param("C([2H])([2H])([2H])Oc1ccc(cc1)C(=O)O", "4-(2H3)methoxybenzoic acid", id="labelled_methoxy_beside_a_ring_principal_group"),
        pytest.param("CC(=O)Nc1ccccc1OC[13CH3]", "N-{2-[(2-13C)ethoxy]phenyl}acetamide", id="labelled_ethoxy_on_a_ring_substituent"),
        pytest.param("OC(=O)CC[13CH2]Cl", "4-chloro(4-13C)butanoic acid", id="descriptor_after_prefixes_of_an_acid"),
        pytest.param("OC(=O)C[13CH](C)C", "3-methyl(3-13C)butanoic acid", id="labelled_branch_point_of_an_acid"),
        pytest.param("CC(=O)[14CH3]", "(1-14C)propan-2-one", id="lowest_locant_to_the_modified_atom"),
        pytest.param("FCC[2H]", "1-fluoro(2-2H1)ethane", id="every_locant_cited_for_a_modified_chain"),
        pytest.param("OC(=O)c1ccc(cc1)C([2H])([2H])CCl", "4-[2-chloro(1,1-2H2)ethyl]benzoic acid", id="labelled_chain_beside_a_prefix"),
        pytest.param("OC(=O)c1ccc(cc1)C(C)[13CH3]", "4-[(1-13C)propan-2-yl]benzoic acid", id="lowest_locant_on_a_branch_point_group"),
        pytest.param("OC(=O)c1ccc(cc1)C([2H])([2H])C([2H])([2H])[2H]", "4-(2H5)ethylbenzoic acid", id="every_position_of_a_substituent_modified"),
        pytest.param("OC(=O)CCc1c([2H])cccc1", "3-(2-2H)phenylpropanoic acid", id="labelled_phenyl_group"),
        pytest.param("OC(=O)CCc1c([2H])c([2H])c([2H])c([2H])c1[2H]", "3-(2H5)phenylpropanoic acid", id="fully_labelled_phenyl_group"),
        pytest.param("OC(=O)c1ccc(cc1)O[2H]", "4-(2H)hydroxybenzoic acid", id="deuterated_hydroxy_takes_no_subscript"),
        pytest.param("OC(=O)c1ccc(cc1)N([2H])[2H]", "4-(2H2)aminobenzoic acid", id="deuterated_amino_counts_its_atoms"),
    ],
)
def test_isotope_on_substituents_and_prefixed_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[13CH3]C[14CH2]C", "(1-13C,3-14C)butane", id="two_carbon_nuclides_in_one_series"),
        pytest.param("[13CH]([2H])1CCCCC1", "(1-13C,1-2H1)cyclohexane", id="carbon_and_deuterium_on_one_ring_atom"),
        pytest.param("CC=C[2H]", "(1-2H1)prop-1-ene", id="deuterated_alkene"),
        pytest.param("C=C([2H])C", "(2-2H)prop-1-ene", id="deuterium_on_an_alkene_carbon"),
        pytest.param("[2H]C#CC", "(1-2H)prop-1-yne", id="deuterated_alkyne"),
        pytest.param("C/C=C/C[2H]", "(2E)-(1-2H1)but-2-ene", id="alkene_configuration_with_a_label"),
        pytest.param("C[C@H]([2H])O", "(1S)-(1-2H1)ethan-1-ol", id="stereocentre_created_by_deuterium"),
        pytest.param("C[C@@H]([2H])O", "(1R)-(1-2H1)ethan-1-ol", id="opposite_isotopic_stereocentre"),
        pytest.param("[13CH3][C@H](Cl)C(=O)O", "(2S)-2-chloro(3-13C)propanoic acid", id="stereocentre_beside_a_labelled_acid"),
        pytest.param("[13CH3]C(=O)[C@H](Cl)C", "(3R)-3-chloro(1-13C)butan-2-one", id="stereocentre_beside_a_labelled_ketone"),
        pytest.param("C(=C\\C)/[2H]", "(1E)-(1-2H1)prop-1-ene", id="alkene_stereo_from_a_deuterium"),
        pytest.param("C(=C/C)/[2H]", "(1Z)-(1-2H1)prop-1-ene", id="alkene_stereo_from_a_deuterium_z"),
        pytest.param("C1(=CC=CC=C1)[13C]#N", "benzene(13C)carbonitrile", id="labelled_nitrile_carbon_on_a_ring"),
        pytest.param("C1(=CC=CC=C1)[13C](=O)O", "benzene(13C)carboxylic acid", id="labelled_carboxyl_carbon_on_a_ring"),
        pytest.param("C1(=CC=CC=C1)[13CH2]NN", "[phenyl(13C)methyl]hydrazine", id="labelled_benzyl_carbon_on_a_hydrazine"),
        pytest.param("[14C](=O)(O)C1(CCCCC1)[14C](=O)O", "cyclohexane-1,1-di[(14C)carboxylic acid]", id="identical_labelled_acids_multiplied"),
        pytest.param("[13C](=O)(O)C1(CCCCC1)[14C](=O)O", "1-(13C)carboxycyclohexane-1-(14C)carboxylic acid", id="different_labelled_acids_senior_nuclide_is_suffix"),
        pytest.param("OC(=O)C1(CCCCC1)[13C](=O)O[2H]", "1-carboxycyclohexane-1-(13C,2H)carboxylic acid", id="acid_with_more_modifications_is_suffix"),
        pytest.param("[2H]OC(=O)C1(CCCCC1)[13C](=O)O", "1-(2H)carboxycyclohexane-1-(13C)carboxylic acid", id="labelled_carboxy_prefix"),
        pytest.param("C(OC([2H])([2H])SCOO)[2H]", "{[(2H1)methoxy(2H2)methyl]sulfanyl}methaneperoxol", id="labelled_ether_in_a_sulfanyl_peroxol"),
    ],
)
def test_isotopes_with_unsaturation_and_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CCC(=O)[18O]CC", "O-ethyl propan(18O1)oate", id="labelled_ester_oxygen"),
        pytest.param("CCC(=[18O])OCC", "O-ethyl propan(18O1)oate", id="labelled_carbonyl_oxygen_of_an_ester"),
        pytest.param("COC(=O)O[2H]", "methyl (2H)hydrogen carbonate", id="labelled_hydrogen_of_a_carbonate_half_ester"),
        pytest.param("COC(=O)[18OH]", "O-methyl hydrogen (18O1)carbonate", id="labelled_hydroxy_oxygen_of_a_carbonate_half_ester"),
        pytest.param("CC[18O]C(=O)OC", "O-ethyl O-methyl (18O1)carbonate", id="labelled_bridging_oxygen_of_a_carbonate"),
        pytest.param("CCOC(=[18O])[18O]C", "O-ethyl 18O-methyl (18O2)carbonate", id="two_labelled_oxygens_of_a_carbonate"),
        pytest.param("CCOC(=[18O])OCC", "O,O-diethyl (18O1)carbonate", id="labelled_carbonyl_oxygen_of_a_symmetric_carbonate"),
        pytest.param("CC(=[18O])OC", "O-methyl acet(18O1)ate", id="labelled_acetate_carbonyl"),
        pytest.param("O=C([18O]C)c1ccccc1", "O-methyl benz(18O1)oate", id="labelled_benzoate_ester_oxygen"),
        pytest.param("OC(=O)c1ccc(cc1)C(C)(C)[13CH3]", "4-[2-methyl(1-13C)propan-2-yl]benzoic acid", id="modified_tert_butyl_loses_its_retained_name"),
        pytest.param("OC(=O)c1ccc(cc1)[13C](C)(C)C", "4-[2-methyl(2-13C)propan-2-yl]benzoic acid", id="modified_quaternary_carbon_of_tert_butyl"),
        pytest.param("[2H]C([2H])([2H])C", "(1,1,1-2H3)ethane", id="locants_kept_when_isomers_exist"),
        pytest.param("C(OC1=CC=CC=C1)([2H])([2H])[2H]", "(2H3)methoxybenzene", id="lone_modified_alkoxy_is_not_enclosed"),
        pytest.param("N[14CH2]C1(CCCC1)O", "1-[amino(14C)methyl]cyclopentan-1-ol", id="descriptor_inside_a_compound_prefix"),
        pytest.param("C(C([2H])[2H])C(CO)C(CCC)CC", "2-(2,2-2H2)ethyl-3-ethylhexan-1-ol", id="modified_prefix_cited_first"),
        pytest.param("[15NH]1C=CC2=CC=CC=C12", "(15N)-1H-indole", id="sole_heteroatom_needs_no_locant"),
        pytest.param("[15N]1=C(C=C(C=C1)[2H])[2H]", "(2,4-2H2,15N)pyridine", id="sole_heteroatom_after_hydrogen_locants"),
        pytest.param("CC([14CH2]C)([2H])[2H]", "(3-14C,2,2-2H2)butane", id="lowest_locants_to_all_nuclides_together"),
    ],
)
def test_isotope_on_ester_oxygens_and_tert_butyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CCC(CC)CC([2H])([2H])[2H]", "3-ethyl(1,1,1-2H3)pentane", id="chain_through_the_modified_branch"),
        pytest.param("CCC(C)(CC)CC([2H])([2H])[2H]", "3-ethyl-3-methyl(1,1,1-2H3)pentane", id="modified_chain_with_substituents"),
        pytest.param("[2H]C(CC)(CC)CC", "3-ethyl(3-2H)pentane", id="one_modified_atom_on_the_chain"),
        pytest.param("CCC(CC[13CH3])(CC)CC[14CH3]", "4,4-diethyl(7-13C,1-14C)heptane", id="both_chains_modified"),
        pytest.param("CCC(CC)(CC[2H])C[13CH3]", "3,3-diethyl(1-13C,5-2H1)pentane", id="more_nuclides_of_higher_atomic_number_first"),
    ],
)
def test_parent_chain_with_more_isotopic_modifications_is_senior(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[13CH3]C1=C(C=CC=C1)[13CH3]", "1,2-di[(13C)methyl]benzene", id="multiplied_isotopic_prefix"),
        pytest.param("[13CH3]c1ccc(cc1)[13CH3]", "1,4-di[(13C)methyl]benzene", id="multiplied_isotopic_prefix_para"),
        pytest.param("[13CH3]C1CCCCC1[13CH3]", "1,2-di[(13C)methyl]cyclohexane", id="multiplied_isotopic_prefix_on_a_ring"),
        pytest.param("[2H]C([2H])([2H])c1ccccc1C([2H])([2H])[2H]", "1,2-di[(2H3)methyl]benzene", id="multiplied_deuterated_prefix"),
        pytest.param("C(C)S[34S]SCCC(=O)O", "3-[ethyl(2-34S)trisulfanyl]propanoic acid", id="modified_atom_in_a_sulfanyl_chain"),
        pytest.param("CC[34S]SCCC(=O)O", "3-[ethyl(2-34S)disulfanyl]propanoic acid", id="modified_atom_in_a_disulfanyl_group"),
        pytest.param("OC(=O)CCSS[34S]CC", "3-[ethyl(3-34S)trisulfanyl]propanoic acid", id="modified_terminal_atom_of_a_sulfanyl_chain"),
    ],
)
def test_isotopic_descriptors_in_multiplied_prefixes_and_chalcogen_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[2H]c1ccc(Cc2ccc(cc2)[2H])c([2H])c1", "1-[(4-2H)phenylmethyl](2,4-2H2)benzene"),
        ("[2H]c1ccc(Cc2ccc([3H])cc2)cc1", "1-[(4-2H)phenylmethyl](4-3H)benzene"),
    ],
)
def test_ring_with_more_or_heavier_nuclides_is_the_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[2H]c1ccc(cc1)-c1ccccc1", "(4-2H)-1,1'-biphenyl"),
        ("[2H]c1ccc(cc1)-c1ccc([2H])cc1", "(4,4'-2H2)-1,1'-biphenyl"),
        ("OC(=O)c1ccc(cc1)-c1ccccc1[2H]", "(2'-2H)[1,1'-biphenyl]-4-carboxylic acid"),
        ("[2H]c1ccc2ccccc2c1-c1ccc2ccccc2c1", "(2-2H)-1,2'-binaphthalene"),
        ("C(c1ccc2ccccc2c1)c1ccc2cc([2H])ccc2c1", "2-[(naphthalen-2-yl)methyl](6-2H)naphthalene"),
    ],
)
def test_isotopic_modification_of_ring_assemblies_and_fused_parent_choice(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[2H]CC(C)CC[SiH3]", "[3-methyl(4-2H1)butyl]silane"),
        ("[SiH3]CCC([14CH3])C[2H]", "[3-(2H1)methyl(4-14C)butyl]silane"),
        ("CC(=O)NC(CO)C[18OH]", "N-[1-(18O)hydroxy-3-hydroxypropan-2-yl]acetamide"),
        ("CCCC(Br)COCC([81Br])CCC", "2-bromo-1-{[2-(81Br)bromopentyl]oxy}pentane"),
        ("OCC(CC(C[14CH3])CCC)CC(C[13CH3])CCC", "4-(2-13C)ethyl-2-[2-(2-14C)ethylpentyl]heptan-1-ol"),
        ("[13CH3]OCCNCC[18O]C", "2-(13C)methoxy-N-[2-(18O)methoxyethyl]ethan-1-amine"),
        ("[2H]C1CCCCC1(C)CC1CCCCC1", "1-(cyclohexylmethyl)-1-methyl(2-2H1)cyclohexane"),
        ("[14CH2]1CCCCC1CC1([2H])CCCCC1", "1-[(1-2H)cyclohexylmethyl](2-14C)cyclohexane"),
    ],
)
def test_substituent_chains_and_prefixes_are_chosen_by_isotopic_modification(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_deuterium_alkene_stereodescriptor_writes_no_rdkit_logs_to_stderr(capfd):
    assert smiles_to_iupac("C(=C\\C)/[2H]") == "(1E)-(1-2H1)prop-1-ene"
    assert capfd.readouterr().err == ""
