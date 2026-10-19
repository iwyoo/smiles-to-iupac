import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


# P-82.2.1, P-82.6.1.2: the carbonyl oxygen, the only oxygen of the molecule, is cited in front of the acyl compound without a locant
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=[18O])NN", "(18O)acetohydrazide", id="hydrazide"),
        pytest.param("CC(=[18O])N=[N+]=[N-]", "(18O)acetyl azide", id="acyl_azide"),
        pytest.param("CC(=[18O])SC", "S-methyl (18O)ethanethioate", id="thioester"),
        pytest.param("CC=[18O]", "(18O)acetaldehyde", id="retained_aldehyde_name"),
        pytest.param("[18O]=Cc1ccc(Cl)cc1", "4-chloro(18O)benzaldehyde", id="substituted_retained_aldehyde_name"),
    ],
)
def test_carbonyl_oxygen_of_an_acyl_compound(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[13CH2]=O", "(13C)formaldehyde", id="skeletal_carbon"),
        pytest.param("[2H]C=S", "(2H1)methanethial", id="one_hydrogen_of_two"),
        pytest.param("[2H]C#N", "(2H)formonitrile", id="only_hydrogen_takes_no_subscript"),
        pytest.param("[13CH]#[15N]", "(13C,15N)formonitrile", id="two_elements_in_alphabetical_order"),
    ],
)
def test_isotope_descriptor_before_a_one_carbon_retained_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[2H]CC", "(2H1)ethane", id="one_deuterium_on_ethane"),
        pytest.param("[3H]C", "(3H1)methane", id="tritiated_methane"),
        pytest.param("[1H]C([2H])([2H])[2H]", "(1H1,2H3)methane", id="protium_cited_before_deuterium"),
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
        ("[18OH]CCC[18OH]", "propane-1,3-di[(18O)ol]"),
    ],
)
def test_isotope_alcohol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[13CH2][18OH]", "(1-13C)ethan-1-(18O)ol", id="combined_carbon_and_oxygen_isotope"),
        pytest.param("[18OH]CCC", "propan-1-(18O)ol", id="longer_chain"),
        pytest.param("[15OH]C", "methan(15O)ol", id="short_lived_oxygen_nuclide"),
        pytest.param("CC(O[2H])CC", "butan-2-(2H)ol", id="deuterated_hydroxyl_on_a_chain"),
        pytest.param("[18OH]CC(CO)C1CCCCC1", "2-cyclohexyl-3-hydroxypropan-1-(18O)ol", id="labelled_hydroxyl_stays_the_suffix"),
        pytest.param("OCC([18OH])CO", "1,3-dihydroxypropan-2-(18O)ol", id="labelled_hydroxyl_between_two_unmodified_ones"),
        pytest.param("OC1CCC([18OH])CC1", "4-hydroxycyclohexan-1-(18O)ol", id="labelled_hydroxyl_on_a_ring"),
    ],
)
def test_isotope_alcohol_oxygen_inserted_before_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_isotope_alcohol_17o_name():
    assert smiles_to_iupac("[17OH]C") == "methan(17O)ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=[18O])O", "(O-18O)acetic acid"),
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
        ("[2H]CC(N)CNC", "N1-methyl(3-2H1)propane-1,2-diamine"),
        ("[2H]CC(CN)NC", "N2-methyl(3-2H1)propane-1,2-diamine"),
        ("[2H]CC(N)CN(C)C", "N1,N1-dimethyl(3-2H1)propane-1,2-diamine"),
        ("[2H]c1ccc(NC)c(N)c1", "N1-methyl(4-2H)benzene-1,2-diamine"),
        ("[2H]COc1cccc(NCC)c1N", "N1-ethyl-3-(2H1)methoxybenzene-1,2-diamine"),
        ("CC[NH+](C)[2H]", "N-methyl(N-2H)ethanaminium"),
        ("C[NH2+][2H]", "(N-2H1)methanaminium"),
        ("C[15NH3+]", "(15N)methanaminium"),
        ("[2H]C([2H])([2H])[N+](C)(C)CC[NH3+]", "N1-(2H3)methyl-N1,N1-dimethylethane-1,2-bis(aminium)"),
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
        ("[18O-][n+]1ccccc1", "pyridine 1-(18O)oxide"),
        ("[18O-][n+]1cc[n+]([18O-])cc1", "pyrazine 1,4-di[(18O)oxide]"),
        ("[O-][n+]1cc[n+]([O-])c([2H])c1", "(2-2H)pyrazine 1,4-dioxide"),
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
        pytest.param("[2H]c1ccccc1", "(2H1)benzene", id="one_deuterium_on_benzene"),
        pytest.param("[2H]c1ccccc1[2H]", "(1,2-2H2)benzene", id="two_deuteriums_need_locants"),
        pytest.param("[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]", "(2H6)benzene", id="every_position_modified"),
        pytest.param("Cc1ccccc1[2H]", "1-methyl(2-2H)benzene", id="benzene_locants_cited_for_a_ring_label"),
        pytest.param("[2H]c1ccc(cc1)C(=O)O", "(4-2H)benzoic acid", id="ring_label_with_a_principal_group"),
        pytest.param("[2H]C([2H])([2H])C#N", "(2H3)acetonitrile", id="every_hydrogen_modified_beside_a_position_without_hydrogen"),
        pytest.param("[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]", "(2H6)propan-2-one", id="every_hydrogen_modified_keeps_the_locant_of_the_name"),
        pytest.param("[2H]C([2H])([2H])C(=O)O", "(2,2,2-2H3)acetic acid", id="locants_kept_when_a_hydroxy_hydrogen_remains"),
    ],
)
def test_isotope_on_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]C([2H])(C(=O)O)c1ccccc1", "phenyl(2,2-2H2)acetic acid", id="acid"),
        pytest.param("[2H][C@@H](C(=O)O)C1CCCC1", "(2R)-cyclopentyl(2-2H1)acetic acid", id="stereo_prefix_stays_in_front"),
        pytest.param("[2H]C([2H])(C(=O)Cl)c1ccccc1", "phenyl(2,2-2H2)acetyl chloride", id="acyl_halide"),
        pytest.param("[2H]C([2H])(C(=O)OC)c1ccccc1", "methyl phenyl(2,2-2H2)acetate", id="ester"),
        pytest.param("[2H]C([2H])(C(=O)O)c1ccc(Cl)cc1", "(4-chlorophenyl)(2,2-2H2)acetic acid", id="compound_prefix"),
    ],
)
def test_descriptor_follows_the_prefixes_of_a_substituted_retained_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)[15NH]CC(=O)O", "[(15N)acetamido]acetic acid", id="nitrogen_of_an_acylamino_prefix"),
        pytest.param("CC(=O)[15NH][C@@H](C)C(=O)O", "(2S)-2-[(15N)acetamido]propanoic acid", id="stereocentre_beside_the_labelled_prefix"),
        pytest.param("CC(=O)[15N](C)CC(=O)O", "[N-methyl(15N)acetamido]acetic acid", id="substituted_nitrogen_of_the_prefix"),
    ],
)
def test_nuclide_on_the_nitrogen_of_an_acylamino_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]C([2H])([2H])C(=O)Nc1cccc(C(=O)O)c1", "3-[(2H3)acetamido]benzoic acid", id="acylamino_prefix_beside_an_acid"),
        pytest.param("[2H]C([2H])([2H])C(=O)Nc1ccc(C(=O)OC)cc1", "methyl 4-[(2H3)acetamido]benzoate", id="acylamino_prefix_beside_an_ester"),
        pytest.param("[2H]C([2H])([2H])C(=O)OCC(=O)O", "[(2H3)acetyloxy]acetic acid", id="acyloxy_prefix_beside_an_acid"),
    ],
)
def test_nuclides_on_the_acyl_group_of_a_prefix_beside_an_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=[18O])N", "(18O)acetamide", id="amide_oxygen"),
        pytest.param("CC(=[18O])NC", "N-methyl(18O)acetamide", id="substituted_amide_oxygen"),
        pytest.param("CC(=[18O])Cl", "(18O)acetyl chloride", id="acid_halide_oxygen"),
        pytest.param("CC(=[18O])CC(=[18O])C", "pentane-2,4-di[(18O)one]", id="identically_labelled_ketones"),
        pytest.param("CC(=[18O])N[2H]", "(N-2H1,18O)acetamide", id="amide_oxygen_and_nitrogen_hydrogen"),
        pytest.param("CC(=O)[15NH][2H]", "(N-2H1,15N)acetamide", id="nitrogen_nuclide_and_hydrogen_nuclide"),
        pytest.param("OC(=[18O])CC(=O)[18OH]", "(O,O'-18O2)propanedioic acid", id="one_oxygen_in_each_acid_group"),
    ],
)
def test_nuclide_on_the_carbonyl_oxygen_of_an_amide_acid_halide_or_several_ketones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]C([2H])([2H])C(=O)NC(C)C(=O)O", "N-[(2H3)acetyl]alanine", id="every_hydrogen_of_the_acyl_group_modified"),
        pytest.param("[3H]CC(=O)NC(C)C(=O)O", "N-[(2-3H1)acetyl]alanine", id="locant_kept_beside_remaining_hydrogens"),
        pytest.param(
            "[2H]C([2H])([2H])C(=O)N[C@H](CC(C)C)C(=O)O", "N-[(2H3)acetyl]-D-leucine", id="configuration_unchanged_by_the_acyl_label"
        ),
    ],
)
def test_nuclide_on_the_acyl_group_of_an_acylated_amino_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC=[18O]", "(O-18O)formic acid"),
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
        pytest.param("[2H]C([2H])([2H])OC(=O)/C=C/C(C)C", "(2H3)methyl (2E)-4-methylpent-2-enoate", id="alkyl_label_stays_off_the_acid_branch"),
        pytest.param("[2H]CC(=O)OC(C)=O", "(2-2H1)acetic acetic anhydride", id="anhydride_word_with_the_modified_acid_cited_first"),
        pytest.param("[2H]CC(=O)SC", "S-methyl (2-2H1)ethanethioate", id="thioester_acyl_part"),
        pytest.param("[2H]CSC(C)=O", "S-(2H1)methyl ethanethioate", id="thioester_organyl_part"),
        pytest.param("[2H][C@H](C)C(=O)OC", "methyl (2R)-(2-2H1)propanoate", id="isotopic_stereocentre_of_the_acid_part"),
        pytest.param("[2H]C(=O)CNC(=O)OC", "methyl [2-oxo(2-2H)ethyl]carbamate", id="carbamate_substituent_on_nitrogen"),
    ],
)
def test_isotopic_descriptor_follows_the_unmodified_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_name_of_an_unlabelled_molecule_ignores_earlier_labelled_ones():
    smiles_to_iupac("[2H]C([2H])([2H])OC(=O)/C=C/C(C)C")
    assert smiles_to_iupac("CCOC(=O)/C=C/C(C)C") == "ethyl (2E)-4-methylpent-2-enoate"


@pytest.mark.parametrize("smiles", ["[2H]I(C)C"])
def test_unplaceable_isotope_label_is_never_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)N[2H]", "(N-2H1)acetamide", id="amide_nitrogen_deuterium"),
        pytest.param("[2H]Nc1ccccn1", "(N-2H1)pyridin-2-amine", id="amino_group_on_an_azine"),
        pytest.param("[2H]NC(=S)CC", "(N-2H1)propanethioamide", id="thioamide_nitrogen"),
        pytest.param("[2H]NS(C)(=O)=O", "(N-2H1)methanesulfonamide", id="sulfonamide_nitrogen"),
        pytest.param("[2H]SC(C)=O", "(S-2H)ethanethioic S-acid", id="thioacid_sulfur"),
        pytest.param("[2H]NNC(C)=O", "(N'-2H1)acetohydrazide", id="terminal_hydrazide_nitrogen_takes_the_primed_locant"),
        pytest.param("CS(=O)(=O)N([2H])N", "(N-2H)methanesulfonohydrazide", id="sulfonohydrazide_nitrogen_beside_the_centre_is_unprimed"),
        pytest.param("CS(=N[2H])(=N)NC", "N-methyl(N'-2H)methanesulfonodiimidamide", id="modified_imido_nitrogen_of_a_diimidamide_takes_the_lower_prime"),
        pytest.param("CS(=NC)(=N[2H])N", "N'-methyl(N''-2H)methanesulfonodiimidamide", id="substituent_decides_the_prime_before_the_modification"),
        pytest.param("CS(=N[2H])(=N[2H])N", "(N',N''-2H2)methanesulfonodiimidamide", id="both_imido_nitrogens_of_a_diimidamide"),
        pytest.param("[2H]N=C(C)N", "(N'-2H)ethanimidamide", id="imino_nitrogen_of_an_amidine_takes_the_primed_locant"),
        pytest.param("[2H]N=C(C)N[2H]", "(N,N'-2H2)ethanimidamide", id="both_nitrogens_of_an_amidine_in_one_descriptor"),
        pytest.param("CC(=O)N[15NH2]", "(N'-15N)acetohydrazide", id="terminal_hydrazide_nitrogen_nuclide_takes_the_primed_locant"),
        pytest.param("[2H]N1Cc2ccccc2C1", "2,3-dihydro(2-2H)-1H-isoindole", id="ring_nitrogen_that_carries_a_hydro_prefix"),
        pytest.param("[2H]OP(=O)(O)O", "(O-2H)phosphoric acid", id="phosphoric_acid_oxygen"),
        pytest.param("[2H]OC(=O)NC", "methyl(O-2H)carbamic acid", id="carbamic_acid_oxygen"),
        pytest.param("[2H]OB(O)CC", "ethyl(O-2H)boronic acid", id="boronic_acid_oxygen"),
        pytest.param("[2H]OC(=O)O", "(O-2H)carbonic acid", id="carbonic_acid_oxygen"),
        pytest.param("[2H]OP(=O)(O[2H])O[2H]", "(O,O',O''-2H3)phosphoric acid", id="phosphoric_acid_three_oxygens"),
        pytest.param("[2H]OP(=O)(O[2H])OP(=O)(O)O", "(O1,O1'-2H2)diphosphoric acid", id="diphosphoric_acid_two_oxygens_on_one_centre"),
        pytest.param("[2H]OC#N", "(O-2H)cyanic acid", id="cyanic_acid_oxygen"),
        pytest.param("[2H]NC#N", "(N-2H1)cyanamide", id="cyanamide_nitrogen"),
        pytest.param("[2H]NNC", "1-methyl(2-2H1)hydrazine", id="substituted_hydrazine_cites_every_locant"),
        pytest.param("[2H]P([2H])C", "methyl(2H2)phosphane", id="substituted_mononuclear_hydride_needs_no_locant"),
        pytest.param("CC(=O)[15NH2]", "(15N)acetamide", id="amide_nitrogen_nuclide"),
        pytest.param("[2H]N([2H])c1ccccc1", "(N-2H2)aniline", id="amine_nitrogen_two_deuteriums"),
        pytest.param("CCN([2H])[2H]", "(N-2H2)ethanamine", id="amine_on_a_chain"),
        pytest.param("CCCOC(=O)N[2H]", "propyl (N-2H1)carbamate", id="carbamate_nitrogen_deuterium"),
        pytest.param("[2H]N(C)C(=O)OC", "methyl N-methyl(N-2H)carbamate", id="carbamate_nitrogen_deuterium_beside_a_substituent"),
        pytest.param("[2H]NC(N)=N", "(N-2H1)guanidine", id="guanidine_amino_nitrogen"),
        pytest.param("[2H]N=C(N)N", "(N''-2H)guanidine", id="guanidine_imino_nitrogen_takes_the_double_primed_locant"),
        pytest.param("CNC(=N)N[2H]", "N-methyl(N'-2H1)guanidine", id="guanidine_labelled_amino_nitrogen_beside_a_substituent"),
        pytest.param("[2H]NC(N)=N[2H]", "(N,N''-2H2)guanidine", id="guanidine_two_nitrogens_in_one_descriptor"),
        pytest.param("[2H]O[Si]([2H])(C)C", "dimethyl(1,O-2H2)silanol", id="silanol_hydrogen_of_the_centre_cited_with_the_hydroxy_hydrogen"),
        pytest.param("[2H]OOS(=O)(=O)O", "(O-2H)sulfuro(peroxoic) acid", id="hydroperoxy_oxygen_cited_through_the_peroxy_group"),
        pytest.param("C[Si](C)(C)N[2H]", "1,1,1-trimethyl(N-2H1)silanamine", id="amino_group_of_a_silane"),
        pytest.param("[2H]NC(=O)NC", "N-methyl(N'-2H1)urea", id="urea_labelled_nitrogen_is_the_primed_one_when_substituted"),
        pytest.param("[2H]NC(N)=O", "(N-2H1)urea", id="urea_labelled_nitrogen_is_unprimed_when_equivalent"),
        pytest.param("CN[13C](=[18O])N", "N-methyl(13C,18O)urea", id="urea_carbon_and_oxygen_need_no_locant"),
        pytest.param("CNC(=O)[15NH2]", "N-methyl(N'-15N)urea", id="urea_nitrogen_nuclide_is_cited_like_its_hydrogen"),
        pytest.param("CNC(=[34S])[15NH2]", "N-methyl(N'-15N,34S)thiourea", id="thiourea_carbon_nitrogen_and_sulfur"),
        pytest.param("NC(=O)N[15NH2]", "(2-15N)hydrazinecarboxamide", id="semicarbazide_nitrogen"),
        pytest.param("CCN(C)[13C](=O)O", "ethyl(methyl)(13C)carbamic acid", id="carbamic_acid_carbon_nuclide_after_prefixes"),
        pytest.param("N#[13C]O", "(13C)cyanic acid", id="cyanic_acid_carbon_nuclide"),
        pytest.param("CC(=O)O[2H]", "(O-2H)acetic acid", id="acid_hydroxyl_deuterium"),
        pytest.param("CCC(=O)O[2H]", "(O-2H)propanoic acid", id="acid_hydroxyl_deuterium_on_a_chain"),
        pytest.param("Clc1ccc(N([2H])[2H])cc1", "4-chloro(N-2H2)aniline", id="group_label_after_prefixes"),
        pytest.param("OC(=O)CCSc1ccccc1N([2H])[2H]", "3-{[2-(2H2)aminophenyl]sulfanyl}propanoic acid", id="amino_prefix_on_a_ring_of_a_substituent"),
        pytest.param(
            "[2H]N([2H])C1C=C(C=C(C1)C)SC1=C(C(CC=C1)N([2H])[2H])C",
            "3-{[3-(2H2)amino-5-methylcyclohexa-1,5-dien-1-yl]sulfanyl}-2-methyl(N-2H2)cyclohexa-2,4-dien-1-amine",
            id="amine_of_the_parent_beside_an_amino_prefix_of_a_substituent",
        ),
        pytest.param("[2H]N(C)CCO", "2-[methyl(2H)amino]ethan-1-ol", id="substituted_amino_prefix_beside_a_hydroxy_group"),
        pytest.param("[2H]Nc1ccc(N[2H])cc1", "(N1,N4-2H2)benzene-1,4-diamine", id="nitrogens_of_a_multiplied_amine_suffix_take_the_locant_of_their_carbon"),
        pytest.param("[2H]NC(=O)CCC(=O)N[2H]", "(N1,N4-2H2)butanediamide", id="nitrogens_of_a_diamide_take_the_locant_of_their_acyl_carbon"),
        pytest.param("CN([2H])C(=O)C1CCC(CC1)C(=O)N[2H]", "N1-methyl(N1,N4-2H2)cyclohexane-1,4-dicarboxamide", id="nitrogens_of_a_dicarboxamide_take_the_locant_of_the_ring_atom"),
        pytest.param("[2H]N(C)CCN", "2-amino-N-methyl(N-2H)ethan-1-amine", id="modified_substituted_amine_is_the_suffix_beside_an_unmodified_amine"),
        pytest.param("[2H]NCCN(C)C", "2-(dimethylamino)(N-2H1)ethan-1-amine", id="modified_amine_is_the_suffix_beside_an_unmodified_substituted_amine"),
        pytest.param("[2H]N(C)CCN(C)[2H]", "N1,N2-dimethyl(N1,N2-2H2)ethane-1,2-diamine", id="substituted_amines_modified_alike_are_multiplied"),
        pytest.param("[2H]NCCN([2H])C", "N1-methyl(N1,N2-2H2)ethane-1,2-diamine", id="primary_and_substituted_amine_modified_alike_are_multiplied"),
        pytest.param("CC[C@H](N(C)[2H])C", "(2R)-N-methyl(N-2H)butan-2-amine", id="stereocentre_of_a_modified_substituted_amine_parent"),
        pytest.param("CC[C@@]([2H])(N(C)[2H])C", "(2S)-N-methyl(2,N-2H2)butan-2-amine", id="nuclide_makes_the_centre_of_a_substituted_amine_parent_stereogenic"),
        pytest.param("[2H]NC[C@H](N(C)[2H])C", "(2R)-N2-methyl(N1,N2-2H2)propane-1,2-diamine", id="stereocentre_beside_amines_modified_alike"),
        pytest.param("CN([3H])CC(=O)O", "[methyl(3H)amino]acetic acid", id="substituted_amino_prefix_beside_a_carboxy_group"),
        pytest.param("OCC[15NH]c1ccc(C)cc1", "2-[4-methyl(15N)anilino]ethan-1-ol", id="skeletal_nitrogen_of_an_anilino_prefix"),
        pytest.param("OCCN([2H])c1ccccc1", "2-[(N-2H)anilino]ethan-1-ol", id="hydrogen_of_an_anilino_prefix_takes_the_letter_locant"),
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
        pytest.param("CC(=O)Nc1ccccc1C(=O)[13CH3]", "N-{2-[(2-13C)acetyl]phenyl}acetamide", id="labelled_acyl_on_a_ring_substituent"),
        pytest.param("OC(=O)CC[13CH2]Cl", "4-chloro(4-13C)butanoic acid", id="descriptor_after_prefixes_of_an_acid"),
        pytest.param("OC(=O)C[13CH](C)C", "3-methyl(3-13C)butanoic acid", id="labelled_branch_point_of_an_acid"),
        pytest.param("CC(=O)[14CH3]", "(1-14C)propan-2-one", id="lowest_locant_to_the_modified_atom"),
        pytest.param("FCC[2H]", "1-fluoro(2-2H1)ethane", id="every_locant_cited_for_a_modified_chain"),
        pytest.param("[2H]CS(=O)(=O)CC", "[(2H1)methanesulfonyl]ethane", id="alkane_of_a_sulfonyl_prefix"),
        pytest.param("[2H]CS(C)=O", "(methanesulfinyl)(2H1)methane", id="labelled_alkyl_of_a_symmetrical_sulfoxide_is_the_parent"),
        pytest.param("[2H]CC(=O)OOC", "methyl (2-2H1)ethaneperoxoate", id="acid_part_of_a_peroxy_ester"),
        pytest.param("OC(=O)c1ccc(cc1)C([2H])([2H])CCl", "4-[2-chloro(1,1-2H2)ethyl]benzoic acid", id="labelled_chain_beside_a_prefix"),
        pytest.param("OC(=O)c1ccc(cc1)C(C)[13CH3]", "4-[(1-13C)propan-2-yl]benzoic acid", id="lowest_locant_on_a_branch_point_group"),
        pytest.param("OC(=O)c1ccc(cc1)C([2H])([2H])C([2H])([2H])[2H]", "4-(2H5)ethylbenzoic acid", id="every_position_of_a_substituent_modified"),
        pytest.param("OC(=O)CCc1c([2H])cccc1", "3-(2-2H)phenylpropanoic acid", id="labelled_phenyl_group"),
        pytest.param("OC(=O)CCc1c([2H])c([2H])c([2H])c([2H])c1[2H]", "3-(2H5)phenylpropanoic acid", id="fully_labelled_phenyl_group"),
        pytest.param("C[C@H]([15NH2])C(=O)O", "(2S)-2-(15N)aminopropanoic acid", id="labelled_amino_acid_leaves_the_retained_name"),
        pytest.param("COC(=O)C(C)[15NH2]", "methyl 2-(15N)aminopropanoate", id="labelled_amino_acid_ester"),
        pytest.param("c1cnccc1C1SCC[14CH2]1", "4-[(3-14C)thiolan-2-yl]pyridine", id="labelled_heterocyclic_group"),
        pytest.param("OC(=O)Cc1c(Cl)ccc2c1cccc2[2H]", "[2-chloro(5-2H)naphthalen-1-yl]acetic acid", id="labelled_fused_group_after_its_prefixes"),
        pytest.param("OC(=O)Cc1cc[15n]cc1[2H]", "[(3-2H,15N)pyridin-4-yl]acetic acid", id="sole_ring_heteroatom_takes_no_locant"),
        pytest.param("OC(=O)Cc1cc2cc([2H])ccc2[nH]1", "[(5-2H)-1H-indol-2-yl]acetic acid", id="hyphen_before_indicated_hydrogen_of_a_group"),
        pytest.param(
            "[2H]c1c([2H])c([2H])c(Cc2ccccn2)c2c([2H])c([2H])c([2H])c([2H])c12",
            "2-{[(2H7)naphthalen-1-yl]methyl}pyridine",
            id="every_position_of_a_fused_group_modified",
        ),
        pytest.param(
            "N[C@@H](CCCN[14C](N)=[15NH])C(=O)O",
            "(2S)-2-amino-5-{[(14C,N′-15N)carbamimidoyl]amino}pentanoic acid",
            id="carbamimidoyl_carbon_unlocated_nitrogens_take_letter_locants",
        ),
        pytest.param("OC(=O)CCCNC(N[2H])=N", "4-{[(N-2H1)carbamimidoyl]amino}butanoic acid", id="deuterated_carbamimidoyl_nitrogen"),
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
        pytest.param("[2H][C@H](c1ccccc1)OC(C)=O", "(R)-phenyl(2H1)methyl acetate", id="isotopic_stereocentre_of_a_benzyl_group"),
        pytest.param("[2H][C@H]1CCCC[C@@H]1NC(N)=O", "[(1S,2S)-(2-2H1)cyclohexyl]urea", id="isotopic_stereocentre_of_a_ring_group"),
        pytest.param("C[C@H]([2H])O", "(1S)-(1-2H1)ethan-1-ol", id="stereocentre_created_by_deuterium"),
        pytest.param("C[C@@H]([2H])O", "(1R)-(1-2H1)ethan-1-ol", id="opposite_isotopic_stereocentre"),
        pytest.param("[13CH3][C@H](Cl)C(=O)O", "(2S)-2-chloro(3-13C)propanoic acid", id="stereocentre_beside_a_labelled_acid"),
        pytest.param("[13CH3]C(=O)[C@H](Cl)C", "(3R)-3-chloro(1-13C)butan-2-one", id="stereocentre_beside_a_labelled_ketone"),
        pytest.param("C[C@@H](C[C@@H](C)[2H])[3H]", "(2S,4R)-(4-2H1,2-3H1)pentane", id="higher_nuclide_numbered_before_the_descriptors"),
        pytest.param("C[C@H](C[C@H](C)[2H])[3H]", "(2R,4S)-(4-2H1,2-3H1)pentane", id="mirror_image_keeps_the_isotope_positions"),
        pytest.param("C(=C\\C)/[2H]", "(1E)-(1-2H1)prop-1-ene", id="alkene_stereo_from_a_deuterium"),
        pytest.param("C(=C/C)/[2H]", "(1Z)-(1-2H1)prop-1-ene", id="alkene_stereo_from_a_deuterium_z"),
        pytest.param("[2H]C=CCCOC", "4-methoxy(1-2H1)but-1-ene", id="ether_on_a_labelled_unsaturated_chain"),
        pytest.param("[2H]C#CCOC", "3-methoxy(1-2H)prop-1-yne", id="ether_on_a_labelled_alkyne"),
        pytest.param("[2H]C=CCCOC(C)=O", "(4-2H1)but-3-en-1-yl acetate", id="labelled_unsaturated_alkyl_of_an_ester"),
        pytest.param("C1(=CC=CC=C1)[13C]#N", "benzene(13C)carbonitrile", id="labelled_nitrile_carbon_on_a_ring"),
        pytest.param("C1(=CC=CC=C1)[13C](=O)O", "benzene(13C)carboxylic acid", id="labelled_carboxyl_carbon_on_a_ring"),
        pytest.param("[2H]C(=O)c1ccccc1", "benzene(2H1)carbaldehyde", id="formyl_hydrogen_of_a_ring_aldehyde"),
        pytest.param("[2H]C(=O)c1ccncc1", "pyridine-4-(2H1)carbaldehyde", id="formyl_hydrogen_of_a_heteroring_aldehyde"),
        pytest.param("[2H]C(=O)c1ccccc1C(=O)O", "2-(2H)formylbenzoic acid", id="formyl_hydrogen_of_a_prefix_beside_an_acid"),
        pytest.param("C1(=CC=CC=C1)[13CH2]NN", "[phenyl(13C)methyl]hydrazine", id="labelled_benzyl_carbon_on_a_hydrazine"),
        pytest.param("[14C](=O)(O)C1(CCCCC1)[14C](=O)O", "cyclohexane-1,1-di[(14C)carboxylic acid]", id="identical_labelled_acids_multiplied"),
        pytest.param("[13C](=O)(O)C1(CCCCC1)[14C](=O)O", "1-(13C)carboxycyclohexane-1-(14C)carboxylic acid", id="different_labelled_acids_senior_nuclide_is_suffix"),
        pytest.param("OC(=O)C1(CCCCC1)[13C](=O)O[2H]", "1-carboxycyclohexane-1-(13C,2H)carboxylic acid", id="acid_with_more_modifications_is_suffix"),
        pytest.param("[2H]OC(=O)C1(CCCCC1)[13C](=O)O", "1-(2H)carboxycyclohexane-1-(13C)carboxylic acid", id="labelled_carboxy_prefix"),
        pytest.param("C(OC([2H])([2H])SCOO)[2H]", "{[(2H1)methoxy(2H2)methyl]sulfanyl}methaneperoxol", id="labelled_ether_in_a_sulfanyl_peroxol"),
        pytest.param("OC[C@]1([2H])CC=C[C@@H]1CC(=O)O", "[(1S,5R)-5-(hydroxymethyl)(5-2H)cyclopent-2-en-1-yl]acetic acid", id="labelled_stereocentre_of_a_cycloalkenyl_group"),
        pytest.param("OC[C@]1([2H])CCCC[C@@H]1CC(=O)O", "[(1R,2R)-2-(hydroxymethyl)(2-2H)cyclohexyl]acetic acid", id="labelled_stereocentre_of_a_group_named_in_the_second_candidate_parent"),
    ],
)
def test_isotopes_with_unsaturation_and_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CCC(=O)[18O]CC", "18O-ethyl propan(18O1)oate", id="labelled_ester_oxygen"),
        pytest.param("CCC(=[18O])OCC", "O-ethyl propan(18O1)oate", id="labelled_carbonyl_oxygen_of_an_ester"),
        pytest.param("COC(=O)O[2H]", "methyl (2H)hydrogen carbonate", id="labelled_hydrogen_of_a_carbonate_half_ester"),
        pytest.param("COC(=O)[18OH]", "O-methyl hydrogen (18O1)carbonate", id="labelled_hydroxy_oxygen_of_a_carbonate_half_ester"),
        pytest.param("CC[18O]C(=O)OC", "18O-ethyl O-methyl (18O1)carbonate", id="labelled_bridging_oxygen_of_a_carbonate"),
        pytest.param("CCOC(=[18O])[18O]C", "O-ethyl 18O-methyl (18O2)carbonate", id="two_labelled_oxygens_of_a_carbonate"),
        pytest.param("CCOC(=[18O])OCC", "O,O-diethyl (18O1)carbonate", id="labelled_carbonyl_oxygen_of_a_symmetric_carbonate"),
        pytest.param("CC(=[18O])OC", "O-methyl acet(18O1)ate", id="labelled_acetate_carbonyl"),
        pytest.param("C[18O]C(=O)[13CH3]", "18O-methyl (2-13C)acet(18O1)ate", id="labelled_ester_oxygen_beside_a_labelled_acyl_carbon"),
        pytest.param("CC[C@H](C)C(=O)[18O]C", "18O-methyl (2S)-2-methylbutan(18O1)oate", id="labelled_ester_oxygen_beside_a_stereocentre"),
        pytest.param("COC(=[18O])OC(=O)OC", "dimethyl (O1-18O)dicarbonate", id="carbonyl_oxygen_of_a_dicarbonate"),
        pytest.param("COC(=O)OC(=O)[18O]C", "18O-methyl O-methyl (18O1)dicarbonate", id="ester_oxygen_of_a_dicarbonate"),
        pytest.param("COC(=O)OC(=[18O])O", "3-methyl hydrogen (O1-18O)dicarbonate", id="nuclide_decides_the_locant_of_the_ester"),
        pytest.param("COC(=[18O])OC(N)=O", "methyl amido(O3-18O)dicarbonate", id="amido_end_takes_locant_one"),
        pytest.param("NC(=O)[18O]C(N)=O", "(2-18O)dicarbonic diamide", id="bridging_oxygen_of_a_dicarbonic_diamide"),
        pytest.param("OC(=[18O])OC(=O)O", "(O1-18O)dicarbonic acid", id="oxygen_of_a_dicarbonic_acid"),
        pytest.param("O=C([18O]C)c1ccccc1", "18O-methyl benz(18O1)oate", id="labelled_benzoate_ester_oxygen"),
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
    "smiles,expected",
    [
        pytest.param("CC(=[18O])OC(=O)C", "(O-18O)acetic acetic anhydride", id="carbonyl_oxygen_of_one_acyl_group"),
        pytest.param("CC(=O)[18O]C(=O)C", "acetic (18O)anhydride", id="bridging_oxygen_cited_before_the_class_word"),
        pytest.param("CCC(=O)OC(=[18O])c1ccccc1", "(O-18O)benzoic propanoic anhydride", id="modified_acid_ordered_by_its_letters"),
        pytest.param("CC(=O)OC(N)=[18O]", "acetic (O-18O)carbamic anhydride", id="carbonyl_oxygen_of_a_carbamic_component"),
        pytest.param("CC(=O)OC(=[18O])OC", "acetic (O-methyl hydrogen (18O1)carbonate) anhydride", id="carbonyl_oxygen_of_a_carbonic_component"),
        pytest.param("CC(=O)OC(=O)[18O]C", "acetic (18O-methyl hydrogen (18O1)carbonate) anhydride", id="ester_oxygen_of_a_carbonic_component"),
        pytest.param("NC(=[18O])O", "(O-18O)carbamic acid", id="labelled_carbamic_acid"),
    ],
)
def test_isotope_on_anhydride_oxygens(smiles, expected):
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
        ("[2H]c1cccc(Cc2ccccc2[2H])c1", "1-[(3-2H)phenylmethyl](2-2H)benzene"),
        ("c1ccc(Cc2cccc([2H])c2)c([2H])c1", "1-[(3-2H)phenylmethyl](2-2H)benzene"),
        ("[2H]c1cccc(Cc2ccc[13cH]c2[2H])[13cH]1", "1-[(3-13C,2-2H)phenylmethyl](2-13C,3-2H)benzene"),
        ("c1c[14cH][13cH]c(Cc2ccc[13cH][14cH]2)c1", "1-[(2-13C,3-14C)phenylmethyl](3-13C,2-14C)benzene"),
        ("c1ccc(Cc2ccc[14cH][13cH]2)[14cH][13cH]1", "1-[(2-13C,3-14C)phenylmethyl](3-13C,2-14C)benzene"),
        ("[2H]c1cc(Cc2cccc([2H])n2)ccn1", "4-{[(6-2H)pyridin-2-yl]methyl}(2-2H)pyridine"),
    ],
)
def test_ring_with_lower_locants_for_the_nuclides_is_the_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[2H]c1ccc(cc1)-c1ccccc1", "(4-2H)-1,1'-biphenyl"),
        ("c1ccc(cc1)-c1ccc(cc1)[2H]", "(4-2H)-1,1'-biphenyl"),
        ("c1ccc(cc1)-c1ccccc1[2H]", "(2-2H)-1,1'-biphenyl"),
        ("[2H]c1ccc(cc1)-c1ccc([2H])cc1", "(4,4'-2H2)-1,1'-biphenyl"),
        ("OC(=O)c1ccc(cc1)-c1ccccc1[2H]", "(2'-2H)[1,1'-biphenyl]-4-carboxylic acid"),
        ("OC(=O)Cc1cccc(-c2ccc([2H])cc2)c1", "[(4'-2H)[1,1'-biphenyl]-3-yl]acetic acid"),
        ("OC(=O)Cc1ccccc1-c1c(Cl)cccc1[2H]", "[2'-chloro(6'-2H)[1,1'-biphenyl]-2-yl]acetic acid"),
        ("[2H]c1ccc2ccccc2c1-c1ccc2ccccc2c1", "(2-2H)-1,2'-binaphthalene"),
        ("C(c1ccc2ccccc2c1)c1ccc2cc([2H])ccc2c1", "2-[(naphthalen-2-yl)methyl](6-2H)naphthalene"),
        ("[2H]c1cc2ccccc2cc1Cc1cccc2cccc([2H])c12", "2-{[(8-2H)naphthalen-1-yl]methyl}(3-2H)naphthalene"),
    ],
)
def test_isotopic_modification_of_ring_assemblies_and_fused_parent_choice(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)Cc1ccccc1-c1ccccc1-c1ccccc1[2H]", "[(3²-2H)[1¹,2¹:2²,3¹-terphenyl]-1²-yl]acetic acid", id="chain_of_rings"),
        pytest.param("OC(=O)Cc1ccc2ccccc2c1-c1ccc2ccccc2c1[2H]", "[(1'-2H)[1,2'-binaphthalen]-2-yl]acetic acid", id="fused_components"),
    ],
)
def test_isotopic_modification_of_chained_and_fused_ring_assembly_groups(smiles, expected):
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
        ("CC[18O]c1ccccc1", "(18O)ethoxybenzene"),
        ("CCC[18O]CC", "1-(18O)ethoxypropane"),
        ("CC[18O]C=C", "(18O)ethoxyethene"),
        ("CCC[18O]C=C", "1-[(18O)ethenyloxy]propane"),
        ("[2H]C1CCCCC1(C)CC1CCCCC1", "1-(cyclohexylmethyl)-1-methyl(2-2H1)cyclohexane"),
        ("[14CH2]1CCCCC1CC1([2H])CCCCC1", "1-[(1-2H)cyclohexylmethyl](2-14C)cyclohexane"),
    ],
)
def test_substituent_chains_and_prefixes_are_chosen_by_isotopic_modification(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_deuterium_alkene_stereodescriptor_writes_no_rdkit_logs_to_stderr(capfd):
    assert smiles_to_iupac("C(=C\\C)/[2H]") == "(1E)-(1-2H1)prop-1-ene"
    assert capfd.readouterr().err == ""


# P-82.2.1, P-82.3.1, P-82.5.2, P-82.6.1.3: nuclides on parent hydrides other than hydrocarbons
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]N([2H])[2H]", "(2H3)azane", id="every_hydrogen_of_a_mononuclear_hydride"),
        pytest.param("[18OH2]", "(18O)oxidane", id="skeletal_nuclide_of_a_mononuclear_hydride"),
        pytest.param("[15NH2][2H]", "(2H1,15N)azane", id="nuclide_symbols_alphabetical"),
        pytest.param("[2H]N([2H])N([2H])[2H]", "(2H4)hydrazine", id="chain_hydride_modified_everywhere"),
        pytest.param("[2H]N([2H])N", "(1,1-2H2)hydrazine", id="chain_hydride_modified_in_part"),
        pytest.param("[15NH2]CC[13CH2][13CH2][15NH2]", "(1,2-13C2,15N2)butane-1,4-diamine", id="two_elements_one_in_the_amine_groups"),
        pytest.param("[2H]OOC(=C)C", "prop-1-ene-2-(2H)peroxol", id="nuclide_on_a_peroxol_of_an_unsaturated_chain"),
        pytest.param("CC[18O]O", "ethane(1-18O)peroxol", id="oxygen_bonded_to_the_parent_of_a_peroxol"),
        pytest.param("CO[18OH]", "methane(2-18O)peroxol", id="hydroxy_oxygen_of_a_peroxol_cited_with_its_position"),
        pytest.param("CC(C)S[2H]", "propane-2-(2H)thiol", id="nuclide_on_the_sulfur_of_a_thiol"),
        pytest.param("[2H][N+]([2H])([2H])[2H]", "(2H4)azanium", id="ion_of_a_mononuclear_hydride"),
        pytest.param("[13NH2]", "(13N)azanyl", id="radical_of_a_mononuclear_hydride"),
    ],
)
def test_nuclides_on_other_parent_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.5.1, P-82.5.2, P-82.6.1.1, P-82.6.3.3: nuclides on the nitrogens of C=N, N=N and azoxy groups
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("c1ccc2cc(ccc2c1)[15N]=Nc1ccccc1", "1-(naphthalen-2-yl)-2-phenyl(1-15N)diazene", id="diazene_substituents_fix_the_locants"),
        pytest.param("c1ccc2cc(ccc2c1)N=[15N]c1ccccc1", "1-(naphthalen-2-yl)-2-phenyl(2-15N)diazene", id="diazene_nuclide_on_the_second_nitrogen"),
        pytest.param("c1ccccc1N=[15N]c1ccccc1", "1,2-diphenyl(1-15N)diazene", id="identical_substituents_give_the_nuclide_locant_1"),
        pytest.param("CN=N[2H]", "1-methyl(2-2H)diazene", id="deuterium_on_a_diazene_nitrogen"),
        pytest.param("c1ccc([2H])cc1N=[15N]c1ccccc1", "1-(3-2H)phenyl-2-phenyl(2-15N)diazene", id="nuclides_on_the_substituent_and_the_parent"),
        pytest.param("Clc1c(ccc2ccccc12)[15N]=[N+]([O-])c1ccccc1", "1-(1-chloronaphthalen-2-yl)-2-phenyl(1-15N)diazene 2-oxide", id="azoxy"),
        pytest.param("c1ccccc1N=[15N+]([O-])c1ccccc1", "1,2-diphenyl(1-15N)diazene 1-oxide", id="identical_azoxy_groups_give_the_oxide_locant_1"),
        pytest.param("CCC=[15N]N", "1-propylidene(1-15N)hydrazine", id="hydrazone_nitrogen_of_the_ylidene"),
        pytest.param("CCC=N[15NH2]", "1-propylidene(2-15N)hydrazine", id="hydrazone_terminal_nitrogen"),
        pytest.param("C/C=N/[15NH2]", "(E)-1-ethylidene(2-15N)hydrazine", id="hydrazone_configuration"),
        pytest.param("C=[15N]C", "N-methyl(15N)methanimine", id="imine"),
        pytest.param("CC=N[2H]", "(N-2H)ethanimine", id="deuterated_imine"),
        pytest.param("CC=[15N]O", "N-hydroxy(15N)ethanimine", id="oxime"),
    ],
)
def test_nuclides_on_nitrogens_of_double_bonded_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.2.4, P-82.6.1.1, P-82.6.1.2: the nuclide descriptor of an acyl or sulfonyl halide is that of its acid
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[13C](=O)Cl", "(1-13C)acetyl chloride", id="nuclide_on_the_carbon_of_the_halide_group"),
        pytest.param("ClC(=O)CC[13CH2]Cl", "4-chloro(4-13C)butanoyl chloride", id="nuclide_on_the_chain_of_a_substituted_acyl_halide"),
        pytest.param("[2H]C([2H])([2H])C(=O)Cl", "(2H3)acetyl chloride", id="deuterium_on_the_only_hydrogen_bearing_atom_needs_no_locant"),
        pytest.param("[2H]C([2H])CC(=O)Cl", "(3,3-2H2)propanoyl chloride", id="deuterium_beside_another_hydrogen_bearing_atom_keeps_its_locants"),
        pytest.param("O=C(Cl)[13CH2]C(=O)Cl", "(2-13C)propanedioyl dichloride", id="several_halide_groups"),
        pytest.param("[2H][C@H](C(=O)Cl)CCl", "(2S)-3-chloro(2-2H1)propanoyl chloride", id="descriptor_ranks_the_acyl_halide_carbon_not_the_acid"),
        pytest.param("[2H][C@H](C(=O)F)CCl", "(2R)-3-chloro(2-2H1)propanoyl fluoride", id="halogen_of_the_acyl_halide_changes_the_rank"),
        pytest.param("[13CH3]S(=O)(=O)Cl", "(13C)methanesulfonyl chloride", id="sulfonyl_halide"),
        pytest.param("ClC(=O)c1ccccc1[2H]", "(2-2H)benzoyl chloride", id="aroyl_halide"),
    ],
)
def test_nuclides_on_acyl_and_sulfonyl_halides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.5, P-82.6.1.2: the cyano group has one carbon and one nitrogen, so its nuclides take no locant
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("N#[13C]CCC(=O)O", "3-(13C)cyanopropanoic acid", id="cyano_carbon_before_the_substituent_locant"),
        pytest.param("N#[13C]CC(=O)OC", "methyl (13C)cyanoacetate", id="cyano_carbon_of_an_ester"),
        pytest.param("[15N]#CCC(=O)O", "(15N)cyanoacetic acid", id="cyano_nitrogen"),
    ],
)
def test_nuclides_on_the_cyano_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.6.3.2: the carbonyl carbon of a retained or carbonyl-suffixed acyl group is not numbered, so the nuclide is cited with the systematic name
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OP(O)(=O)[14C](=O)c1ccccc1", "[benzene(14C)carbonyl]phosphonic acid", id="retained_benzoyl_becomes_benzenecarbonyl"),
        pytest.param("OP(O)(=O)[14C](=O)c1ccc2ccccc2c1", "[naphthalene-2-(14C)carbonyl]phosphonic acid", id="carbonyl_with_a_ring_locant"),
    ],
)
def test_nuclide_on_the_carbonyl_carbon_of_an_acyl_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-65.6.3.2.1, P-41: a labelled halogen in the organyl group of a sulfonic, sulfinic or phosphinic ester is cited in that
# group, and the ester stays the parent; the same group on a sulfonyl prefix cites it in the acid name it comes from
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[18F]CCCOS(=O)(=O)C", "3-(18F)fluoropropyl methanesulfonate", id="sulfonate_ester"),
        pytest.param("[18F]CCCOS(C)=O", "3-(18F)fluoropropyl methanesulfinate", id="sulfinate_ester"),
        pytest.param("[18F]CCCOS(=O)(=O)c1ccc(Cl)cc1", "3-(18F)fluoropropyl 4-chlorobenzene-1-sulfonate", id="substituted_arenesulfonate_ester"),
        pytest.param("[18F]CCCS(=O)(=O)c1ccccc1", "[3-(18F)fluoropropane-1-sulfonyl]benzene", id="labelled_group_of_a_sulfonyl_prefix"),
    ],
)
def test_nuclide_in_the_organyl_group_of_an_oxoacid_ester_or_acyl_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.4.1, P-82.4.2: steroid parent hydrides take the locants of the steroid numbering; a nuclide that makes a centre stereogenic adds its descriptor
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[C@@]12[C@@H]([C@H](C)CC[C@@H]([2H])C(C)C)CC[C@H]1[C@@H]1CC[C@H]3CCCC[C@]3(C)[C@H]1CC2", "(24R)-5α-(24-2H1)cholestane", id="side_chain_centre_made_stereogenic_by_deuterium"),
        pytest.param("C[C@@]12[C@@H]([C@H](C)CC[C@H]([2H])C(C)C)CC[C@H]1[C@@H]1CC[C@H]3CCCC[C@]3(C)[C@H]1CC2", "(24S)-5α-(24-2H1)cholestane", id="other_configuration_of_the_same_centre"),
        pytest.param("[2H][C@]1(CC)CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@@]21C", "5α-(17-2H)pregnane", id="existing_centre_needs_no_descriptor"),
        pytest.param("C[C@H](CC[C@@H]([2H])C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C", "(24R)-5α-(24-2H1)cholestan-3β-ol", id="hydroxy_suffix_after_the_configuration"),
        pytest.param("C[C@H](CC[C@@H]([2H])C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", "(24R)-(24-2H1)cholest-5-en-3β-ol", id="unsaturated_parent_with_a_hydroxy_suffix"),
        pytest.param("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[13CH2]C[C@@]43C)[C@@H]1CC[C@@H]2O", "17β-hydroxy(2-13C)androst-4-en-3-one", id="nuclide_after_a_prefix_of_the_parent"),
        pytest.param("C[C@]12CC[C@H](O)C[C@H]1CC[13C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12", "3β-hydroxy-5β-(8-13C)androstan-17-one", id="cis_fused_5beta_parent_keeps_the_face_of_a_ring_a_group"),
        pytest.param("C[C@]12CC[C@@H](O)C[C@H]1CC[13C@@H]1[C@@H]2CC[C@]2(C)[C@@H](O)CC[C@@H]12", "5β-(8-13C)androstane-3α,17β-diol", id="5beta_parent_with_two_hydroxy_groups"),
    ],
)
def test_isotopically_modified_steroid_parent_hydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.3: the descriptor follows the hydro prefixes; a hydrogen added by a hydro prefix beyond those of the mancude parent is the modified one
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[15NH]1CCc2ccccc12", "2,3-dihydro(15N)-1H-indole", id="skeletal_atom_after_the_hydro_prefix"),
        pytest.param("[15NH]1C(C(c2ccccc12)[2H])[2H]", "2,3-dihydro(2,3-2H2,15N)-1H-indole", id="one_deuterium_per_hydro_position_stays_in_the_descriptor"),
        pytest.param("C1CCCc2cc[13cH]cc12", "1,2,3,4-tetrahydro(6-13C)naphthalene", id="descriptor_directly_before_the_stem"),
        pytest.param("[15NH]1C([2H])([2H])C([2H])([2H])c2ccccc12", "2,3-di[(2H)hydro](2,3-2H2,15N)-1H-indole", id="modified_hydro_prefixes"),
        pytest.param("OC1C([2H])([2H])C([2H])([2H])Cc2cc(C)ccc12", "6-methyl-2,3-di[(2H)hydro]-1,4-dihydro(2,3-2H2)naphthalen-1-ol", id="modified_before_unmodified_hydro_prefixes"),
        pytest.param("C1CCC2C([2H])CCCC2C1", "decahydro(1-2H)naphthalene", id="hydro_prefix_without_locants"),
        pytest.param("[2H]C1CCC(=O)c2ccccc21", "3,4-dihydro(4-2H)naphthalen-1(2H)-one", id="added_hydrogen_is_not_a_nuclide"),
        pytest.param("[2H]c1ccc2c(c1)CCCN2C(C)=O", "1-[3,4-dihydro(6-2H)quinolin-1(2H)-yl]ethan-1-one", id="descriptor_of_a_ring_substituent_follows_its_hydro_prefixes"),
        pytest.param("C[13C]12CCCCC1CCCC2", "4a-methyldecahydro(4a-13C)naphthalene", id="skeletal_label_on_a_substituted_ring_fusion_carbon"),
        pytest.param("O[13C]12CCCCC1CCCC2", "octahydro(4a-13C)naphthalen-4a(2H)-ol", id="skeletal_label_on_a_fusion_carbon_that_takes_added_hydrogen"),
        pytest.param("C[13C]12CCCC1CCCC2", "3a-methyloctahydro(3a-13C)-1H-indene", id="substituted_fusion_carbon_of_a_five_membered_ring"),
    ],
)
def test_isotopic_descriptor_after_the_hydro_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1: the nuclide of an ether oxygen is cited in the alkoxy prefix, in front of the part that holds the oxygen: the stem of a
# contracted alkoxy, the 'oxy' of a prefix written with enclosing marks
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("ClCC[18O]CCC", "1-[2-chloro(18O)ethoxy]propane", id="substituted_alkoxy"),
        pytest.param("ClCC[18O]CC", "1-chloro-2-(18O)ethoxyethane", id="alkoxy_that_is_not_the_prefix_of_the_larger_chain"),
        pytest.param("CCC[18O]C(C)C", "1-[(propan-2-yl)(18O)oxy]propane", id="branched_alkoxy"),
        pytest.param("CC(C)[18O]C(C)C", "2-[(propan-2-yl)(18O)oxy]propane", id="two_equal_branched_groups"),
        pytest.param("COCC[18O]C=C", "[2-methoxy(18O)ethoxy]ethene", id="alkoxy_with_an_alkoxy_prefix_on_a_vinyl_ether"),
        pytest.param("ClCC[18O]C=C", "[2-chloro(18O)ethoxy]ethene", id="substituted_alkoxy_on_a_vinyl_ether"),
        pytest.param("CCC[18O]C(C)=C", "2-(18O)propoxyprop-1-ene", id="plain_alkoxy_on_a_branched_alkene"),
        pytest.param("CCC[18O]COC", "1-[methoxy(18O)methoxy]propane", id="oxygen_of_the_outer_of_two_alkoxy_groups"),
    ],
)
def test_nuclide_on_the_oxygen_of_a_substituted_or_branched_alkoxy_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.2: hydrogen nuclides on the nitrogens of a hydrazinecarboxamide take N (amide) and 1, 2 (hydrazine) and stand in
# front of 'hydrazine'
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]NC(=O)NN", "(N-2H1)hydrazinecarboxamide", id="amide_nitrogen"),
        pytest.param("[2H]NNC(N)=O", "(2-2H1)hydrazinecarboxamide", id="terminal_hydrazine_nitrogen"),
        pytest.param("[2H]N(N)C(N)=O", "(1-2H)hydrazinecarboxamide", id="substituted_hydrazine_nitrogen"),
        pytest.param("NC(=O)NN([2H])[2H]", "(2,2-2H2)hydrazinecarboxamide", id="two_deuteriums_on_one_nitrogen"),
        pytest.param("[2H]NC(=O)NN[2H]", "(2,N-2H2)hydrazinecarboxamide", id="locants_cited_numbers_before_letters"),
        pytest.param("[2H]NC(=O)NNC", "2-methyl(N-2H1)hydrazine-1-carboxamide", id="descriptor_follows_the_prefixes"),
        pytest.param("[2H]N(C)C(=O)NN", "N-methyl(N-2H)hydrazinecarboxamide", id="amide_nitrogen_with_a_substituent"),
        pytest.param("NC(=O)NNC([2H])([2H])[2H]", "2-(2H3)methylhydrazine-1-carboxamide", id="labelled_substituent_of_the_hydrazine"),
    ],
)
def test_hydrogen_nuclides_on_the_nitrogens_of_a_semicarbazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.2.2.1, P-82.2.4: the acid hydrogen of a partial ester of a noncarbon oxoacid is the word 'hydrogen' that carries the nuclide,
# the modified one cited first; nuclides elsewhere stay on their groups
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]OS(=O)(=O)OC", "methyl (2H)hydrogen sulfate", id="sulfate"),
        pytest.param("[2H]OP(=O)(OC)OC", "dimethyl (2H)hydrogen phosphate", id="phosphate_diester"),
        pytest.param("[2H]OP(=O)(O)OC[C@H](C)O", "(2S)-2-hydroxypropyl (2H)hydrogen hydrogen phosphate", id="phosphate_monoester_beside_an_alcohol"),
        pytest.param("[2H]OP(=O)(O[2H])OC", "methyl (2H2)dihydrogen phosphate", id="both_hydrogens_of_a_monoester"),
        pytest.param("[2H]OP(=O)(OC)c1ccccc1", "methyl (2H)hydrogen phenylphosphonate", id="phosphonate"),
        pytest.param("[2H]OP(=O)(OC([2H])([2H])[2H])OC", "(2H3)methyl methyl (2H)hydrogen phosphate", id="nuclides_on_an_alkyl_group_too"),
        pytest.param("[2H]OS(=O)(=O)Cl", "(O-2H)sulfurochloridic acid", id="sulfurochloridic_acid"),
        pytest.param("[2H]OP(=O)(O)N", "(O-2H)phosphoramidic acid", id="phosphoramidic_acid"),
    ],
)
def test_acid_hydrogen_nuclide_of_a_noncarbon_oxoacid_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.1.1, P-82.6.2, P-44.4.1.11: the carbon of the carboxamide or carbohydrazide suffix stands in front of the suffix (with
# the chalcogen when both are modified), a lone chalcogen in front of 'hydrazine'; the hydrazine unit that is modified is the parent
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("N[13C](=O)NN", "hydrazine(13C)carboxamide", id="suffix_carbon"),
        pytest.param("NC(=[18O])NN", "(18O)hydrazinecarboxamide", id="lone_oxygen"),
        pytest.param("N[13C](=[18O])N[15NH2]", "(2-15N)hydrazine(13C,18O)carboxamide", id="carbon_oxygen_and_nitrogen"),
        pytest.param("NC(=[18O])NN[2H]", "(2-2H1,O-18O)hydrazinecarboxamide", id="oxygen_located_beside_a_located_nuclide"),
        pytest.param("N[13C](=S)NNC", "2-methylhydrazine-1-(13C)carbothioamide", id="suffix_carbon_after_a_locant"),
        pytest.param("NN[13C](=[18O])NN", "hydrazine(13C,18O)carbohydrazide", id="carbohydrazide_carbon_and_oxygen"),
        pytest.param("[15NH2]NC(=O)NN", "(2-15N)hydrazinecarbohydrazide", id="modified_unit_is_the_parent"),
        pytest.param("[2H]NNC(=O)N([2H])N", "(1,N'-2H2)hydrazinecarbohydrazide", id="lower_locant_unit_is_the_parent"),
        pytest.param("CNNC(=O)N[15NH2]", "2-methyl(N'-15N)hydrazine-1-carbohydrazide", id="substituent_decides_before_modification"),
        pytest.param("N#CN[15NH2]", "(N'-15N)cyanohydrazide", id="nitrogen_of_a_cyanohydrazide"),
    ],
)
def test_nuclides_of_the_carboxamide_group_and_of_a_carbonohydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.5, P-82.6.2: the acyl group with the more modified atoms is the hydrazide parent, then the lowest locants (N before N');
# the nuclides of the other acyl group stand in front of its stem
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]CC(=O)NNC(C)=O", "N'-acetyl(2-2H1)acetohydrazide", id="labelled_acyl_is_the_parent"),
        pytest.param("CC(=O)NNC(=O)C[2H]", "N'-acetyl(2-2H1)acetohydrazide", id="labelled_acyl_is_the_parent_in_either_order"),
        pytest.param("[2H]c1ccc(C(=O)NNC=O)cc1", "N'-formyl(4-2H)benzohydrazide", id="different_acyl_groups"),
        pytest.param("[13CH3]C(=O)NNC(C)=O", "N'-acetyl(2-13C)acetohydrazide", id="carbon_nuclide"),
        pytest.param("[2H]N(NC(C)=O)C(C)=O", "N'-acetyl(N-2H)acetohydrazide", id="nitrogen_nuclide_takes_the_lower_locant_N"),
        pytest.param("CC(=O)N([2H])N([2H])C(C)=O", "N'-acetyl(N,N'-2H2)acetohydrazide", id="both_nitrogens"),
        pytest.param("[2H]C([2H])([2H])C(=O)NNC(=O)C([2H])([2H])[2H]", "N'-[(2H3)acetyl](2,2,2-2H3)acetohydrazide", id="both_acyl_groups"),
        pytest.param("[2H]CCC(=O)NNC(=O)CC", "N'-propanoyl(3-2H1)propanehydrazide", id="longer_chain"),
    ],
)
def test_nuclides_on_a_diacylhydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.3.3: the nitrogen of the amide parent takes its nuclide in front of the amide name, the sulfur of the
# sulfonyl group in front of the hydride name inside its prefix
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)[15NH]S(C)(=O)=O", "N-(methanesulfonyl)(15N)acetamide", id="amide_nitrogen_of_an_acylsulfonamide"),
        pytest.param("CC(=O)N[34S](C)(=O)=O", "N-[(34S)methanesulfonyl]acetamide", id="sulfonyl_sulfur_of_an_acylsulfonamide"),
        pytest.param("CC(=O)[15NH]C(C)=O", "N-acetyl(15N)acetamide", id="nitrogen_of_a_diacylamine"),
        pytest.param("OC(=O)CS(=[18O])(=O)Cl", "[chloro(O-18O)sulfonyl]acetic acid", id="oxo_oxygen_of_a_halosulfonyl_prefix"),
        pytest.param("OC(=O)CS(=O)(=O)[15NH]C", "[methyl(N-15N)sulfamoyl]acetic acid", id="nitrogen_of_a_sulfamoyl_prefix"),
        pytest.param("OC(=O)C[18O]C(C)=O", "[acetyl(18O)oxy]acetic acid", id="linking_oxygen_of_an_acyloxy_prefix"),
    ],
)
def test_nuclides_on_an_acylsulfonamide_or_diacylamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.1, P-41: a group that outranks the others is the parent whatever else the chain carries, and the nuclides of a
# chain that is the parent skeleton of a polyol ester stand in front of that skeleton
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]CC(NC(N)=O)CS(=O)(=O)C", "[1-(methanesulfonyl)(3-2H1)propan-2-yl]urea", id="urea_beside_a_sulfone"),
        pytest.param("[2H]CC(NC(=N)N)CS(=O)C", "N-[1-(methanesulfinyl)(3-2H1)propan-2-yl]guanidine", id="guanidine_beside_a_sulfoxide"),
        pytest.param("[2H]CC(OC(=O)C)COC(N)=O", "(3-2H1)propane-1,2-diyl 2-acetate 1-carbamate", id="polyol_skeleton_of_two_esters"),
        pytest.param("[NH3+]C[13CH2]P(=O)([O-])O", "hydrogen [2-azaniumyl(1-13C)ethyl]phosphonate", id="phosphonate_monoanion"),
    ],
)
def test_nuclide_on_a_chain_beside_a_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-68.4.1.1: a nuclide on an end group of a polychalcogen chain is cited in that organyl prefix
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]CSSSC", "methyl[(2H1)methyl]trisulfane", id="trisulfane"),
        pytest.param("[13CH3]SSSC", "methyl[(13C)methyl]trisulfane", id="carbon_nuclide"),
        pytest.param("[2H]CSSSSC", "methyl[(2H1)methyl]tetrasulfane", id="tetrasulfane"),
        pytest.param("[2H]C([2H])([2H])SSSC([2H])([2H])[2H]", "di[(2H3)methyl]trisulfane", id="both_end_groups"),
        pytest.param("C[Se][Se][Se]C[2H]", "methyl[(2H1)methyl]triselane", id="triselane"),
        pytest.param("[2H]CCSSSC", "(2-2H1)ethyl(methyl)trisulfane", id="different_end_groups"),
    ],
)
def test_nuclides_on_the_end_groups_of_a_polychalcogen_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-73.2.2.3: the diazonium suffix stands where the amine suffix does, so the nuclides of the skeleton are cited as for an amine
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]C[N+]#N", "(2H1)methanediazonium", id="methane"),
        pytest.param("[13CH3][N+]#N", "(13C)methanediazonium", id="carbon_nuclide"),
        pytest.param("[2H]CC[N+]#N", "(2-2H1)ethane-1-diazonium", id="locant_cited_with_a_nuclide"),
        pytest.param("[2H]C(C)(C)[N+]#N", "(2-2H)propane-2-diazonium", id="branched_chain"),
        pytest.param("[2H]C=CC[N+]#N", "(3-2H1)prop-2-ene-1-diazonium", id="unsaturated_chain"),
        pytest.param("[2H]C1CCCCC1[N+]#N", "(2-2H1)cyclohexane-1-diazonium", id="ring"),
        pytest.param("[2H]c1ccccc1[N+]#N", "(2-2H)benzenediazonium", id="benzene_ring"),
        pytest.param("ClC([2H])[N+]#N", "1-chloro(2H1)methanediazonium", id="halogen_beside_the_nuclide"),
    ],
)
def test_nuclides_on_the_skeleton_of_a_diazonium_cation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.1.1, P-82.6.3.3: the nuclide of a heteroatom stands in front of the part of the name that denotes that atom
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[34S]C", "[methyl(34S)sulfanyl]methane", id="sulfide_linking_two_carbon_groups"),
        pytest.param("OCC[15N+](=O)[O-]", "2-(15N)nitroethan-1-ol", id="nitro_prefix"),
        pytest.param("CC(=[34S])C", "propane-2-(34S)thione", id="thione_suffix"),
        pytest.param("C[32P](=O)(O)O", "methyl(32P)phosphonic acid", id="phosphonic_acid"),
        pytest.param("C[32P](=O)(C)O", "dimethyl(32P)phosphinic acid", id="phosphinic_acid"),
        pytest.param("C[32P](=O)(Cl)Cl", "methyl(32P)phosphonic dichloride", id="phosphorus_acid_halide"),
        pytest.param("CP(=[18O])(C)C", "trimethyl-λ5-phosphan(18O)one", id="phosphanone_oxygen"),
        pytest.param("C[32P](=O)(C)C", "trimethyl(32P)-λ5-phosphanone", id="phosphanone_centre"),
        pytest.param("CN(C)[15N]=O", "dimethyl(15N)nitrous amide", id="nitrous_amide"),
        pytest.param("CS(=[34S])(=O)C", "1,1-dimethyl(2-34S)-1λ6-disulfen-1-one", id="chalcogen_chain_atom"),
        pytest.param("C[31P](C)C", "trimethyl(31P)phosphane", id="phosphane"),
        pytest.param("ClC[31P](C)C", "(chloromethyl)di(methyl)(31P)phosphane", id="mononuclear_parent_beside_a_halogen"),
        pytest.param("OC[29SiH2]C", "[methyl(29Si)silyl]methanol", id="hydride_atom_as_a_prefix"),
        pytest.param("C[SiH2][29SiH2]C", "1,2-dimethyl(1-29Si)disilane", id="hydride_chain_parent"),
        pytest.param("O=[32P](O)(O)O", "(32P)phosphoric acid", id="mononuclear_oxoacid"),
        pytest.param("CO[32P](=O)(OC)OC", "trimethyl (32P)phosphate", id="phosphate_ester"),
        pytest.param("CO[34S](=O)(=O)OC", "dimethyl (34S)sulfate", id="sulfate"),
        pytest.param("C[34S](=O)(=O)Cl", "(34S)methanesulfonyl chloride", id="sulfonyl_chloride"),
        pytest.param("C[34S](=O)(=O)N", "(34S)methanesulfonamide", id="sulfonamide_sulfur"),
        pytest.param("C[11B](O)O", "methyl(11B)boronic acid", id="boronic_acid"),
        pytest.param("O=C(O)C[31P](=O)(O)O", "(31P)phosphonoacetic acid", id="phosphono_prefix"),
        pytest.param("N[34S](N)(=O)=O", "(34S)sulfuric diamide", id="parent_stem_of_a_lone_atom"),
        pytest.param("OB(O)[10B]", "2,2-dihydroxy(1-10B)diboran(4)-1-ylidene", id="boron_of_a_radical_diborane_chain"),
    ],
)
def test_nuclide_on_a_heteroatom_of_the_parent_or_a_linking_atom(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.4.1, P-92.1.4: a carbon made stereogenic only by its nuclides cites its descriptor in front of the group that holds it
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H][C@@H](C)NC(N)=O", "[(1S)-(1-2H1)ethyl]urea", id="urea"),
        pytest.param("[2H][C@@H](C)NC(=O)OC", "methyl [(1S)-(1-2H1)ethyl]carbamate", id="carbamate_ester"),
        pytest.param("[2H][C@H](C)OC(=O)OC", "(1R)-(1-2H1)ethyl methyl carbonate", id="carbonate_ester_cited_in_alphanumerical_order"),
        pytest.param("[2H][C@H](C)B(O)O", "[(1S)-(1-2H1)ethyl]boronic acid", id="boron_ranks_below_carbon"),
    ],
)
def test_stereocentre_made_by_nuclides_beside_a_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.2: a labelled organyl group of an acid is cited in front of its own name, whichever dedicated namer names the acid
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[13CH3]CNC(=O)O", "(2-13C)ethylcarbamic acid", id="carbamic_acid"),
        pytest.param("[13CH3]CCB(O)O", "(3-13C)propylboronic acid", id="boronic_acid"),
        pytest.param("[13CH3]C([N+](=O)[O-])CB(O)O", "[2-nitro(3-13C)propyl]boronic acid", id="boronic_acid_with_a_nitro_prefix"),
        pytest.param("[13CH3]CP(=O)(O)Cl", "(2-13C)ethylphosphonochloridic acid", id="functional_replacement_phosphonic_acid"),
    ],
)
def test_labelled_organyl_group_of_an_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.5.1: a nuclide on an acid centre of a chain takes the chain locant; on a sulfur imide it stands in front of the imino prefix or suffix
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=P(O)(O)O[31P](=O)(O)O", "(1-31P)diphosphoric acid", id="polynuclear_acid_centre"),
        pytest.param("CS(C)(=[15NH])=O", "(15N)iminodi(methyl)-λ6-sulfanone", id="imino_prefix_of_a_sulfoximide"),
        pytest.param("CCS(CC)(=N)=[18O]", "diethyl(imino)-λ6-sulfan(18O)one", id="one_suffix_of_a_sulfoximide"),
        pytest.param("CS(C)=[15NH]", "dimethyl-λ4-sulfan(15N)imine", id="imine_suffix_of_a_sulfimide"),
    ],
)
def test_nuclide_on_a_polynuclear_acid_centre_or_a_sulfur_imide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.3.2, P-66.1.2: the junior acyl group of a diacylamine is an N-prefix that cites the nuclides of its own carbonyl carbon
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[13C](=O)NC(=O)c1ccccc1", "N-[(1-13C)acetyl]benzamide", id="carbonyl_carbon_of_acetyl"),
        pytest.param("ClC[13C](=O)NC(=O)c1ccccc1", "N-[chloro(1-13C)acetyl]benzamide", id="substituted_acetyl"),
        pytest.param("NC[13C](=O)NC(=O)c1ccccc1", "N-[amino(1-13C)acetyl]benzamide", id="amino_acid_acyl_without_its_retained_name"),
        pytest.param("C[13C](=O)N(C(C)=O)C(=O)c1ccccc1", "N-[(1-13C)acetyl]-N-acetylbenzamide", id="labelled_acyl_beside_an_identical_unlabelled_one"),
        pytest.param("CC[13C](=O)NC(=O)c1ccccc1", "N-[(1-13C)propanoyl]benzamide", id="acyl_group_other_than_acetyl"),
        pytest.param("C[13C](=[18O])NNC(=O)c1ccccc1", "N'-[(1-13C,18O)acetyl]benzohydrazide", id="carbonyl_oxygen_with_the_carbon_on_a_hydrazide"),
        pytest.param("CC(=[18O])NC(=O)c1ccccc1", "N-(18O)acetylbenzamide", id="carbonyl_oxygen_of_the_acyl_prefix"),
        pytest.param("C[13C](=O)NC(=O)NC(C)=O", "N-(acetylcarbamoyl)(1-13C)acetamide", id="acyl_inside_a_carbamoyl_prefix"),
        pytest.param("C[13C](=O)NC(=O)NC(=O)c1ccccc1", "N-{[(1-13C)acetyl]carbamoyl}benzamide", id="labelled_acyl_inside_the_carbamoyl_prefix_of_a_senior_amide"),
        pytest.param("C[13C](=O)NNC(=O)[13c]1ccccc1", "N'-[(1-13C)acetyl](1-13C)benzohydrazide", id="labelled_parent_and_labelled_acyl_prefix"),
    ],
)
def test_nuclide_on_the_carbonyl_carbon_of_an_n_acyl_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.5.1, P-62.5: the nuclide of the oxidized nitrogen stands in front of the amine name as for the neutral amine, the class word follows
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[15N+](C)(C)[O-]", "N,N-dimethyl(15N)methanamine N-oxide", id="acyclic_amine_oxide"),
        pytest.param("C[15N+]1([O-])CCCCC1", "1-methyl(15N)piperidine 1-oxide", id="ring_amine_oxide"),
        pytest.param("[O-][15n+]1ccncc1", "(1-15N)pyrazine 1-oxide", id="oxidized_nitrogen_of_a_diazine"),
        pytest.param("[O-][n+]1cc[15n]cc1", "(4-15N)pyrazine 1-oxide", id="other_nitrogen_of_a_diazine_oxide"),
        pytest.param("C[N+](C)([O-])C[N+](C)([13CH3])[O-]", "N-{[dimethyl(oxo)-λ5-azanyl]methyl}-N-methyl(13C)methanamine N-oxide", id="more_modified_oxide_is_the_parent"),
        pytest.param("C[N+](C)([O-])C[15N+](C)(C)[O-]", "1-[dimethyl(oxo)-λ5-azanyl]-N,N-dimethyl(15N)methanamine N-oxide", id="oxidized_nitrogen_of_the_parent_of_a_dioxide"),
        pytest.param("[O-][N+](C)(C)[13CH2]CN(C)C", "2-(dimethylamino)-N,N-dimethyl(1-13C)ethan-1-amine N-oxide", id="backbone_carbon_beside_an_amino_prefix"),
        pytest.param("[O-][N+](C)(C)CC[15N](C)C", "2-[dimethyl(15N)amino]-N,N-dimethylethan-1-amine N-oxide", id="nitrogen_of_an_amino_prefix"),
        pytest.param("[O-][N+](C)(C)CCN(C)CC[15N](C)C", "2-({2-[dimethyl(15N)amino]ethyl}(methyl)amino)-N,N-dimethylethan-1-amine N-oxide", id="nitrogen_of_an_amino_prefix_inside_another"),
        pytest.param("C[N+]([O-])(C)CC[15N]1CCCC1", "N,N-dimethyl-2-[(15N)pyrrolidin-1-yl]ethan-1-amine N-oxide", id="nitrogen_of_a_ring_amino_prefix"),
    ],
)
def test_nuclide_on_the_nitrogen_of_an_n_oxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-74.1.3: the descriptor of a zwitterion whose anion is the parent stands in front of the parent hydride, a nuclide of the cationic nitrogen in front of its prefix stem
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[2H]C([NH3+])C[O-]", "2-azaniumyl(2-2H1)ethan-1-olate", id="chain_of_an_ammonium_olate"),
        pytest.param("[2H]C([N+](C)(C)C)C[O-]", "2-(N,N-dimethylmethanaminiumyl)(2-2H1)ethan-1-olate", id="chain_of_a_quaternary_ammonium_olate"),
        pytest.param("C[15N+](C)(C)CC[O-]", "2-[trimethyl(15N)azaniumyl]ethan-1-olate", id="cationic_nitrogen_of_the_prefix"),
    ],
)
def test_nuclide_in_a_zwitterion_with_an_olate_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.2: a nuclide on an oxygen of an acid group takes the letter locant O and, with the nuclide of the centre, shares its descriptor
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CP(=O)(O)[18OH]", "methyl(O-18O)phosphonic acid", id="phosphonic_acid"),
        pytest.param("O=[31P](O)(O)[18OH]", "(O-18O,31P)phosphoric acid", id="oxygen_and_centre_in_one_descriptor"),
        pytest.param("OB([18OH])c1ccccc1", "phenyl(O-18O)boronic acid", id="boronic_acid"),
        pytest.param("CP(=[18O])(Cl)Cl", "methyl(O-18O)phosphonic dichloride", id="acid_halide"),
        pytest.param("O=C(O)CP(O)(O)=[18O]", "(O-18O)phosphonoacetic acid", id="phosphono_prefix"),
        pytest.param("CS(N)(=O)=[18O]", "(O-18O)methanesulfonamide", id="sulfonamide"),
        pytest.param("NS(=O)(=O)[18OH]", "(O-18O)sulfamic acid", id="sulfamic_acid"),
    ],
)
def test_nuclide_on_an_oxygen_of_an_acid_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.3.3: a nuclide on an atom of a pseudohalide group stands in front of the part of the name that denotes the atom
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CN=C=[18O]", "iso(18O)cyanatomethane", id="isocyanato_group"),
        pytest.param("C[N+]#[15N]", "methane(15N)diazonium", id="diazonium_group"),
        pytest.param("N#C[18O-]", "(18O)cyanate", id="cyanate_anion"),
    ],
)
def test_nuclide_on_an_atom_of_a_pseudohalide_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.2.4: a nuclide on an oxygen of a peroxy linkage stands in front of 'peroxy' or of the class word
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("ClCO[18O]C", "chloro[methyl(2-18O)peroxy]methane", id="peroxy_prefix"),
        pytest.param("CC(=O)O[18O]C(=O)c1ccccc1", "acetic benzoic (2-18O)peroxyanhydride", id="peroxyanhydride_bridge"),
        pytest.param("CP(=O)(O)O[18O]C", "O18O-methyl hydrogen methylphosphono(perox(18O1)oate)", id="peroxy_ester_of_a_phosphorus_acid"),
    ],
)
def test_nuclide_on_an_oxygen_of_a_peroxy_linkage(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.2.4, P-82.6.4: nuclides on the nitrogen of a carbamate, the oxygens of a nitro or nitrate group and the oxygens of a sulfonate
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[15NH]C(=O)OC", "methyl N-methyl(N-15N)carbamate", id="carbamate_nitrogen"),
        pytest.param("CO[N+](=[18O])[O-]", "methyl (O-18O)nitrate", id="nitrate_oxygen"),
        pytest.param("C[18O]S(=O)(=O)c1ccccc1", "18O-methyl benzene(18O1)sulfonate", id="sulfonate_ester_oxygen"),
    ],
)
def test_nuclide_on_the_nitrogen_or_oxygens_of_a_carbamate_nitrate_and_sulfonate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.4, P-82.2.1: a halogen nuclide of a halide stands in front of the word or the infix that holds it
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C([36Cl])[36Cl]", "carbonyl (36Cl2)dichloride", id="acyl_dihalide"),
        pytest.param("O=P([36Cl])(Cl)Cl", "phosphoryl (36Cl)chloride dichloride", id="oxoacid_halide"),
        pytest.param("COC(=O)[36Cl]", "methyl carbono(36Cl)chloridate", id="halocarbonic_ester"),
    ],
)
def test_nuclide_on_the_halogen_of_an_acid_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.2.1, P-82.6.4, P-62.5: a nuclide on the oxide oxygen of an amine or imine oxide and on one oxygen of oxalic acid
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[N+](C)([18O-])CC", "N,N-dimethylethanamine N-(18O)oxide", id="acyclic_amine_oxide"),
        pytest.param("CC=[N+](C)[18O-]", "N-methylethanimine N-(18O)oxide", id="imine_oxide"),
        pytest.param("OC(=O)C(=O)[18OH]", "(O-18O)oxalic acid", id="oxalic_acid_one_oxygen"),
    ],
)
def test_nuclide_on_the_oxide_oxygen_and_on_oxalic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-82.6.3.3, P-82.2.1, P-82.6.4: nuclides on the nitrogens of an azine, a hydrazide-type group and an azoxy or bis-azo compound
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC=N[15N]=C(C)C", "1-ethylidene-2-(propan-2-ylidene)(2-15N)hydrazine", id="azine"),
        pytest.param("C[N+]([18O-])=NC", "1,2-dimethyldiazene 1-(18O)oxide", id="azoxy_oxide_oxygen"),
        pytest.param("c1ccccc1[15N]=Nc1ccc(cc1)N=Nc1ccccc1", "1-phenyl-2-[4-(phenyldiazenyl)phenyl](1-15N)diazene", id="two_azo_groups"),
    ],
)
def test_nuclide_on_the_nitrogens_of_an_azine_azoxy_and_bis_azo_compound(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
