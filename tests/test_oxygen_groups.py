import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(OCC)OCC", "1,1-diethoxypropane"),
    ],
)
def test_acetal_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_acetal():
    assert smiles_to_iupac("C1CCC(OC)(OC)CC1") == "1,1-dimethoxycyclohexane"


def test_branched_alkoxy():
    assert smiles_to_iupac("CCC(OC(C)C)OCC") == "1-ethoxy-1-(propan-2-yloxy)propane"


def test_unsaturated_acetal():
    assert smiles_to_iupac("C=CC(OCC)OCC") == "3,3-diethoxyprop-1-ene"


def test_hemiacetal_named_as_alkoxy_alcohol():
    assert smiles_to_iupac("OC(OCC)CC") == "1-ethoxypropan-1-ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(OC)OC", "(2S)-1,1-dimethoxy-2-methylbutane"),
    ],
)
def test_acetal_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_acyl_halide():
    assert smiles_to_iupac("C=CC(=O)Cl") == "prop-2-enoyl chloride"


def test_branched_substituent():
    assert smiles_to_iupac("CC(C)CC(=O)Cl") == "3-methylbutanoyl chloride"


def test_branch_stereocenter_cites_the_specified_elements():
    assert smiles_to_iupac("CC(C[C@H](C)Cl)C(=O)Cl") == "(4S)-4-chloro-2-methylpentanoyl chloride"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC(=O)Cl", "phenylacetyl chloride"),
    ],
)
def test_phenyl_chain_acyl_halide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("C[C@H](O)/C=C/C") == "(2S,3E)-pent-3-en-2-ol"


def test_ring_stereocenter_cites_the_specified_elements():
    assert smiles_to_iupac("O[C@H]1CCCCC1C") == "(1S)-2-methylcyclohexan-1-ol"


def test_ring_branch_stereocenter_supported():
    assert smiles_to_iupac("OC1(CCCCC1)[C@@H](C)CC") == "1-[(2S)-butan-2-yl]cyclohexan-1-ol"


def test_ring_branch_two_stereocenters():
    assert (
        smiles_to_iupac("OC1(CCCCC1)[C@@H]([C@H](C)Cl)CC")
        == "1-[(2S,3S)-2-chloropentan-3-yl]cyclohexan-1-ol"
    )


def test_polycyclic_ring_stereocenter_cites_the_specified_elements():
    assert smiles_to_iupac("O[C@H]1CCCC2CCCC12") == "(4S)-octahydro-1H-inden-4-ol"


def test_ring_hydroxyl_count_exceeds_chain():
    assert smiles_to_iupac("OC1C(O)CC(CO)CC1") == "4-(hydroxymethyl)cyclohexane-1,2-diol"


def test_ring_with_hydroxyl_and_branched_diol_chain():
    assert smiles_to_iupac("OC1CCCCC1C(CO)CO") == "2-(2-hydroxycyclohexyl)propane-1,3-diol"


def test_phenyl_chain_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CO") == "(2-ethenylphenyl)methanol"


def test_heteroaromatic_substituent_alcohol_directly_on_ring():
    assert smiles_to_iupac("Oc1cccnc1") == "pyridin-3-ol"


def test_phenol_multiple_hydroxyls():
    assert smiles_to_iupac("Oc1ccccc1O") == "benzene-1,2-diol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1CCCCC1c1ccc[nH]1", "2-(1H-pyrrol-2-yl)cyclohexan-1-ol"),
    ],
)
def test_two_ring_aromatic_substituent_alcohol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_aromatic_substituent_alcohol_chain_hydroxyl():
    assert smiles_to_iupac("OCC1CCCCC1c1ccccc1") == "(2-phenylcyclohexyl)methanol"


def test_ether_on_ring():
    assert smiles_to_iupac("OC1CCCCC1OC") == "2-methoxycyclohexan-1-ol"


def test_enol():
    assert smiles_to_iupac("OC=CC") == "prop-1-en-1-ol"


def test_unsaturated_ring_alcohol_with_substituent():
    assert smiles_to_iupac("OC1CCCC=C1C") == "2-methylcyclohex-2-en-1-ol"


def test_unsaturated_ring_alcohol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCC#C1")


def test_amine_hetero_mix_dispatches_to_alcohol_amine():
    assert smiles_to_iupac("OCCN") == "2-aminoethan-1-ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(c1ccccc1)(c1ccccc1)c1ccccc1", "triphenylmethanol"),
    ],
)
def test_chain_alcohol_with_several_aromatic_ring_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formaldehyde_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=O")


def test_ketene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCCCCCCCCCCC=C=O")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=C=O")


def test_benzaldehyde():
    assert smiles_to_iupac("O=Cc1ccccc1") == "benzaldehyde"
    assert smiles_to_iupac("O=Cc1ccc(C)cc1") == "4-methylbenzaldehyde"
    assert smiles_to_iupac("O=Cc1ccccc1Cl") == "2-chlorobenzaldehyde"


def test_phenyl_chain_aldehyde_with_hydroxyl():
    assert smiles_to_iupac("OCc1ccccc1CC=O") == "2-[2-(hydroxymethyl)phenyl]acetaldehyde"


def test_phenyl_chain_aldehyde_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC=O") == "2-(2-ethenylphenyl)acetaldehyde"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CC1CCCCC1c1ccccc1", "2-phenylcyclohexane-1-carbaldehyde"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_aldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_aromatic_substituent_aldehyde_substituted_ring():
    assert (
        smiles_to_iupac("O=CC1CCCC(c2ccc(C)cc2)C1")
        == "3-(4-methylphenyl)cyclohexane-1-carbaldehyde"
    )


def test_two_ring_aromatic_substituent_aldehyde_chain_aldehyde():
    assert smiles_to_iupac("O=CCC1CCCCC1c1ccccc1") == "2-(2-phenylcyclohexyl)acetaldehyde"


def test_heteroaromatic_ring_directly_attached_aldehyde():
    assert smiles_to_iupac("O=Cc1cccnc1") == "pyridine-3-carbaldehyde"


def test_aldehyde_enol_mix():
    assert smiles_to_iupac("OC=CC=O") == "3-hydroxyprop-2-enal"


def test_ring_aldehyde_multiple_groups():
    assert smiles_to_iupac("O=CC1CCC(C=O)CC1") == "cyclohexane-1,4-dicarbaldehyde"


def test_acyclic_aldehyde_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("O=C[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enal"
    assert smiles_to_iupac("O=C[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enal"


def test_cyclic_aldehyde_ring_stereocenter():
    assert smiles_to_iupac("O=C[C@H]1CCCC[C@@H]1Cl") == "(1R,2S)-2-chlorocyclohexane-1-carbaldehyde"


def test_aldehyde_ring_branch_stereocenter():
    assert (
        smiles_to_iupac("O=C[C@H]1CC[C@H](C[C@@H](C)CC)C1")
        == "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carbaldehyde"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=COC=O", "formic anhydride"),
    ],
)
def test_symmetric_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_symmetric_anhydride():
    assert smiles_to_iupac("CC(C)C(=O)OC(=O)C(C)C") == "2-methylpropanoic anhydride"


def test_halogen_substituent():
    assert smiles_to_iupac("ClCC(=O)OC(=O)CCl") == "chloroacetic anhydride"


def test_cyclic_anhydride_names_as_ring_dione():
    assert smiles_to_iupac("O=C1CCC(=O)O1") == "oxolane-2,5-dione"


def test_carbonic_acid():
    assert smiles_to_iupac("OC(=O)O") == "carbonic acid"


def test_dimethyl_carbonate():
    assert smiles_to_iupac("COC(=O)OC") == "dimethyl carbonate"


def test_methyl_hydrogen_carbonate():
    assert smiles_to_iupac("COC(=O)O") == "methyl hydrogen carbonate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # -oate combined with existing unsaturation support.
        ("C=CC(=O)[O-]", "prop-2-enoate"),
    ],
)
def test_carboxylate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_dicarboxylate_is_named():
    assert smiles_to_iupac("[O-]C(=O)CC(=O)[O-]") == "propanedioate"


def test_ring_carboxylate_is_named():
    assert smiles_to_iupac("[O-]C(=O)C1CCCCC1") == "cyclohexanecarboxylate"


def test_amine_coexisting_is_named():
    assert smiles_to_iupac("NCC(=O)[O-]") == "aminoacetate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(=O)[O-]", "(2S)-2-methylbutanoate"),
    ],
)
def test_carboxylate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_carboxylate_is_named():
    assert smiles_to_iupac("c1ccccc1C(=O)[O-]") == "benzoate"


def test_phenyl_chain_carboxylate_ring_halogen():
    assert smiles_to_iupac("Clc1ccc(CCC(=O)[O-])cc1") == "3-(4-chlorophenyl)propanoate"


def test_phenyl_chain_carboxylate_unsaturation_is_named():
    assert smiles_to_iupac("C=Cc1ccccc1CC(=O)[O-]") == "(2-ethenylphenyl)acetate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCCCC1C(=O)O", "cycloheptanecarboxylic acid"),
    ],
)
def test_ring_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_carboxylic_acid_standalone_hydroxyl():
    assert smiles_to_iupac("OC1CCCCC1C(=O)O") == "2-hydroxycyclohexane-1-carboxylic acid"


def test_heteroaromatic_substituent_carboxylic_acid_directly_on_ring():
    assert smiles_to_iupac("OC(=O)c1cccnc1") == "pyridine-3-carboxylic acid"


def test_phenyl_substituent_carboxylic_acid_ring_branched_alkyl():
    assert smiles_to_iupac("CC(C)c1ccc(cc1)CC(=O)O") == "[4-(propan-2-yl)phenyl]acetic acid"


def test_phenyl_substituent_carboxylic_acid_ring_tert_butyl():
    assert smiles_to_iupac("CC(C)(C)c1ccc(cc1)CC(=O)O") == "(4-tert-butylphenyl)acetic acid"


def test_carboxylic_acid_enol_mix():
    assert smiles_to_iupac("OC=CC(=O)O") == "3-hydroxyprop-2-enoic acid"


def test_cyclic_carboxylic_acid_ring_stereocenter():
    assert (
        smiles_to_iupac("OC(=O)[C@H]1CCCC[C@@H]1Cl")
        == "(1R,2S)-2-chlorocyclohexane-1-carboxylic acid"
    )


def test_carboxylic_acid_ring_branch_stereocenter():
    assert (
        smiles_to_iupac("OC(=O)[C@H]1CC[C@H](C[C@@H](C)CC)C1")
        == "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carboxylic acid"
    )


def test_phenyl_substituent_carboxylic_acid_branch_tie_prefers_more_substituents():
    assert smiles_to_iupac("c1ccccc1CC(C)C(=O)O") == "2-methyl-3-phenylpropanoic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC(=O)OCC", "ethyl phenylacetate"),
        ("BrC=CCCCCCCC=CCC#CC#CCCC(=O)OC", "methyl 18-bromooctadeca-9,17-dien-4,6-diynoate"),
    ],
)
def test_ester_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_cyclyl_alcohol_ester():
    assert smiles_to_iupac("CC(=O)OC1CCCCC1") == "cyclohexyl acetate"
    assert smiles_to_iupac("CC(=O)OC1CCCC1") == "cyclopentyl acetate"
    assert smiles_to_iupac("CCC(=O)OC1CCCCC1") == "cyclohexyl propanoate"


def test_substituted_cyclyl_alcohol_ester():
    assert smiles_to_iupac("CC(=O)OC1CCC(C)CC1") == "4-methylcyclohexyl acetate"


def test_amine_coexisting_now_supported_via_ester_amine():
    assert smiles_to_iupac("NCC(=O)OC") == "methyl aminoacetate"


def test_benzoate_ester_substituted():
    assert smiles_to_iupac("Cc1ccccc1C(=O)OC") == "methyl 2-methylbenzoate"


def test_ring_acyl_ester():
    assert smiles_to_iupac("O=C(OC)C1CCCCC1") == "methyl cyclohexanecarboxylate"


def test_ring_acyl_ester_unsaturated_ring():
    assert smiles_to_iupac("O=C(OC)C1CCCC=C1") == "methyl cyclohex-2-ene-1-carboxylate"


def test_ring_acyl_ester_substituent_on_acyl_ring_atom():
    assert smiles_to_iupac("O=C(OC)C1(C)CCCCC1") == "methyl 1-methylcyclohexane-1-carboxylate"


def test_ring_acyl_chain_ester():
    assert smiles_to_iupac("O=C(OC)CC1CCCCC1") == "methyl cyclohexylacetate"


def test_ring_acyl_chain_ester_ring_with_substituent():
    assert smiles_to_iupac("O=C(OC)CC1CCC(C)CC1") == "methyl (4-methylcyclohexyl)acetate"


def test_ring_acyl_chain_ester_unsaturated_ring():
    assert smiles_to_iupac("O=C(OC)CC1CCCC=C1") == "methyl (cyclohex-2-en-1-yl)acetate"


def test_phenyl_ester_oxygen_side():
    assert smiles_to_iupac("CC(=O)Oc1ccccc1") == "phenyl acetate"
    assert smiles_to_iupac("O=COc1ccccc1") == "phenyl formate"
    assert smiles_to_iupac("CCC(=O)Oc1ccccc1") == "phenyl propanoate"


def test_phenyl_ester_oxygen_side_substituted_ring():
    assert smiles_to_iupac("CC(=O)Oc1ccccc1C") == "2-methylphenyl acetate"


def test_phenyl_acyl_chain_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(=O)OC") == "methyl (2-ethenylphenyl)acetate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(=O)OCC", "ethyl (2S)-2-methylbutanoate"),
    ],
)
def test_ester_acyl_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC(=O)C(c1ccccc1)c1ccccc1", "methyl diphenylacetate"),
        ("O=C(OC)c1ccc(Cc2ccccc2)cc1", "methyl 4-benzylbenzoate"),
        ("CC(=O)OCc1ccc(Cl)cc1", "(4-chlorophenyl)methyl acetate"),
    ],
)
def test_ester_named_from_its_alkyl_and_acid_parts(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OCC(=O)C", "2-oxopropyl acetate"),
        ("O=C(OCCN(C)C)c1ccccc1", "2-(dimethylamino)ethyl benzoate"),
        ("CC(=O)OCC(=O)O", "(acetyloxy)acetic acid"),
    ],
)
def test_ester_alkyl_part_with_heteroatom_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCOC(=O)CC(O)(CC(=O)OCC)C(=O)OCC", "triethyl 2-hydroxypropane-1,2,3-tricarboxylate"),
    ],
)
def test_esters_of_one_polyacid_cite_every_alkyl_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCOCC", "ethoxyethane"),
        ("COCC", "methoxyethane"),
    ],
)
def test_ether(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)OC(C)C", "2-(propan-2-yloxy)propane"),
        ("CC(C)COC(C)CC", "2-(2-methylpropyloxy)butane"),
    ],
)
def test_ether_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)OCC", "(2S)-2-ethoxybutane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_substituent_branch():
    assert smiles_to_iupac("CCCCCO[C@H](C)CC") == "1-{[(2R)-butan-2-yl]oxy}pentane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1OC", "methoxybenzene"),
        ("c1ccccc1OC(C)C", "(propan-2-yloxy)benzene"),
        ("c1ccccc1COCC", "(ethoxymethyl)benzene"),
        ("c1ccccc1COC(C)C", "[(propan-2-yloxy)methyl]benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_stereocenter():
    assert smiles_to_iupac("c1ccccc1O[C@H](C)CC") == "{[(2R)-butan-2-yl]oxy}benzene"


def test_thiophosphoric_acid():
    assert smiles_to_iupac("OP(=S)(O)O") == "phosphorothioic O,O,O-acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CS(=O)C=C1", "thiophene 1-oxide"),  # PubChem CID 9548690
    ],
)
def test_hetero_ring_oxide_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_dihydrate_adduct():
    # Blue Book P-14.8: the proportion (1/n) is always cited.
    assert smiles_to_iupac("OC(=O)C(=O)O.O.O") == "oxalic acid—water (1/2)"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COO", "methaneperoxol"),
        ("ClCCOO", "2-chloroethane-1-peroxol"),
    ],
)
def test_hydroperoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCOO")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOC1CCCCC1")


def test_phenyl_chain_hydroperoxide_internal_locant():
    assert smiles_to_iupac("C(c1ccccc1)C(C)OO") == "1-phenylpropane-2-peroxol"


def test_hydroperoxybenzene():
    assert smiles_to_iupac("c1ccccc1OO") == "hydroperoxybenzene"


def test_phenyl_substituted_benzene_ring_hydroperoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CCOO")


def test_phenyl_chain_hydroperoxide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCOO")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OOc1ccco1", "2-hydroperoxyfuran"),
    ],
)
def test_heteroaromatic_ring_direct_attachment_hydroperoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_heteroaromatic_ring_chain_hydroperoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOCCc1cccnc1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCO", "ethanol"),
        ("OCCCl", "2-chloroethan-1-ol"),
        ("NCCS", "2-aminoethane-1-thiol"),
    ],
)
def test_two_carbon_suffix_locant_is_cited_once_a_prefix_is_present(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aryl_ketone_cites_the_ring_as_a_substituent():
    assert smiles_to_iupac("CC(=O)c1ccccc1") == "1-phenylethan-1-one"
    assert smiles_to_iupac("O=C(c1ccccc1)c1ccccc1") == "diphenylmethanone"
    assert smiles_to_iupac("O=C(c1ccc(Cl)cc1)c1ccccc1") == "(4-chlorophenyl)phenylmethanone"


def test_heteroaromatic_substituent_ketone_directly_on_ring():
    assert smiles_to_iupac("CC(=O)c1cccnc1") == "1-(pyridin-3-yl)ethan-1-one"


def test_ether_now_supported_via_ether_ketone():
    assert smiles_to_iupac("CCOCC(=O)C") == "1-ethoxypropan-2-one"


def test_hetero_ring_ketone_non_alkyl_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNCCCN(C)S(=O)(=O)N1CCC(OC)C1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1COCCN1", "morpholin-3-one"),
        ("O=C1COC(=O)CO1", "1,4-dioxane-2,5-dione"),
    ],
)
def test_two_hetero_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_hetero_ring_ketone_substituted_heteroatom():
    assert smiles_to_iupac("O=C1COCCN1C") == "4-methylmorpholin-3-one"


def test_two_hetero_ring_ketone_wrong_element_pair():
    assert smiles_to_iupac("O=C1COCC[Se]1") == "1,4-oxaselenan-3-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CNCCNC1=O", "1,4-diazepan-5-one"),
        ("C1COCCNC1=O", "1,4-oxazepan-5-one"),
    ],
)
def test_seven_membered_1_4_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCNC(=O)NC1", "1,3-diazepan-2-one"),
    ],
)
def test_seven_membered_1_3_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_membered_1_3_ring_ketone_second_ketone_on_far_arc():
    assert smiles_to_iupac("O=C1NC(=O)CCCN1") == "1,3-diazepane-2,4-dione"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("N1NC(=O)CCCC1", "1,2-diazepan-3-one"),
    ],
)
def test_seven_membered_1_2_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_membered_1_2_ring_ketone_second_ketone_on_far_arc():
    assert smiles_to_iupac("N1NC(=O)CCC(=O)C1") == "1,2-diazepane-3,6-dione"


def test_seven_membered_1_2_ring_ketone_wrong_element_pair():
    assert smiles_to_iupac("O=C1CCCC[Se]N1") == "1,2-selenazepan-3-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1OCCN1", "1,3-oxazolidin-2-one"),
        ("O=C1SCCO1", "1,3-oxathiolan-2-one"),
    ],
)
def test_five_membered_1_3_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_3_ring_ketone_branched_n_substituent():
    assert smiles_to_iupac("O=C1N(C(C)C)C(=O)CN1") == "3-(propan-2-yl)imidazolidine-2,4-dione"


def test_five_membered_1_3_ring_ketone_substituted_n_substituent():
    assert smiles_to_iupac("O=C1N(CCO)C(=O)CN1") == "3-(2-hydroxyethyl)imidazolidine-2,4-dione"


def test_five_membered_1_3_ring_ketone_single_ketone_ring_carbon_substituent_tie_break():
    assert smiles_to_iupac("O=C1NC(C)CN1") == "4-methylimidazolidin-2-one"


def test_five_membered_1_3_ring_ketone_wrong_element_pair():
    assert smiles_to_iupac("O=C1OCC[Se]1") == "1,3-oxaselenolan-2-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1CCNN1", "pyrazolidin-3-one"),
    ],
)
def test_five_membered_1_2_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_2_ring_ketone_wrong_element_pair():
    assert smiles_to_iupac("O=C1CC[Se]N1") == "1,2-selenazolidin-3-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1C(C)C(=O)ON1", "4-methyl-1,2-oxazolidine-3,5-dione"),
    ],
)
def test_five_membered_1_2_ring_dione_carbon_substituent_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_alcohol_hetero_mix_names_hydroxy_prefix():
    assert smiles_to_iupac("OCC(=O)C") == "1-hydroxypropan-2-one"


def test_cyclic_ketone_alcohol_mix_names_hydroxy_prefix():
    assert smiles_to_iupac("OC1CCC(=O)CC1") == "4-hydroxycyclohexan-1-one"


def test_ketone_enol_mix():
    assert smiles_to_iupac("OC=CC(=O)C") == "4-hydroxybut-3-en-2-one"


def test_unsaturated_ring_ketone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC#C1")


def test_ring_substituent_chain_ketone_unsaturated_ring():
    assert smiles_to_iupac("CC(=O)C1CCCC=C1") == "1-(cyclohex-2-en-1-yl)ethan-1-one"


def test_ring_with_ketone_chain_ketone_tie():
    assert smiles_to_iupac("O=C1CCCCC1CC(C)=O") == "2-(2-oxopropyl)cyclohexan-1-one"


def test_ring_with_ketone_chain_ketone_chain_wins():
    assert smiles_to_iupac("O=C1CCCCC1C(=O)CC(C)=O") == "1-(2-oxocyclohexyl)butane-1,3-dione"


def test_cyclic_ketone_stereocenter():
    assert smiles_to_iupac("O=C1CCCC[C@H]1Cl") == "(2R)-2-chlorocyclohexan-1-one"


def test_cyclic_ketone_branch_stereocenter():
    assert smiles_to_iupac("O=C1CCC(CC1)[C@@H](C)CC") == "4-[(2S)-butan-2-yl]cyclohexan-1-one"


def test_acyclic_ketone_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("CC(=O)[C@H](C)/C=C/C") == "(3R,4E)-3-methylhex-4-en-2-one"
    assert smiles_to_iupac("CC(=O)[C@@H](C)/C=C/C") == "(3S,4E)-3-methylhex-4-en-2-one"


def test_two_ring_aromatic_substituent_ketone_ring_bond():
    assert smiles_to_iupac("O=C1C=CCCC1c1ccccc1") == "6-phenylcyclohex-2-en-1-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COOCC", "methylperoxyethane"),
    ],
)
def test_peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported__peroxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOC1CCCCC1")


def test_unsaturated_not_supported__peroxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=COOC")


def test_three_oxygens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COOCOC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)COOC(C)(C)C", "1-tert-butylperoxy-2-methylpropane"),
        ("CC(C)COOC(C)CC", "2-(2-methylpropyl)peroxybutane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)OOCC", "(2S)-2-ethylperoxybutane"),
    ],
)
def test_stereocenter_on_parent_chain__peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCOO[C@H](C)CC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1OOC", "methylperoxybenzene"),
        ("c1ccccc1OOC(C)C", "(propan-2-yl)peroxybenzene"),
        ("c1ccccc1COOCC", "(ethylperoxymethyl)benzene"),
        ("c1ccccc1COOC(C)C", "[(propan-2-yl)peroxymethyl]benzene"),
    ],
)
def test_benzene_ring_parent__peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_stereocenter_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1OO[C@H](C)CC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NC(=O)C(=O)O", "oxamic acid"),
        ("OC(=O)C(O)(c1ccccc1)c1ccccc1", "hydroxydi(phenyl)acetic acid"),
        ("O=C(O)[N+](=O)[O-]", "nitroformic acid"),
        ("OC(=O)N1CCCC1", "pyrrolidine-1-carboxylic acid"),
        ("NNC(=O)O", "hydrazinecarboxylic acid"),
        ("CNNC(=O)O", "2-methylhydrazine-1-carboxylic acid"),
        ("OC(=O)[SiH2]O[SiH2]C(=O)O", "disiloxane-1,3-dicarboxylic acid"),
        ("OC(=O)[SiH2][SiH2][SiH3]", "trisilane-1-carboxylic acid"),
        ("OC(=O)[SiH2]C(=O)O", "silanedicarboxylic acid"),
        ("OC(=O)N(C)NC", "1,2-dimethylhydrazine-1-carboxylic acid"),
        ("OC(=O)N(C(=O)O)N(C)C(=O)O", "2-methylhydrazine-1,1,2-tricarboxylic acid"),
        ("OC(=O)C(=O)C1CCCCC1C(=O)O", "2-oxalocyclohexane-1-carboxylic acid"),
    ],
)
def test_acid_retained_names_and_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCCCCC(=O)SCC", "S-ethyl hexanethioate"),
        ("CC(=N)OC", "methyl ethanimidate"),
        ("CC(=Nc1ccccc1)OC", "methyl N-phenylethanimidate"),
        ("N=C(SC)c1ccccc1", "methyl benzenecarboximidothioate"),
        ("CC(=O)O[Si](C)(C)C", "trimethylsilyl acetate"),
        ("CCOP(C)(=O)C", "ethyl dimethylphosphinate"),
        ("CCOP(=O)(C)O", "ethyl hydrogen methylphosphonate"),
        ("CCOB(OCC)OCC", "triethyl borate"),
        ("CCOP(c1ccccc1)OCC", "diethyl phenylphosphonite"),
        ("CCOP(C)C", "ethyl dimethylphosphinite"),
        ("CC(=O)OP(C)(C)=O", "acetic dimethylphosphinic anhydride"),
        ("CC(=O)OOC(C)(C)C", "tert-butyl ethaneperoxoate"),
        ("COC(=S)OC", "O,O-dimethyl carbonothioate"),
        ("CSC(=S)OC", "O,S-dimethyl carbonodithioate"),
        ("CN=C(O)OC", "methyl hydrogen N-methylcarbonimidate"),
        ("CCOC(=O)NCC", "ethyl ethylcarbamate"),
        ("COC(=O)OC(=O)OC", "dimethyl dicarbonate"),
        ("CS(=O)(=O)SCC", "S-ethyl methanesulfonothioate"),
    ],
)
def test_ester_class_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=S)OC(=S)C", "ethanethioic anhydride"),
        ("CS(=O)(=O)OS(=O)(=O)C", "methanesulfonic anhydride"),
        ("CC(=O)OC(=O)CCC(=O)OC(C)=O", "diacetic butanedioic dianhydride"),
        ("CC(=O)OC(=O)CCC(=O)OC(=O)CCC(=O)OC(=O)CC", "acetic butanedioic 4-oxo-4-(propanoyloxy)butanoic dianhydride"),
        ("CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)c(C(=O)OC(=O)CC)cc1", "2,4-diacetic 1-propanoic benzene-1,2,4-tricarboxylic trianhydride"),
        ("O=C(Cl)c1ccc(C(=O)Cl)cc1", "benzene-1,4-dicarbonyl dichloride"),
        ("NC(=O)C(=O)Br", "oxamoyl bromide"),
        ("ClC(=S)C1CCCCC1", "cyclohexanecarbothioyl chloride"),
        ("CCCC(=O)C#N", "butanoyl cyanide"),
        ("O=C(N=C=S)C(=O)N=C=S", "oxalyl diisothiocyanate"),
        ("CC(=O)OP(=O)(O)O", "(acetyloxy)phosphonic acid"),
        ("CCC(=O)OB(O)O", "(propanoyloxy)boronic acid"),
        ("CC(=O)OP(=O)(OC)OC", "acetic (dimethyl hydrogen phosphate) anhydride"),
        ("CC(=O)O[As](C)(C)=O", "acetic dimethylarsinic anhydride"),
        ("CO[Sb](OC)OC", "trimethyl stiborite"),
        ("ClC(=O)Cl", "carbonyl dichloride"),
        ("N#CC(=O)Cl", "carbonocyanidoyl chloride"),
        ("N#CCl", "carbononitridic chloride"),
        ("ClC(=O)OC(=O)Cl", "dicarbonic dichloride"),
    ],
)
def test_anhydride_and_acyl_halide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccccc1C(=O)Cl", "2-carbonochloridoylbenzoic acid"),
        ("OC(=O)c1ccccc1C(=O)OO", "2-carbonoperoxoylbenzoic acid"),
        ("OC(=O)c1ccccc1C(=O)S", "2-(sulfanylcarbonyl)benzoic acid"),
        ("OC(=O)c1ccccc1C(O)=S", "2-(hydroxycarbonothioyl)benzoic acid"),
        ("OC(=O)c1ccccc1C(=N)O", "2-(C-hydroxycarbonimidoyl)benzoic acid"),
        ("OC(=O)c1ccccc1C(=NO)O", "2-(C,N-dihydroxycarbonimidoyl)benzoic acid"),
        ("OC(=O)c1ccccc1C(=O)OS", "2-[(sulfanyloxy)carbonyl]benzoic acid"),
        ("OC(=O)c1ccccc1C(=O)SO", "2-[(hydroxysulfanyl)carbonyl]benzoic acid"),
        ("OC(=O)c1ccccc1C(C)=S", "2-(ethanethioyl)benzoic acid"),
        ("OC(=O)CCC(=N)O", "4-hydroxy-4-iminobutanoic acid"),
        ("OC(=O)CCCC(=NO)O", "5-hydroxy-5-(hydroxyimino)pentanoic acid"),
    ],
)
def test_acid_derivative_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)OO", "carbonoperoxoic acid"),
        ("NC(=N)S", "carbamimidothioic acid"),
        ("[N-]=[N+]=NC(=O)O", "carbonazidic acid"),
        ("OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)O", "3,5,7-trioxo-2,4,6,8-tetraoxanonanedioic acid"),
        ("OC(=S)OC(=S)O", "1,3-dithiodicarbonic O1,O3-acid"),
        ("OC(=O)SC(=O)O", "2-thiodicarbonic acid"),
        ("OC(=N)OC(=O)O", "1-imidodicarbonic acid"),
        ("ClC(=O)OC(=O)O", "chlorodicarbonic acid"),
    ],
)
def test_carbonic_family_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC(=O)CCC(=O)OCCOC(=O)CCC(=O)OC", "dimethyl ethane-1,2-diyl dibutanedioate"),
        ("CC(=O)Oc1ccc(cc1)C(=O)OC", "methyl 4-(acetyloxy)benzoate"),
        ("O=C(OCCOC(C)=O)CC(=O)OCCOC(C)=O", "bis[2-(acetyloxy)ethyl] propanedioate"),
    ],
)
def test_polyester_principal_acid_and_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
