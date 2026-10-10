import re

import pytest
from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(OCC)OCC", "1,1-diethoxypropane"),
    ],
)
def test_acetal_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1CCC(OC)(OC)CC1", "1,1-dimethoxycyclohexane", id="cyclic_acetal"),
        pytest.param("CCC(OC(C)C)OCC", "1-ethoxy-1-[(propan-2-yl)oxy]propane", id="branched_alkoxy"),
        pytest.param("C=CC(OCC)OCC", "3,3-diethoxyprop-1-ene", id="unsaturated_acetal"),
        pytest.param("OC(OCC)CC", "1-ethoxypropan-1-ol", id="hemiacetal_named_as_alkoxy_alcohol"),
    ],
)
def test_cyclic_acetal_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(OC)OC", "(2S)-1,1-dimethoxy-2-methylbutane"),
    ],
)
def test_acetal_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=CC(=O)Cl", "prop-2-enoyl chloride", id="unsaturated_acyl_halide"),
        pytest.param("CC(C)CC(=O)Cl", "3-methylbutanoyl chloride", id="branched_substituent"),
        pytest.param("CC(C[C@H](C)Cl)C(=O)Cl", "(4S)-4-chloro-2-methylpentanoyl chloride", id="branch_stereocenter_cites_the_specified_elements"),
    ],
)
def test_unsaturated_acyl_halide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC(=O)Cl", "phenylacetyl chloride"),
    ],
)
def test_phenyl_chain_acyl_halide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[C@H](O)/C=C/C", "(2S,3E)-pent-3-en-2-ol", id="stereocenter_with_ez_double_bond_coexistence"),
        pytest.param("O[C@H]1CCCCC1C", "(1S)-2-methylcyclohexan-1-ol", id="ring_stereocenter_cites_the_specified_elements"),
        pytest.param("OC1(CCCCC1)[C@@H](C)CC", "1-[(2S)-butan-2-yl]cyclohexan-1-ol", id="ring_branch_stereocenter_supported"),
        pytest.param("OC1(CCCCC1)[C@@H]([C@H](C)Cl)CC", "1-[(2S,3S)-2-chloropentan-3-yl]cyclohexan-1-ol", id="ring_branch_two_stereocenters"),
        pytest.param("O[C@H]1CCCC2CCCC12", "(4S)-octahydro-1H-inden-4-ol", id="polycyclic_ring_stereocenter_cites_the_specified_elements"),
        pytest.param("OC1C(O)CC(CO)CC1", "4-(hydroxymethyl)cyclohexane-1,2-diol", id="ring_hydroxyl_count_exceeds_chain"),
        pytest.param("OC1CCCCC1C(CO)CO", "2-(2-hydroxycyclohexyl)propane-1,3-diol", id="ring_with_hydroxyl_and_branched_diol_chain"),
    ],
)
def test_stereocenter_with_ez_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1C(C)(CO)CO", "2-(2-hydroxycyclohexyl)-2-methylpropane-1,3-diol"),
        ("OC1CCCCC1CC(CO)CO", "2-[(2-hydroxycyclohexyl)methyl]propane-1,3-diol"),
        ("OC1CCCCC1C(CO)(CO)CO", "2-(2-hydroxycyclohexyl)-2-(hydroxymethyl)propane-1,3-diol"),
    ],
)
def test_chain_with_most_hydroxyls_is_parent_over_longer_chain_and_ring_hydroxyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=Cc1ccccc1CO", "(2-ethenylphenyl)methanol", id="phenyl_chain_unsaturation"),
        pytest.param("Oc1cccnc1", "pyridin-3-ol", id="heteroaromatic_substituent_alcohol_directly_on_ring"),
        pytest.param("Oc1ccccc1O", "benzene-1,2-diol", id="phenol_multiple_hydroxyls"),
    ],
)
def test_phenyl_chain_unsaturation_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1CCCCC1c1ccc[nH]1", "2-(1H-pyrrol-2-yl)cyclohexan-1-ol"),
    ],
)
def test_two_ring_aromatic_substituent_alcohol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OCC1CCCCC1c1ccccc1", "(2-phenylcyclohexyl)methanol", id="two_ring_aromatic_substituent_alcohol_chain_hydroxyl"),
        pytest.param("OC1CCCCC1OC", "2-methoxycyclohexan-1-ol", id="ether_on_ring"),
        pytest.param("OC=CC", "prop-1-en-1-ol", id="enol"),
        pytest.param("OC1CCCC=C1C", "2-methylcyclohex-2-en-1-ol", id="unsaturated_ring_alcohol_with_substituent"),
    ],
)
def test_two_ring_aromatic_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C1C#CCC1", "cyclopent-2-yn-1-one", id="ketone_beside_a_ring_triple_bond"),
        pytest.param("O=CC1C#CC1", "cyclobut-2-yne-1-carbaldehyde", id="carbaldehyde_suffix_beside_a_ring_triple_bond"),
        pytest.param("ClC1CCCC#CCC1O", "8-chlorocyclooct-3-yn-1-ol", id="suffix_before_triple_bond_before_prefix"),
        pytest.param("O=C1CCCC#CCCC#CCC1", "cyclododeca-4,8-diyn-1-one", id="several_ring_triple_bonds"),
        pytest.param("OC1CCC=CCC#CC1", "cyclonon-6-en-3-yn-1-ol", id="ene_and_yne_lowest_locants_together"),
        pytest.param("O=C1C=CCC#C1", "cyclohex-2-en-5-yn-1-one", id="double_bond_lower_when_the_set_ties"),
    ],
)
def test_ring_triple_bond_beside_a_suffix_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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




def test_benzaldehyde():
    assert smiles_to_iupac("O=Cc1ccccc1") == "benzaldehyde"
    assert smiles_to_iupac("O=Cc1ccc(C)cc1") == "4-methylbenzaldehyde"
    assert smiles_to_iupac("O=Cc1ccccc1Cl") == "2-chlorobenzaldehyde"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OCc1ccccc1CC=O", "[2-(hydroxymethyl)phenyl]acetaldehyde", id="with_hydroxyl"),
        pytest.param("C=Cc1ccccc1CC=O", "(2-ethenylphenyl)acetaldehyde", id="unsaturation"),
    ],
)
def test_phenyl_chain_aldehyde_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CC1CCCCC1c1ccccc1", "2-phenylcyclohexane-1-carbaldehyde"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_aldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=CC1CCCC(c2ccc(C)cc2)C1", "3-(4-methylphenyl)cyclohexane-1-carbaldehyde", id="two_ring_aromatic_substituent_aldehyde_substituted_ring"),
        pytest.param("O=CCC1CCCCC1c1ccccc1", "(2-phenylcyclohexyl)acetaldehyde", id="two_ring_aromatic_substituent_aldehyde_chain_aldehyde"),
        pytest.param("O=Cc1cccnc1", "pyridine-3-carbaldehyde", id="heteroaromatic_ring_directly_attached_aldehyde"),
        pytest.param("OC=CC=O", "3-hydroxyprop-2-enal", id="aldehyde_enol_mix"),
        pytest.param("O=CC1CCC(C=O)CC1", "cyclohexane-1,4-dicarbaldehyde", id="ring_aldehyde_multiple_groups"),
    ],
)
def test_two_ring_aromatic_and_related_2(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_aldehyde_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("O=C[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enal"
    assert smiles_to_iupac("O=C[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enal"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C[C@H]1CCCC[C@@H]1Cl", "(1R,2S)-2-chlorocyclohexane-1-carbaldehyde", id="cyclic_aldehyde_ring_stereocenter"),
        pytest.param("O=C[C@H]1CC[C@H](C[C@@H](C)CC)C1", "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carbaldehyde", id="aldehyde_ring_branch_stereocenter"),
    ],
)
def test_cyclic_aldehyde_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=COC=O", "formic anhydride"),
    ],
)
def test_symmetric_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)C(=O)OC(=O)C(C)C", "2-methylpropanoic anhydride", id="branched_symmetric_anhydride"),
        pytest.param("ClCC(=O)OC(=O)CCl", "chloroacetic anhydride", id="halogen_substituent"),
        pytest.param("O=C1CCC(=O)O1", "oxolane-2,5-dione", id="cyclic_anhydride_names_as_ring_dione"),
        pytest.param("OC(=O)O", "carbonic acid", id="carbonic_acid"),
        pytest.param("COC(=O)OC", "dimethyl carbonate", id="dimethyl_carbonate"),
        pytest.param("COC(=O)O", "methyl hydrogen carbonate", id="methyl_hydrogen_carbonate"),
    ],
)
def test_branched_symmetric_anhydride_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # -oate combined with existing unsaturation support.
        ("C=CC(=O)[O-]", "prop-2-enoate"),
    ],
)
def test_carboxylate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[O-]C(=O)CC(=O)[O-]", "propanedioate", id="dicarboxylate_is_named"),
        pytest.param("[O-]C(=O)C1CCCCC1", "cyclohexanecarboxylate", id="ring_carboxylate_is_named"),
        pytest.param("NCC(=O)[O-]", "glycinate", id="amine_coexisting_is_named"),
    ],
)
def test_dicarboxylate_is_named_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C(=O)[O-]", "(2S)-2-methylbutanoate"),
    ],
)
def test_carboxylate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1C(=O)[O-]", "benzoate", id="directly_attached_carboxylate_is_named"),
        pytest.param("Clc1ccc(CCC(=O)[O-])cc1", "3-(4-chlorophenyl)propanoate", id="chain_carboxylate_ring_halogen"),
        pytest.param("C=Cc1ccccc1CC(=O)[O-]", "(2-ethenylphenyl)acetate", id="chain_carboxylate_unsaturation_is_named"),
    ],
)
def test_phenyl_directly_attached_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCCCC1C(=O)O", "cycloheptanecarboxylic acid"),
    ],
)
def test_ring_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OC1CCCCC1C(=O)O", "2-hydroxycyclohexane-1-carboxylic acid", id="ring_carboxylic_acid_standalone_hydroxyl"),
        pytest.param("OC(=O)c1cccnc1", "pyridine-3-carboxylic acid", id="heteroaromatic_substituent_carboxylic_acid_directly_on_ring"),
        pytest.param("CC(C)c1ccc(cc1)CC(=O)O", "[4-(propan-2-yl)phenyl]acetic acid", id="phenyl_substituent_carboxylic_acid_ring_branched_alkyl"),
        pytest.param("CC(C)(C)c1ccc(cc1)CC(=O)O", "(4-tert-butylphenyl)acetic acid", id="phenyl_substituent_carboxylic_acid_ring_tert_butyl"),
        pytest.param("OC=CC(=O)O", "3-hydroxyprop-2-enoic acid", id="carboxylic_acid_enol_mix"),
        pytest.param("OC(=O)[C@H]1CCCC[C@@H]1Cl", "(1R,2S)-2-chlorocyclohexane-1-carboxylic acid", id="cyclic_carboxylic_acid_ring_stereocenter"),
        pytest.param("OC(=O)[C@H]1CC[C@H](C[C@@H](C)CC)C1", "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carboxylic acid", id="carboxylic_acid_ring_branch_stereocenter"),
        pytest.param("c1ccccc1CC(C)C(=O)O", "2-methyl-3-phenylpropanoic acid", id="phenyl_substituent_carboxylic_acid_branch_tie_prefers_more_substituents"),
    ],
)
def test_ring_carboxylic_acid_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)OC1CCC(C)CC1", "4-methylcyclohexyl acetate", id="substituted_cyclyl_alcohol_ester"),
        pytest.param("NCC(=O)OC", "methyl glycinate", id="amine_coexisting_now_supported_via_ester_amine"),
        pytest.param("Cc1ccccc1C(=O)OC", "methyl 2-methylbenzoate", id="benzoate_ester_substituted"),
        pytest.param("O=C(OC)C1CCCCC1", "methyl cyclohexanecarboxylate", id="ring_acyl_ester"),
        pytest.param("O=C(OC)C1CCCC=C1", "methyl cyclohex-2-ene-1-carboxylate", id="ring_acyl_ester_unsaturated_ring"),
        pytest.param("O=C(OC)C1(C)CCCCC1", "methyl 1-methylcyclohexane-1-carboxylate", id="ring_acyl_ester_substituent_on_acyl_ring_atom"),
        pytest.param("O=C(OC)CC1CCCCC1", "methyl cyclohexylacetate", id="ring_acyl_chain_ester"),
        pytest.param("O=C(OC)CC1CCC(C)CC1", "methyl (4-methylcyclohexyl)acetate", id="ring_acyl_chain_ester_ring_with_substituent"),
        pytest.param("O=C(OC)CC1CCCC=C1", "methyl (cyclohex-2-en-1-yl)acetate", id="ring_acyl_chain_ester_unsaturated_ring"),
    ],
)
def test_substituted_cyclyl_alcohol_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_ester_oxygen_side():
    assert smiles_to_iupac("CC(=O)Oc1ccccc1") == "phenyl acetate"
    assert smiles_to_iupac("O=COc1ccccc1") == "phenyl formate"
    assert smiles_to_iupac("CCC(=O)Oc1ccccc1") == "phenyl propanoate"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)Oc1ccccc1C", "2-methylphenyl acetate", id="ester_oxygen_side_substituted_ring"),
        pytest.param("C=Cc1ccccc1CC(=O)OC", "methyl (2-ethenylphenyl)acetate", id="acyl_chain_unsaturation"),
    ],
)
def test_phenyl_ester_oxygen_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
        pytest.param("CC(NCC(=O)OCc1ccccc1)C(=O)OCc1ccccc1", "benzyl N-[2-(benzyloxy)-2-oxoethyl]-alaninate", id="ester_on_substituent_acid"),
        pytest.param(
            "CC(C)CC(NCC(=O)NC(C)C(=O)OCc1ccccc1)C(=O)OCc1ccccc1",
            "benzyl N-(2-{[1-(benzyloxy)-1-oxopropan-2-yl]amino}-2-oxoethyl)-leucinate",
            id="ester_on_substituent_beyond_amide",
        ),
        pytest.param("O=C(OCc1ccccc1)CCC(=O)OCc1ccccc1", "dibenzyl butanedioate", id="both_esters_of_one_acid"),
    ],
)
def test_ester_outside_the_parent_acid_is_a_prefix(smiles, expected):
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
        ("CC(C)OC(C)C", "2-[(propan-2-yl)oxy]propane"),
        ("CC(C)COC(C)CC", "2-(2-methylpropoxy)butane"),
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
        ("c1ccccc1OC", "anisole"),
        ("c1ccc(Cl)cc1OC", "1-chloro-3-methoxybenzene"),
        ("c1ccccc1OC(C)C", "[(propan-2-yl)oxy]benzene"),
        ("c1ccccc1COCC", "(ethoxymethyl)benzene"),
        ("c1ccccc1COC(C)C", "{[(propan-2-yl)oxy]methyl}benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1O[C@H](C)CC", "{[(2R)-butan-2-yl]oxy}benzene", id="benzene_ring_stereocenter"),
        pytest.param("OP(=S)(O)O", "phosphorothioic O,O,O-acid", id="thiophosphoric_acid"),
    ],
)
def test_benzene_ring_stereocenter_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CS(=O)C=C1", "1H-1λ4-thiophen-1-one"),  # PubChem CID 9548690 (thiophene 1-oxide)
    ],
)
def test_hetero_ring_oxide_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_dihydrate_adduct():
    # Blue Book P-14.8: the proportion (1/n) is always cited.
    assert smiles_to_iupac("OC(=O)C(=O)O.O.O") == "oxalic acid—water (1/2)"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCO.c1ccncc1", "ethanol—pyridine (1/1)"),
        ("CC(=O)O.CCN", "acetic acid—ethanamine (1/1)"),
        ("CC.CCOCC", "ethane—ethoxyethane (1/1)"),
        ("F.O=c1cccc[nH]1", "pyridin-2(1H)-one—hydrogen fluoride (1/1)"),
        pytest.param(
            "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61.[O-][N+](=O)c1cc([N+]([O-])=O)cc([N+]([O-])=O)c1",
            "coronene—1,3,5-trinitrobenzene (1/1)",
            id="nitro_groups_leave_a_carbon_compound_ranked_by_its_ring_system",
        ),
        pytest.param(
            "OC(=O)c1cc2ccccc2[nH]1.OC(=O)c1cc(cc(c1)[N+]([O-])=O)[N+]([O-])=O",
            "1H-indole-2-carboxylic acid—3,5-dinitrobenzoic acid (1/1)",
            id="heterocycle_before_carbocycle_within_one_class",
        ),
    ],
)
def test_adduct_components_ordered_by_class_seniority(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CO.[Pd]", "methanol—palladium (1/1)", id="metal_atom"),
        pytest.param("C.[Xe]", "methane—xenon (1/1)", id="noble_gas_atom"),
        pytest.param("C.[HH]", "methane—dihydrogen (1/1)", id="dihydrogen"),
        pytest.param("C.[H+]", "methane—hydron (1/1)", id="hydron"),
        pytest.param("C.[H-]", "methane—hydride (1/1)", id="hydride"),
        pytest.param("C.O.[Xe]", "methane—xenon—water (1/1/1)", id="water_cited_last"),
    ],
)
def test_adduct_with_inorganic_components(smiles, expected):
    with pytest.warns(NonPreferredNameWarning, match="P-14.8.2"):
        assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COO", "methaneperoxol"),
        ("ClCCOO", "2-chloroethane-1-peroxol"),
    ],
)
def test_hydroperoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=CCOO", "prop-2-ene-1-peroxol", id="unsaturated_chain"),
    ],
)
def test_unsaturated_peroxol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCC(=O)[Se][Se][Se]C", "1-(methyltriselanyl)propan-1-one"),
        ("CCC(=O)OOO", "1-trioxidanylpropan-1-one"),
        ("CSSSC(=O)C", "1-(methyltrisulfanyl)ethan-1-one"),
    ],
)
def test_acyl_on_a_homogeneous_chalcogen_chain_is_a_pseudoketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C(c1ccccc1)C(C)OO", "1-phenylpropane-2-peroxol", id="phenyl_chain_hydroperoxide_internal_locant"),
        pytest.param("c1ccccc1OO", "benzeneperoxol", id="ring_peroxol"),
        pytest.param("OOC1CCC(N)CC1", "4-aminocyclohexane-1-peroxol", id="saturated_ring_peroxol"),
    ],
)
def test_phenyl_chain_hydroperoxide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Cc1ccccc1CCOO", "2-(2-methylphenyl)ethane-1-peroxol", id="substituted_benzene_ring"),
        pytest.param("C=Cc1ccccc1CCOO", "2-(2-ethenylphenyl)ethane-1-peroxol", id="unsaturated_ring_substituent"),
    ],
)
def test_peroxol_on_a_chain_beside_a_substituted_benzene_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OOc1ccco1", "furan-2-peroxol"),
    ],
)
def test_heteroaromatic_ring_direct_attachment_hydroperoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_heteroaromatic_ring_chain_hydroperoxide():
    assert smiles_to_iupac("OOCCc1cccnc1") == "2-(pyridin-3-yl)ethane-1-peroxol"


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
    assert smiles_to_iupac("O=C(c1ccc(Cl)cc1)c1ccccc1") == "(4-chlorophenyl)(phenyl)methanone"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)c1cccnc1", "1-(pyridin-3-yl)ethan-1-one", id="heteroaromatic_substituent_ketone_directly_on_ring"),
        pytest.param("CCOCC(=O)C", "1-ethoxypropan-2-one", id="ether_now_supported_via_ether_ketone"),
    ],
)
def test_heteroaromatic_substituent_ketone_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_sulfonamide_with_two_n_substituents():
    assert (
        smiles_to_iupac("CCNCCCN(C)S(=O)(=O)N1CCC(OC)C1")
        == "N-[3-(ethylamino)propyl]-3-methoxy-N-methylpyrrolidine-1-sulfonamide"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1COCCN1", "morpholin-3-one"),
        ("O=C1COC(=O)CO1", "1,4-dioxane-2,5-dione"),
    ],
)
def test_two_hetero_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C1COCCN1C", "4-methylmorpholin-3-one", id="substituted_heteroatom"),
        pytest.param("O=C1COCC[Se]1", "1,4-oxaselenan-3-one", id="wrong_element_pair"),
    ],
)
def test_two_hetero_ring_ketone_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("N1NC(=O)CCC(=O)C1", "1,2-diazepane-3,6-dione", id="second_ketone_on_far_arc"),
        pytest.param("O=C1CCCC[Se]N1", "1,2-selenazepan-3-one", id="wrong_element_pair"),
    ],
)
def test_seven_membered_1_2_ring_ketone_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1OCCN1", "1,3-oxazolidin-2-one"),
        ("O=C1SCCO1", "1,3-oxathiolan-2-one"),
    ],
)
def test_five_membered_1_3_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C1N(C(C)C)C(=O)CN1", "3-(propan-2-yl)imidazolidine-2,4-dione", id="branched_n_substituent"),
        pytest.param("O=C1N(CCO)C(=O)CN1", "3-(2-hydroxyethyl)imidazolidine-2,4-dione", id="substituted_n_substituent"),
        pytest.param("O=C1NC(C)CN1", "4-methylimidazolidin-2-one", id="single_ketone_ring_carbon_substituent_tie_break"),
        pytest.param("O=C1OCC[Se]1", "1,3-oxaselenolan-2-one", id="wrong_element_pair"),
    ],
)
def test_five_membered_1_3_ring_ketone_cases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OCC(=O)C", "1-hydroxypropan-2-one", id="alcohol_hetero_mix_names_hydroxy_prefix"),
        pytest.param("OC1CCC(=O)CC1", "4-hydroxycyclohexan-1-one", id="cyclic_ketone_alcohol_mix_names_hydroxy_prefix"),
        pytest.param("OC=CC(=O)C", "4-hydroxybut-3-en-2-one", id="ketone_enol_mix"),
    ],
)
def test_alcohol_hetero_mix_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)C1CCCC=C1", "1-(cyclohex-2-en-1-yl)ethan-1-one", id="ring_substituent_chain_ketone_unsaturated_ring"),
        pytest.param("O=C1CCCCC1CC(C)=O", "2-(2-oxopropyl)cyclohexan-1-one", id="ring_with_ketone_chain_ketone_tie"),
        pytest.param("O=C1CCCCC1C(=O)CC(C)=O", "1-(2-oxocyclohexyl)butane-1,3-dione", id="ring_with_ketone_chain_ketone_chain_wins"),
        pytest.param("O=C1CCCC[C@H]1Cl", "(2R)-2-chlorocyclohexan-1-one", id="cyclic_ketone_stereocenter"),
        pytest.param("O=C1CCC(CC1)[C@@H](C)CC", "4-[(2S)-butan-2-yl]cyclohexan-1-one", id="cyclic_ketone_branch_stereocenter"),
    ],
)
def test_ring_substituent_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_ketone_stereocenter_with_ez_double_bond_coexistence():
    assert smiles_to_iupac("CC(=O)[C@H](C)/C=C/C") == "(3R,4E)-3-methylhex-4-en-2-one"
    assert smiles_to_iupac("CC(=O)[C@@H](C)/C=C/C") == "(3S,4E)-3-methylhex-4-en-2-one"


def test_two_ring_aromatic_substituent_ketone_ring_bond():
    assert smiles_to_iupac("O=C1C=CCCC1c1ccccc1") == "6-phenylcyclohex-2-en-1-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COOCC", "(methylperoxy)ethane"),
    ],
)
def test_peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
    ],
)
def test_ring_not_supported_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)COOC(C)(C)C", "1-(tert-butylperoxy)-2-methylpropane"),
        ("CC(C)COOC(C)CC", "2-[(2-methylpropyl)peroxy]butane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)OOCC", "(2S)-2-(ethylperoxy)butane"),
    ],
)
def test_stereocenter_on_parent_chain__peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1OOC", "(methylperoxy)benzene"),
        ("c1ccccc1OOC(C)C", "[(propan-2-yl)peroxy]benzene"),
        ("c1ccccc1COOCC", "[(ethylperoxy)methyl]benzene"),
        ("c1ccccc1COOC(C)C", "{[(propan-2-yl)peroxy]methyl}benzene"),
    ],
)
def test_benzene_ring_parent__peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected



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
        pytest.param("NNS(=O)(=O)O", "hydrazinesulfonic acid", id="bare_sulfonic_acid"),
        pytest.param("CNNS(=O)(=O)O", "2-methylhydrazine-1-sulfonic acid", id="substituted_sulfonic_acid"),
        pytest.param("NNS(=O)(=S)O", "hydrazinesulfonothioic O-acid", id="functional_replacement"),
        pytest.param("NNC(=O)[O-]", "hydrazinecarboxylate", id="anion"),
        pytest.param(
            "C[N+](C)(C)NS(=O)(=O)[O-]", "1,1,1-trimethylhydrazin-1-ium-2-sulfonate", id="zwitterion_ionic_centre_first"
        ),
        pytest.param("NNS(=O)(=O)CC(=O)O", "(hydrazinesulfonyl)acetic acid", id="prefix_sulfonyl"),
        pytest.param("CNNS(=O)CC(=O)O", "(2-methylhydrazine-1-sulfinyl)acetic acid", id="prefix_substituted_sulfinyl"),
    ],
)
def test_hydrazine_parent_with_acid_group_suffix_or_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cationic_hydrazine_outranks_a_neutral_acid():
    # P-41: a cation is senior to an acid, so the carboxy group is a prefix
    assert smiles_to_iupac("C[N+](C)(C)NC(=O)O") == "2-carboxy-1,1,1-trimethylhydrazin-1-ium"


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
        ("CC(=O)OB1CCCCC1", "borinan-1-yl acetate"),
        ("O=C(c1ccccc1)OP1CCCCC1", "phosphinan-1-yl benzoate"),
        ("CC(=O)ON1CCCC1", "pyrrolidin-1-yl acetate"),
        ("CC(=O)ON(C)C", "1-[(dimethylamino)oxy]ethan-1-one"),
        ("CCCC(=O)ONC", "1-[(methylamino)oxy]butan-1-one"),
        ("S(C#N)CCC(=O)O", "3-(thiocyanato)propanoic acid"),
        ("S(C#N)CCC(SCC)=O", "S-ethyl 3-(thiocyanato)propanethioate"),
        ("OC(=O)CCOC#N", "3-cyanatopropanoic acid"),
        ("OC(=O)CCC[Se]C#N", "4-(selenocyanato)butanoic acid"),
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
        ("CC(=O)OC(=O)CC(=O)OC(=O)CCC(=O)OC(C)=O", "acetic 3-(acetyloxy)-3-oxopropanoic butanedioic dianhydride"),
        (
            "CCC(=O)OC(=O)C(C)CC(=O)OC(=O)CCC(=O)OC(C)=O",
            "4-[4-(acetyloxy)-4-oxobutanoic] 2-methylbutanedioic 1-propanoic dianhydride",
        ),
        ("CC(=O)OC(=O)c1cc(C(=O)OC(C)=O)c(C(=O)OC(=O)CC)cc1", "2,4-diacetic 1-propanoic benzene-1,2,4-tricarboxylic trianhydride"),
        ("O=C(Cl)c1ccc(C(=O)Cl)cc1", "benzene-1,4-dicarbonyl dichloride"),
        ("NC(=O)C(=O)Br", "oxamoyl bromide"),
        ("ClC(=S)C1CCCCC1", "cyclohexanecarbothioyl chloride"),
        ("CCCC(=O)C#N", "butanoyl cyanide"),
        ("O=C(N=C=S)C(=O)N=C=S", "oxalyl diisothiocyanate"),
        ("CC(=O)OP(=O)(O)O", "(acetyloxy)phosphonic acid"),
        ("O=COP(O)O", "(formyloxy)phosphonous acid"),
        ("CC(=O)OP(C)O", "(acetyloxy)(methyl)phosphinous acid"),
        ("CCC(=O)OB(O)O", "(propanoyloxy)boronic acid"),
        ("CC(=O)OP(=O)(OC)OC", "acetic (dimethyl hydrogen phosphate) anhydride"),
        ("CC(=O)O[As](C)(C)=O", "acetic dimethylarsinic anhydride"),
        ("CO[Sb](OC)OC", "trimethyl stiborite"),
        ("ClC(=O)Cl", "carbonyl dichloride"),
        ("N#CC(=O)Cl", "carbonocyanidoyl chloride"),
        ("N#CCl", "carbononitridic chloride"),
        ("ClC(=O)OC(=O)Cl", "dicarbonic dichloride"),
        ("C(=O)(NC(=O)Cl)Cl", "2-imidodicarbonic dichloride"),
        ("C(N)(=O)Cl", "carbonochloridic amide"),
        ("CNC(=O)Cl", "N-methylcarbonochloridic amide"),
        ("C(C)OC(=O)OC(C(=O)OCC)C(C(C)(C)C)=O", "ethyl 2-[(ethoxycarbonyl)oxy]-4,4-dimethyl-3-oxopentanoate"),
        ("NC(CC1=CC=C(O1)C(=O)O)=O", "5-(2-amino-2-oxoethyl)furan-2-carboxylic acid"),
        ("C(C)(=S)C(C(C)=S)C(C)=S", "3-(ethanethioyl)pentane-2,4-dithione"),
        ("S(=S)(=S)(S)C1=C(C(=O)O)C=CC=C1", "2-(trithiosulfo)benzoic acid"),
        ("CC(=O)OS(=O)(=O)CS(=O)(=O)OC(=O)CC", "acetic methanedisulfonic propanoic dianhydride"),
        ("CC(=O)SC(=O)CCC(=O)SC(C)=O", "diacetic butanedioic bis(thioanhydride)"),
        ("CC(=O)SS(=O)(=O)CS(=O)(=O)SC(C)=O", "diacetic methanedisulfonic bis(thioanhydride)"),
        ("CC(=O)OC(=O)CCCC(=O)SC(C)=O", "acetic 5-(acetylsulfanyl)-5-oxopentanoic anhydride"),
        ("CC(=O)OC(=O)CCC(=O)OC(=O)CCC(=O)SC(=O)CC", "acetic butanedioic 4-oxo-4-(propanoylsulfanyl)butanoic dianhydride"),
        ("CC(=S)OC(=O)CCC(=O)OC(C)=O", "acetic butanedioic ethanethioic dianhydride"),
        ("CC(=O)OC(=S)CCC(=S)OC(C)=O", "diacetic butanebis(thioic) dianhydride"),
        ("CC(=O)OP(=O)(OC(C)=O)OC(C)=O", "triacetic phosphoric trianhydride"),
        ("CC(=O)OP(=O)(OC(C)=O)P(=O)(OC(C)=O)OC(C)=O", "tetraacetic hypodiphosphoric tetraanhydride"),
        ("CC(=O)OS(=O)S(=O)OC(=O)CC", "acetic dithionous propanoic dianhydride"),
        ("CC(=O)OS(=O)(=O)OC(C)=O", "diacetic sulfuric dianhydride"),
        ("CC(=O)OS(=O)(=O)O", "acetic sulfuric monoanhydride"),
        ("CC(=O)OP(=O)(OC(C)=O)CP(=O)(OC(C)=O)OC(C)=O", "tetraacetic methylenebis(phosphonic) tetraanhydride"),
        ("CC(=O)OB(OC(C)=O)OC(C)=O", "triacetic boric trianhydride"),
        ("CC(=O)OP(=O)(OC(=O)CC)OC(=O)CC", "acetic dipropanoic phosphoric trianhydride"),
        (
            "O=C(OC(C)=O)CCC(=O)OP(=O)(OC(=O)CCC(=O)OC(C)=O)OC(=O)CCCCC(=O)OC(=O)CCC",
            "bis[4-(acetyloxy)-4-oxobutanoic] 6-(butanoyloxy)-6-oxohexanoic phosphoric trianhydride",
        ),
        ("OS(=O)(=O)CS(=O)(=O)O", "methanedisulfonic acid"),
        ("OS(=O)(=O)C(C)S(=O)(=O)O", "ethane-1,1-disulfonic acid"),
        ("CC(=O)OOOC(C)=O", "1,1'-trioxidanediyldi(ethan-1-one)"),
        ("CC(=O)SSSSC(C)=O", "1,1'-tetrasulfanediyldi(ethan-1-one)"),
        ("CC(=O)SOSC(C)=O", "1,1'-dithioxanediyldi(ethan-1-one)"),
        ("CCC(=O)SSSSC(C)=O", "1-(acetyltetrasulfanyl)propan-1-one"),
        ("CC(=O)OOSC(C)=O", "1-[(acetylperoxy)sulfanyl]ethan-1-one"),
        ("CC(=O)SOOC(C)=O", "1-[(acetylperoxy)sulfanyl]ethan-1-one"),
        ("CC(=O)O[Te][Se]SC(C)=O", "3-oxa-6-thia-5-selena-4-telluraoctane-2,7-dione"),
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
        ("N#COO", "peroxycyanic acid"),
        ("N#COS", "thioperoxycyanic OS-acid"),
        ("NC(=N)S", "carbamimidothioic acid"),
        ("[N-]=[N+]=NC(=O)O", "carbonazidic acid"),
        ("OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)O", "3,5,7-trioxo-2,4,6,8-tetraoxanonanedioic acid"),
        ("OC(=S)OC(=S)O", "1,3-dithiodicarbonic O1,O3-acid"),
        ("OC(=O)SC(=O)O", "2-thiodicarbonic acid"),
        ("C(=O)(SC(=O)O)N", "1-amido-2-thiodicarbonic acid"),
        ("C(=O)(NC(=O)NC(=O)O)NN", "1-hydrazido-2,4-diimidotricarbonic acid"),
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
        ("CCCC[Sn](CCCC)(OC(=O)CCC(=O)OC)OC(=O)CCC(=O)OC", "dimethyl dibutylstannanediyl dibutanedioate"),
        ("COC(=O)CCC(=O)O[Sn](CCCC)(CCCC)OC(=O)CCC(=O)OC", "dimethyl dibutylstannanediyl dibutanedioate"),
        ("CC(=O)Oc1ccc(cc1)C(=O)OC", "methyl 4-(acetyloxy)benzoate"),
        ("O=C(OCCOC(C)=O)CC(=O)OCCOC(C)=O", "bis[2-(acetyloxy)ethyl] propanedioate"),
        (
            "COC(=O)c1cccc(OC(=O)CCC(=O)OCC(=O)OC)c1",
            "3-(methoxycarbonyl)phenyl 2-methoxy-2-oxoethyl butanedioate",
        ),
        (
            "C(CCC(=O)OC1=C(C=CC=C1)OC(CCC(=O)OC1=C(C=CC=C1)OC(CCC(=O)OC)=O)=O)(=O)OC",
            "dimethyl butanedioylbis(oxy-2,1-phenylene)dibutanedioate",
        ),
    ],
)
def test_polyester_principal_acid_and_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CCCCCOCCCCCCC", "1-(pentyloxy)heptane"),
        ("CC(C)(C)COCCCCC", "1-(2,2-dimethylpropoxy)pentane"),
        ("CCCCCOCCO", "2-(pentyloxy)ethan-1-ol"),
        ("OCCOCC1CCCCC1", "2-(cyclohexylmethoxy)ethan-1-ol"),
        ("CC(C)OCC=O", "[(propan-2-yl)oxy]acetaldehyde"),
        ("c1ccccc1OCCCCC", "(pentyloxy)benzene"),
        ("CC(C)(C)OCCCCC", "1-tert-butoxypentane"),
        ("CCOCCCCCCC", "1-ethoxyheptane"),
    ],
)
def test_alkoxy_prefix_enclosing_marks(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("COOC", "(methylperoxy)methane"),
        ("CCOOCC", "(ethylperoxy)ethane"),
        ("CC(C)OOC", "2-(methylperoxy)propane"),
        ("c1ccccc1OOCC", "(ethylperoxy)benzene"),
    ],
)
def test_peroxy_prefix_enclosing_marks(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CCCCOC(OCCCC)CC", "1,1-dibutoxypropane"),
        ("CCCCCOC(OCCCCC)CC", "1,1-bis(pentyloxy)propane"),
        ("CCCCOC(OCC)CC", "1-butoxy-1-ethoxypropane"),
        ("CCCCCOC(OC)CCCCC", "1-methoxy-1-(pentyloxy)hexane"),
        ("CC(C)(C)OC(C)OC", "1-tert-butoxy-1-methoxyethane"),
        ("CO[C@H](OCC)CC", "(1R)-1-ethoxy-1-methoxypropane"),
    ],
)
def test_acetal_parent_chain_is_the_acetal_carbons_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


_ACETAL_STEMS = {"tert-but": 4, "meth": 1, "eth": 2, "prop": 3, "but": 4, "pent": 5, "hex": 6, "hept": 7}
_ACETAL_STEM_RE = re.compile(r"(bis\(|di|tri)?(tert-but|meth|eth|prop|but|pent|hex|hept)(?=an|yl|oxy|ane|-)")


def _carbons_in_acyclic_name(name):
    total = 0
    for m in _ACETAL_STEM_RE.finditer(name):
        times = {None: 1, "bis(": 2, "di": 2, "tri": 3}[m.group(1)]
        total += times * _ACETAL_STEMS[m.group(2)]
    return total


@pytest.mark.parametrize("alkoxy_a", ["C", "CC", "CCCC", "CCCCC", "CC(C)(C)"])
@pytest.mark.parametrize("alkoxy_b", ["C", "CC", "CCCC", "CCCCC"])
@pytest.mark.parametrize("acyl", ["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"])
def test_acetal_name_accounts_for_every_carbon(alkoxy_a, alkoxy_b, acyl):
    from rdkit import Chem

    smiles = f"{alkoxy_a}OC({acyl})O{alkoxy_b}"
    carbons = sum(1 for atom in Chem.MolFromSmiles(smiles).GetAtoms() if atom.GetAtomicNum() == 6)
    assert _carbons_in_acyclic_name(smiles_to_iupac(smiles)) == carbons


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=C=O", "ethenone", id="ketene_P-64.2.1.2"),
        pytest.param("BrC(Br)=C=O", "dibromoethenone", id="locants_omitted_P-14.3.4.5"),
        pytest.param("ClC=C=O", "chloroethenone", id="locants_omitted_P-14.3.4.4"),
        pytest.param("CCCCC(=C=O)CCCC", "2-butylhex-1-en-1-one", id="chain_through_ketene_carbon"),
        pytest.param("O=C=C1CCCCC1", "cyclohexylidenemethanone", id="ring_ylidene_methanone"),
        pytest.param("CC=C=O", "prop-1-en-1-one", id="longer_ketene"),
        pytest.param("O=C=C=C", "propa-1,2-dien-1-one", id="cumulated_diene"),
    ],
)
def test_ketenes_are_named_as_ene_ones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles", ["O=C=CC(=O)O", "O=C=CC(N)C", "O=C=CCC(C)=O"])
def test_ketene_with_other_groups_is_rejected_not_misnamed(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=C(O)c1ccccc1S(=O)(=O)c1ccccc1C(=O)O", "2,2'-sulfonyldibenzoic acid"),
        ("O=C([O-])c1ccccc1S(=O)(=O)c1ccccc1C(=O)O.[Na+]", "sodium 2-[(2-carboxyphenyl)sulfonyl]benzoate"),
        ("O=C(O)c1ccccc1S(=O)(=O)c1ccccc1C(N)=O", "2-[(2-carbamoylphenyl)sulfonyl]benzoic acid"),
    ],
)
def test_diaryl_sulfone_with_carboxylic_groups_keeps_the_sulfonyl_bridge(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COOC1CCCCC1", "(methylperoxy)cyclohexane", id="peroxy_prefix_on_a_ring"),
        pytest.param("COOCOC", "methoxy(methylperoxy)methane", id="peroxy_and_alkoxy_on_methane"),
        pytest.param("C=COOC", "(methylperoxy)ethene", id="peroxy_prefix_on_a_double_bond"),
    ],
)
def test_peroxy_prefixes_without_a_principal_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=O", "formaldehyde"),
        ("C#N", "formonitrile"),
        ("C(C)=S", "ethanethial"),
        ("CCC=S", "propanethial"),
        ("CCC=[Se]", "propaneselenal"),
        ("CCC=[Te]", "propanetellanal"),
        ("C1(=CC=CC=C1)C=S", "benzenecarbothialdehyde"),
        ("C(CCCCC)=[Se]", "hexaneselenal"),
        ("CC=[Te]", "ethanetellanal"),
        ("C(CCCC=S)=S", "pentanedithial"),
        ("C=CC=S", "prop-2-enethial"),
        ("S=C1CCC(CC1)C=[Se]", "4-sulfanylidenecyclohexane-1-carboselenaldehyde"),
        ("OC(=O)c1ccc(C=S)cc1", "4-(methanethioyl)benzoic acid"),
        ("OC(=O)C1CCC(CC1)C=[Se]", "4-(methaneselenoyl)cyclohexane-1-carboxylic acid"),
        ("N#CN1CCCCC1", "piperidine-1-carbonitrile"),
        ("N#CN1CCOCC1", "morpholine-4-carbonitrile"),
        ("N#Cn1cccc1", "1H-pyrrole-1-carbonitrile"),
        ("N#CN1CCCCC1C", "2-methylpiperidine-1-carbonitrile"),
    ],
)
def test_chalcogen_aldehydes_and_ring_nitrogen_nitriles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("S=CCCC=[Se]", "4-selanylidenebutanethial"),
        ("S=CC1CCC(C(N)=O)CC1", "4-(methanethioyl)cyclohexane-1-carboxamide"),
    ],
)
def test_mixed_chalcogen_aldehydes_and_methanethioyl_beside_an_amide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C(C)(=O)OC=1C=NC(=CC1)C1=CC=C(C=C1)OC(C)=O", "6-[4-(acetyloxy)phenyl]pyridin-3-yl acetate"),
        ("C(C)(=O)OCCC1=C(C=CC=C1)COC(C)=O", "2-{2-[(acetyloxy)methyl]phenyl}ethyl acetate"),
        ("OC(=O)OOC(=O)O", "2-peroxydicarbonic acid"),
        ("[O-]C(=S)OC(=O)[O-]", "1-thiodicarbonate"),
        ("O=C(CCC(=O)N=C=S)[N+]#[C-]", "butanedioyl isocyanide isothiocyanate"),
        ("CCC(=N)[Se]c1ccccc1C(O)=N", "2-(propanimidoylselanyl)benzene-1-carboximidic acid"),
    ],
)
def test_polyester_seniority_peroxydicarbonic_and_acyl_pseudohalide_gaps(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC(=O)CC1CCCCC1C(=O)O", "2-(2-methoxy-2-oxoethyl)cyclohexane-1-carboxylic acid"),
        ("NC(=O)CC1CCCCC1C(=O)O", "2-(2-amino-2-oxoethyl)cyclohexane-1-carboxylic acid"),
        ("CCOC(=O)CC(Cc1ccccc1)C(=O)O", "2-benzyl-4-ethoxy-4-oxobutanoic acid"),
    ],
)
def test_ester_end_carbon_joins_the_substituent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OCCCC1(C(CCCC1)CCCOC=O)CCCOC(C)=O", "2-[3-(formyloxy)propyl]cyclohexane-1,1-diyldi(propane-3,1-diyl) diacetate"),
        ("CC(=O)OCc1ccc(COC(C)=O)cc1", "1,4-phenylenebis(methylene) diacetate"),
    ],
)
def test_diesters_of_a_ring_joined_to_identical_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCN(CC)C(=S)S[Sb](C(C)C)C(C)C", "diethylcarbamothioic di(propan-2-yl)stibinous thioanhydride"),
        ("O=C(c1ccccc1)OI", "benzoic hypoiodous anhydride"),
        ("CB(C)OOB(C)C", "dimethylborinic peroxyanhydride"),
    ],
)
def test_anhydrides_of_stibinous_halogen_and_boron_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-63.2.2.1.1, P-29.6.1: no multiple bond or heteroatom of a substituent group may be lost from its prefix
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=COCCO", "2-(ethenyloxy)ethan-1-ol", id="unsaturated_alkoxy_group_on_an_alkanol"),
        pytest.param("C#COC(C)O", "1-(ethynyloxy)ethan-1-ol", id="triple_bond_in_an_alkoxy_group"),
        pytest.param("CC(O)(O)P=O", "(1,1-dihydroxyethyl)phosphanone", id="tert_butyl_needs_three_carbons"),
    ],
)
def test_substituent_prefixes_keep_multiple_bonds_and_heteroatoms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
