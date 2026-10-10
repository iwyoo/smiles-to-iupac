import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._coexisting_groups import name_via_senior_acyclic
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._seniority import SUFFIX_CLASS_RANK


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NCC=CCO", "4-aminobut-2-en-1-ol", id="unsaturated_chain"),
        pytest.param("C[C@H](N)CO", "(2S)-2-aminopropan-1-ol", id="specified_stereocenter"),
        pytest.param("NCC(O)COCC", "1-amino-3-ethoxypropan-2-ol", id="ether_coexisting"),
    ],
)
def test_unsaturated_chain_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C",
            "5α-cholestan-3β-ol",
        ),
    ],
)
def test_steroid_alcohol_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCC=O", "aminoacetaldehyde"),
    ],
)
def test_smiles_to_iupac_aldehyde_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC1CCC(C=O)CC1", "4-aminocyclohexane-1-carbaldehyde", id="ring"),
        pytest.param("NCC=CC=O", "4-aminobut-2-enal", id="unsaturated_chain__aldehyde_amine"),
        pytest.param("C[C@H](N)C=O", "(2S)-2-aminopropanal", id="specified_stereocenter__aldehyde_amine"),
        pytest.param("NC(CO)C=O", "2-amino-3-hydroxypropanal", id="hydroxyl_coexisting"),
        pytest.param("O=CC(C)C(=O)O", "2-methyl-3-oxopropanoic acid", id="branch_substituent"),
        pytest.param("CC(=O)C(Cl)C=O", "2-chloro-3-oxobutanal", id="halogen_substituent"),
        pytest.param("O=CC(=O)CC=O", "oxobutanedial", id="two_aldehydes_with_ketone"),
        pytest.param("C=CC(=O)CC=O", "3-oxopent-4-enal", id="unsaturated_chain__aldehyde_ketone"),
        pytest.param("O=CC1CCC(=O)C1", "3-oxocyclopentane-1-carbaldehyde", id="ring__aldehyde_ketone"),
        pytest.param("OCC(=O)CC=O", "4-hydroxy-3-oxobutanal", id="hydroxyl_coexistence"),
    ],
)
def test_ring_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_aldehyde_ketone():
    assert smiles_to_iupac("c1ccccc1CC(=O)CC=O") == "3-oxo-4-phenylbutanal"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C=O") == "2-oxo-4-phenylbutanal"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Cc1ccccc1CC(=O)CC=O", "4-(2-methylphenyl)-3-oxobutanal", id="substituted_benzene_ring_aldehyde_ketone"),
        pytest.param("C=Cc1ccccc1CC(=O)CC=O", "4-(2-ethenylphenyl)-3-oxobutanal", id="chain_aldehyde_ketone_unsaturation"),
    ],
)
def test_phenyl_substituted_benzene_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCC(N)=O", "3-aminopropanamide"),
    ],
)
def test_smiles_to_iupac_amide_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NCC(N)C(N)=O", "2,3-diaminopropanamide", id="two_amines"),
        pytest.param("NC1CCCCC1C(N)=O", "2-aminocyclohexane-1-carboxamide", id="ring__amide_amine"),
        pytest.param("NCC=CC(N)=O", "4-aminobut-2-enamide", id="unsaturated_chain__amide_amine"),
        pytest.param("N[C@@H](C)C(N)=O", "(2S)-2-aminopropanamide", id="specified_stereocenter__amide_amine"),
        pytest.param("NCC(O)C(N)=O", "3-amino-2-hydroxypropanamide", id="hydroxyl_coexisting__amide_amine"),
        pytest.param("C[N-][NH+](C)C", "1,2,2-trimethylhydrazin-2-ium-1-ide", id="trimethylhydrazinium_ide_amine_imide"),
        pytest.param("OC(=O)CC(=O)NC", "3-(methylamino)-3-oxopropanoic acid", id="n_substituted_amide"),
        pytest.param("OC(=O)C1CCCCC1C(=O)N", "2-carbamoylcyclohexane-1-carboxylic acid", id="amide_on_ring"),
    ],
)
def test_two_amines_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)(C(=O)O)N", "2-amino-2-methylpropanoic acid"),
    ],
)
def test_smiles_to_iupac_carboxylic_acid_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_carboxylic_acid_amine():
    assert smiles_to_iupac("c1ccccc1CC(N)C(=O)O") == "phenylalanine"
    assert smiles_to_iupac("c1ccccc1CCC(N)C(=O)O") == "2-amino-4-phenylbutanoic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccccc1C(N)C(=O)O", "amino(phenyl)acetic acid", id="directly_attached_to_amine_carbon"),
        pytest.param("Cc1ccccc1CC(N)C(=O)O", "2-amino-3-(2-methylphenyl)propanoic acid", id="substituted_benzene_ring_carboxylic_acid_amine"),
        pytest.param("C=Cc1ccccc1CC(N)C(=O)O", "2-amino-3-(2-ethenylphenyl)propanoic acid", id="chain_carboxylic_acid_amine_unsaturation"),
    ],
)
def test_phenyl_directly_attached_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)C[Se](=O)O", "seleninoacetic acid"),
    ],
)
def test_carboxylic_acid_seleninic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_seleninic_acid_still_routes_normally():
    assert smiles_to_iupac("C[Se](=O)O") == "methaneseleninic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)S(=O)O", "sulfinoformic acid"),
    ],
)
def test_carboxylic_acid_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)S(=O)(=O)N", "sulfamoylformic acid"),
    ],
)
def test_carboxylic_acid_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("OC(=O)CC(S(=O)(=O)N)C(=O)O", "sulfamoylbutanedioic acid", id="multiple_carboxylic_acids"),
        pytest.param("OC(=O)C(O)CS(=O)(=O)N", "2-hydroxy-3-sulfamoylpropanoic acid", id="other_heteroatom"),
        pytest.param("OC(=O)C=CCS(=O)(=O)N", "4-sulfamoylbut-2-enoic acid", id="unsaturated_chain__carboxylic_acid_sulfonamide"),
        pytest.param("OC(=O)C1CCC(S(=O)(=O)N)CC1", "4-sulfamoylcyclohexane-1-carboxylic acid", id="ring__carboxylic_acid_sulfonamide"),
        pytest.param("OC(=O)CC(S(=O)(=O)O)C(=O)O", "sulfobutanedioic acid", id="multiple_carboxylic_acids__carboxylic_acid_sulfonic_acid"),
        pytest.param("OC(=O)C(O)CS(=O)(=O)O", "2-hydroxy-3-sulfopropanoic acid", id="other_heteroatom__carboxylic_acid_sulfonic_acid"),
        pytest.param("OC(=O)C=CCS(=O)(=O)O", "4-sulfobut-2-enoic acid", id="unsaturated_chain__carboxylic_acid_sulfonic_acid"),
        pytest.param("OC(=O)C1CCC(S(=O)(=O)O)CC1", "4-sulfocyclohexane-1-carboxylic acid", id="ring__carboxylic_acid_sulfonic_acid"),
    ],
)
def test_multiple_carboxylic_acids_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC(S(=O)(=O)O)C(=O)O", "3-phenyl-2-sulfopropanoic acid"),  # CID 129643411
    ],
)
def test_phenyl_chain_carboxylic_acid_sulfonic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_name_via_senior_acyclic_rejects_reversed_seniority():
    with pytest.raises(AssertionError):
        name_via_senior_acyclic(lambda *a, **k: "unused", "alcohol", "sulfonic_acid", (), {})


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC(=O)OCC(COC(C)=O)(COC(C)=O)COC(C)=O",
            "2,2-bis[(acetyloxy)methyl]propane-1,3-diyl diacetate",
        ),
        ("CC(=O)OCC=CCOC(C)=O", "but-2-ene-1,4-diyl diacetate"),
    ],
)
def test_diester_acyloxy_generalized_backbones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC(=O)Oc1ccc(Cl)cc1OC(=O)c1ccc(Cl)cc1",
            "4-chloro-1,2-phenylene 1-acetate 2-(4-chlorobenzoate)",
        ),
    ],
)
def test_ring_diyl_differing_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(OC1CCCCC1OC(=O)C1CCCCC1)C1CCCCC1", "cyclohexane-1,2-diyl di(cyclohexanecarboxylate)"),
    ],
)
def test_ring_diyl_substituted_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCC(=O)OCCCOC(=O)CCl", "propane-1,3-diyl bis(chloroacetate)"),
    ],
)
def test_acyclic_diyl_with_generalized_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1ccc(N(=O)=O)cc1OC(C)=O", "4-nitro-1,2-phenylene diacetate"),
        ("CC(=O)Oc1ccc(N(C)CC)cc1OC(C)=O", "4-[ethyl(methyl)amino]-1,2-phenylene diacetate"),
        ("CC(=O)Oc1ccc(SC)cc1OC(C)=O", "4-(methylsulfanyl)-1,2-phenylene diacetate"),
        ("CC(=O)Oc1ccc(Oc2ccccc2)cc1OC(C)=O", "4-phenoxy-1,2-phenylene diacetate"),
    ],
)
def test_functional_substituents_on_diyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(OCCOC(=O)c1ccncc1)c1ccncc1", "ethane-1,2-diyl di(pyridine-4-carboxylate)"),
        ("CC(=O)OCCOC(=O)CC#N", "ethane-1,2-diyl acetate cyanoacetate"),
        ("CC(=O)OC1CCCCC1OC(=O)CC(C)=O", "cyclohexane-1,2-diyl acetate 3-oxobutanoate"),
        ("CC(=O)Oc1ccccc1OC(=O)/C=C/C", "1,2-phenylene acetate (2E)-but-2-enoate"),
        (
            "CC(=O)OC1CCCC1OC(=O)[C@H](C)Cl",
            "cyclopentane-1,2-diyl acetate (2S)-2-chloropropanoate",
        ),
    ],
)
def test_generalized_acyl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OC1CCC(COC(C)=O)CC1", "4-[(acetyloxy)methyl]cyclohexyl acetate"),
    ],
)
def test_ring_and_chain_polyol_uses_acyloxy_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)OCC(=C)COC(C)=O",
    ],
)
def test_ring_diyl_out_of_scope(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccsc1Cc1ccsc1", "2-[(thiophen-3-yl)methyl]thiophene"),
    ],
)
def test_disjoint_ring_pair_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC1CCCCC1CC1CCCCC1", "1-(cyclohexylmethyl)-2-methylcyclohexane", id="ring_with_extra_substituent"),
        pytest.param("CC(=O)OCCN", "2-aminoethyl acetate", id="amine_on_alcohol_part"),
        pytest.param("NCC(N)C(=O)OC", "methyl 2,3-diaminopropanoate", id="two_amines__ester_amine"),
        pytest.param("NCC(=O)OCOC(=O)C", "methylene acetate aminoacetate", id="two_esters_names_polyester"),
        pytest.param("NC1CCCCC1C(=O)OC", "methyl 2-aminocyclohexane-1-carboxylate", id="ring__ester_amine"),
        pytest.param("NCC=CC(=O)OC", "methyl 4-aminobut-2-enoate", id="unsaturated_chain__ester_amine"),
    ],
)
def test_ring_with_extra_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NCC(O)C(=O)OC", "methyl 3-amino-2-hydroxypropanoate", id="hydroxyl_coexisting__ester_amine"),
        pytest.param("CC(C)OCC=O", "[(propan-2-yl)oxy]acetaldehyde", id="branched_alkoxy_r_prime"),
        pytest.param("COCC(Cl)C=O", "2-chloro-3-methoxypropanal", id="halogen_on_main_chain_still_works"),
        pytest.param("COC[C@@H](C)C=O", "(2R)-3-methoxy-2-methylpropanal", id="specified_stereocenter__ether_aldehyde"),
        pytest.param("CC(C)OCC(N)=O", "2-[(propan-2-yl)oxy]acetamide", id="branched_alkoxy_r_prime__ether_amide"),
        pytest.param("COCC(Cl)C(N)=O", "2-chloro-3-methoxypropanamide", id="halogen_on_main_chain_still_works__ether_amide"),
        pytest.param("NC(=O)C(COC)C(N)=O", "2-(methoxymethyl)propanediamide", id="two_amides"),
        pytest.param("COCC(OC)C(N)=O", "2,3-dimethoxypropanamide", id="two_ethers"),
        pytest.param("NC(=O)C1CCCCC1COC", "2-(methoxymethyl)cyclohexane-1-carboxamide", id="ring__ether_amide"),
        pytest.param("COC[C@@H](C)C(N)=O", "(2R)-3-methoxy-2-methylpropanamide", id="specified_stereocenter__ether_amide"),
        pytest.param("CC(C)OCCN", "2-[(propan-2-yl)oxy]ethan-1-amine", id="branched_alkoxy_r_prime__ether_amine"),
        pytest.param("COCC(Cl)CN", "2-chloro-3-methoxypropan-1-amine", id="halogen_on_main_chain_still_works__ether_amine"),
        pytest.param("NCC(N)COC", "3-methoxypropane-1,2-diamine", id="two_amines__ether_amine"),
        pytest.param("NCC=CCOC", "4-methoxybut-2-en-1-amine", id="unsaturated_chain__ether_amine"),
        pytest.param("N[C@@H](C)COC", "(2S)-1-methoxypropan-2-amine", id="specified_stereocenter__ether_amine"),
        pytest.param("CC(C)OCC(=O)OC", "methyl [(propan-2-yl)oxy]acetate", id="branched_alkoxy_r_prime__ether_ester"),
        pytest.param("COCC(Cl)C(=O)OC", "methyl 2-chloro-3-methoxypropanoate", id="halogen_on_acyl_chain_still_works"),
        pytest.param("CC(=O)OCCOC", "2-methoxyethyl acetate", id="ether_on_alcohol_part"),
        pytest.param("O=C(OC)C1CCCCC1COC", "methyl 2-(methoxymethyl)cyclohexane-1-carboxylate", id="ring__ether_ester"),
        pytest.param("COC[C@@H](C)C(=O)OC", "methyl (2R)-3-methoxy-2-methylpropanoate", id="specified_stereocenter__ether_ester"),
        pytest.param("CC(C)OCCOO", "2-[(propan-2-yl)oxy]ethane-1-peroxol", id="branched_alkoxy_r_prime__ether_hydroperoxide"),
        pytest.param("COCC(Cl)COO", "2-chloro-3-methoxypropane-1-peroxol", id="halogen_on_main_chain_still_works__ether_hydroperoxide"),
    ],
)
def test_hydroxyl_coexisting__and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COCC(OC)COO", "2,3-dimethoxypropane-1-peroxol", id="two_ethers"),
        pytest.param("OOCC=CCOC", "4-methoxybut-2-ene-1-peroxol", id="unsaturated_chain"),
        pytest.param("OO[C@@H](C)COC", "(2S)-1-methoxypropane-2-peroxol", id="specified_stereocenter_ether_hydroperoxide"),
    ],
)
def test_peroxol_with_ether_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)OCC(C)=O", "1-[(propan-2-yl)oxy]propan-2-one", id="branched_alkoxy_r_prime__ether_ketone"),
        pytest.param("O=CC(=O)COC", "3-methoxy-2-oxopropanal", id="two_ketones"),
        pytest.param("O=C1CCCCC1COC", "2-(methoxymethyl)cyclohexan-1-one", id="ring__ether_ketone"),
        pytest.param("COC[C@@H](C)C(C)=O", "(3R)-4-methoxy-3-methylbutan-2-one", id="specified_stereocenter__ether_ketone"),
        pytest.param("CC(C)OCCS", "2-[(propan-2-yl)oxy]ethane-1-thiol", id="branched_alkoxy_r_prime__ether_thiol"),
        pytest.param("COCC(Cl)CS", "2-chloro-3-methoxypropane-1-thiol", id="halogen_on_main_chain_still_works__ether_thiol"),
        pytest.param("SCC=CCOC", "4-methoxybut-2-ene-1-thiol", id="unsaturated_chain__ether_thiol"),
        pytest.param("S[C@@H](C)COC", "(2S)-1-methoxypropane-2-thiol", id="specified_stereocenter__ether_thiol"),
        pytest.param("ClCC(=O)N1CCCCC1", "2-chloro-1-(piperidin-1-yl)ethan-1-one", id="chloroacetylpiperidine"),
        pytest.param("PC(=O)CCC", "1-phosphanylbutan-1-one", id="acyl_on_group_15_hydride_is_a_pseudoketone"),
        pytest.param("CC(=O)[Si](C)(C)C", "1-(trimethylsilyl)ethan-1-one", id="acyl_on_group_14_atom_is_a_pseudoketone"),
        pytest.param("[SiH3]C(=O)CC(=O)O", "3-oxo-3-silylpropanoic acid", id="senior_acid_keeps_the_pseudoketone_as_a_prefix"),
        pytest.param("O=C(CC)[SiH2]C(=O)CC", "dipropanoylsilane", id="two_acyl_groups_on_one_heteroatom_keep_the_acyl_name"),
        pytest.param("CC(=O)P(=O)(C)C", "acetyldi(methyl)-λ5-phosphanone", id="phosphoryl_acyl_is_not_a_pseudoketone"),
        pytest.param("CNC(=O)c1cc(C)on1", "N,5-dimethyl-1,2-oxazole-3-carboxamide", id="n_substituted_amide_on_heteroaromatic_monocycle"),
        pytest.param("O=C(NC1CCCCC1)c1ccc2ncccc2c1", "N-cyclohexylquinoline-6-carboxamide", id="n_substituted_amide_on_fused_heteroaromatic"),
        pytest.param("COc1ccc(C(=O)N2CCCCC2)cc1", "(4-methoxyphenyl)(piperidin-1-yl)methanone", id="aroyl_ring_nitrogen_pseudoketone"),
        pytest.param("CC(=O)N1CCN(c2ccccc2)CC1", "1-(4-phenylpiperazin-1-yl)ethan-1-one", id="acyl_ring_nitrogen_with_aryl_on_other_nitrogen"),
        pytest.param("CC(C)C(=O)N1CCCC1C", "2-methyl-1-(2-methylpyrrolidin-1-yl)propan-1-one", id="branched_acyl_substituted_ring"),
    ],
)
def test_branched_alkoxy_r_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCOO", "2-aminoethane-1-peroxol"),
    ],
)
def test_smiles_to_iupac_hydroperoxide_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC(=O)COOC", "2-(methylperoxy)acetamide", id="peroxy_prefix_on_amide"),
        pytest.param("NCCSSC", "2-(methyldisulfanyl)ethan-1-amine", id="disulfanyl_prefix_on_amine"),
        pytest.param("OC(=O)CSSSCCO", "[(2-hydroxyethyl)trisulfanyl]acetic acid", id="homonuclear_run_of_three"),
        pytest.param("OC(=O)COOOC", "(methyltrioxidanyl)acetic acid", id="oxygen_run_of_three"),
        pytest.param("OC(=O)CSSOC", "[(methoxysulfanyl)sulfanyl]acetic acid", id="mixed_run_of_three"),
        pytest.param("OC(=O)CSSSCC(=O)O", "2,2'-trisulfanediyldiacetic acid", id="trisulfanediyl_linker"),
    ],
)
def test_chalcogen_chain_prefixes_beside_a_principal_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NCCSS", id="perthiol_outranks_amine"),
    ],
)
def test_chalcogen_chain_beside_a_principal_group_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC1CCCCC1COO", "(2-aminocyclohexyl)methaneperoxol", id="ring_amine"),
        pytest.param("NCC=CCOO", "4-aminobut-2-ene-1-peroxol", id="unsaturated_chain_hydroperoxide_amine"),
        pytest.param("N[C@@H](C)COO", "(2S)-2-aminopropane-1-peroxol", id="specified_stereocenter_hydroperoxide_amine"),
        pytest.param("CC(C)(OO)CCN(C)C", "4-(dimethylamino)-2-methylbutane-2-peroxol", id="tertiary_amine_prefix"),
        pytest.param("OOC1CCCc2ccccc12", "1,2,3,4-tetrahydronaphthalene-1-peroxol", id="fused_ring_peroxol"),
    ],
)
def test_peroxol_with_amino_prefixes_and_on_fused_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(=O)CC(=O)CC(N)=O", "3,5-dioxohexanamide", id="two_ketones__ketone_amide"),
        pytest.param("C=CC(=O)CC(N)=O", "3-oxopent-4-enamide", id="unsaturated_chain__ketone_amide"),
        pytest.param("NC(=O)C1CCC1=O", "2-oxocyclobutane-1-carboxamide", id="ring__ketone_amide"),
        pytest.param("OCC(=O)CC(N)=O", "4-hydroxy-3-oxobutanamide", id="hydroxyl_coexistence__ketone_amide"),
        pytest.param("CNC(=O)CC(=O)C", "N-methyl-3-oxobutanamide", id="n_substituted_amide__ketone_amide"),
    ],
)
def test_two_ketones__and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_ketone_amide():
    assert smiles_to_iupac("c1ccccc1CC(=O)CC(N)=O") == "3-oxo-4-phenylbutanamide"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C(N)=O") == "2-oxo-4-phenylbutanamide"


def test_phenyl_chain_ketone_amide_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(=O)CC(N)=O") == "4-(2-ethenylphenyl)-3-oxobutanamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCC(C)=O", "1-aminopropan-2-one"),
    ],
)
def test_smiles_to_iupac_ketone_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CNCC(C)=O", "1-(methylamino)propan-2-one", id="secondary_amine"),
        pytest.param("NC1CCC(=O)CC1", "4-aminocyclohexan-1-one", id="ring__ketone_amine"),
        pytest.param("NCC=CC(C)=O", "5-aminopent-3-en-2-one", id="unsaturated_chain__ketone_amine"),
        pytest.param("C[C@H](N)C(C)=O", "(3S)-3-aminobutan-2-one", id="specified_stereocenter__ketone_amine"),
        pytest.param("NCC(=O)CO", "1-amino-3-hydroxypropan-2-one", id="hydroxyl_coexisting__ketone_amine"),
        pytest.param("CC(=O)CC(=O)CC(=O)OC", "methyl 3,5-dioxohexanoate", id="two_ketones__ketone_ester"),
        pytest.param("C=CC(=O)CC(=O)OC", "methyl 3-oxopent-4-enoate", id="unsaturated_chain__ketone_ester"),
        pytest.param("COC(=O)C1CCC1=O", "methyl 2-oxocyclobutane-1-carboxylate", id="ring__ketone_ester"),
        pytest.param("OCC(=O)CC(=O)OC", "methyl 4-hydroxy-3-oxobutanoate", id="hydroxyl_coexistence__ketone_ester"),
    ],
)
def test_secondary_amine_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_chain_ketone_ester():
    assert smiles_to_iupac("c1ccccc1CC(=O)CC(=O)OC") == "methyl 3-oxo-4-phenylbutanoate"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C(=O)OC") == "methyl 2-oxo-4-phenylbutanoate"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=Cc1ccccc1CC(=O)CC(=O)OC", "methyl 4-(2-ethenylphenyl)-3-oxobutanoate", id="phenyl_chain_ketone_ester_unsaturation"),
        pytest.param("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34", "17-hydroxyestr-4-en-3-one", id="hydroxyl_alongside_steroid_ketone_is_named"),
    ],
)
def test_phenyl_chain_ketone_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC.CCC[N+](=O)[O-]",  # _nitro.py
    ],
)
def test_multi_fragment_rejected_instead_of_silently_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC(Cl)COCC(Cl)CO", "3,3'-oxybis(2-chloropropan-1-ol)"),
        ("OCCNc1ccccc1", "2-anilinoethan-1-ol"),
        ("CN(c1ccccc1)c1cccc(O)c1", "3-(N-methylanilino)phenol"),
        ("Oc1cccc(Nc2ccc(Cl)cc2)c1", "3-(4-chloroanilino)phenol"),
        ("Oc1cccc(N(CC)c2ccc(Cl)cc2)c1", "3-(4-chloro-N-ethylanilino)phenol"),
        ("Oc1cccc(N(c2ccc(Cl)cc2)c2ccccc2)c1", "3-(4-chloro-N-phenylanilino)phenol"),
        ("Oc1cccc(Nc2ccc3ccccc3c2)c1", "3-[(naphthalen-2-yl)amino]phenol"),
        ("OCC[N+](=O)[O-]", "2-nitroethan-1-ol"),
        ("OCCC#N", "3-hydroxypropanenitrile"),
    ],
)
def test_chain_parent_with_heteroatom_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COc1ccc(O)cc1", "4-methoxyphenol"),
        ("Oc1ccc(cc1)[N+](=O)[O-]", "4-nitrophenol"),
        ("Oc1ccccc1C(O)=O", "2-hydroxybenzoic acid"),
        ("COc1ccc(cc1)C=O", "4-methoxybenzaldehyde"),
    ],
)
def test_ring_parent_with_heteroatom_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COc1ccc(cc1)-c1ccccc1", "4-methoxy-1,1'-biphenyl"),
        ("OC(=O)c1ccc(cc1)-c1ccc(cc1)C(O)=O", "[1,1'-biphenyl]-4,4'-dicarboxylic acid"),
    ],
)
def test_identical_rings_joined_directly_form_a_ring_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC=C", "methoxyethene"),
        ("CSc1ccc(Cl)cc1", "1-chloro-4-(methylsulfanyl)benzene"),
        ("c1ccccc1CN(C)C1CCCCC1", "N-benzyl-N-methylcyclohexanamine"),
    ],
)
def test_parents_without_a_principal_group_and_n_substituted_amines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COC(=O)c1ccccc1C(O)=O", "2-(methoxycarbonyl)benzoic acid"),
    ],
)
def test_known_compounds_through_the_fallback_engines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(Nc1ccccc1)c1ccccc1", "N-phenylbenzamide"),
        ("CC(=O)N(C)CCO", "N-(2-hydroxyethyl)-N-methylacetamide"),
    ],
)
def test_n_substituted_amide_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereodescriptors_on_a_substituent_are_rejected_not_dropped_is_named():
    assert smiles_to_iupac("Oc1ccc(cc1)[C@H](C)Cl") == "4-[(1S)-1-chloroethyl]phenol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)CC(=CC(O)=O)C(O)=O", "prop-1-ene-1,2,3-tricarboxylic acid"),
    ],
)
def test_carbo_suffix_when_the_groups_do_not_fit_one_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC(Cl)c1ccc(C(Cl)CO)cc1", "2,2'-(1,4-phenylene)bis(2-chloroethan-1-ol)"),
    ],
)
def test_chain_units_on_a_ring_linker_are_named_multiplicatively(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC(CO)(CO)CO", "2,2-bis(hydroxymethyl)propane-1,3-diol"),
    ],
)
def test_chains_that_hold_every_group_stay_substitutive(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NC(=O)c1cccc2ccccc12", "naphthalene-1-carboxamide"),
        ("Oc1cc2ccccc2c2ccccc12", "phenanthren-9-ol"),
    ],
)
def test_fused_aromatic_parents_with_functional_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Nc1ccc(cc1)S(N)(=O)=O", "4-aminobenzene-1-sulfonamide"),
        ("OS(=O)(=O)c1ccc2ccccc2c1", "naphthalene-2-sulfonic acid"),
        ("NS(=O)(=O)c1ccc(Cl)cc1C(O)=O", "5-chloro-2-sulfamoylbenzoic acid"),
        ("CS(=O)(=O)CCO", "2-(methanesulfonyl)ethan-1-ol"),
        ("OC(=O)CS(O)(=O)=O", "sulfoacetic acid"),
        ("O=S(=O)(O)c1ccc(cc1)S(O)(=O)=O", "benzene-1,4-disulfonic acid"),
    ],
)
def test_sulfonic_acid_and_sulfonamide_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCN1CCN(CCO)CC1", "2,2'-(piperazine-1,4-diyl)di(ethan-1-ol)"),
        ("Cc1ccc(cc1)N1CCOCC1", "4-(4-methylphenyl)morpholine"),
    ],
)
def test_heterocyclic_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_thioketone_is_named():
    assert smiles_to_iupac("S=C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-thione"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)(C)OC(=O)NCC(=O)O", "[(tert-butoxycarbonyl)amino]acetic acid"),
        ("O=C(O)CNC(=O)c1ccccc1", "benzamidoacetic acid"),
    ],
)
def test_acyl_prefixes_with_substituents_and_alkoxycarbonylamino(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC[C@H]1CCCC[C@@H]1C", "[(1S,2S)-2-methylcyclohexyl]methanol"),
        ("OC(=O)CC[C@H]1CCCC[C@@H]1C", "3-[(1R,2S)-2-methylcyclohexyl]propanoic acid"),
        ("OC[C@H]1CCCC[C@@H]1Cl", "[(1R,2S)-2-chlorocyclohexyl]methanol"),
        ("OC(=O)CNC(=O)[C@H](Cl)CC", "{[(2R)-2-chlorobutanoyl]amino}acetic acid"),
        ("OC(=O)CNC(=O)[C@@H](N)c1ccccc1", "{[(2S)-amino(phenyl)acetyl]amino}acetic acid"),
        ("OC(=O)CCOC(=O)[C@@H]1CCC[C@H]1C", "3-{[(1R,2R)-2-methylcyclopentane-1-carbonyl]oxy}propanoic acid"),
    ],
)
def test_stereodescriptors_inside_substituent_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)COCCOCCOCCOCC(O)=O", "3,6,9,12-tetraoxatetradecane-1,14-dioic acid"),
        ("NCCNCCNCCNCCOCCCC(O)CCCC", "1-amino-12-oxa-3,6,9-triazaicosan-16-ol"),
        ("C[SiH2]C[SiH2]C[SiH2]C[SiH2]CCCl", "10-chloro-2,4,6,8-tetrasiladecane"),
        ("CSC[SiH2]CSC[SiH2]CCN", "2,6-dithia-4,8-disiladecan-10-amine"),
        ("C[SiH2]C[SiH2]C[SiH2]C[SiH2]CC(=O)C", "2,4,6,8-tetrasilaundecan-10-one"),
        ("[SiH3]C[SiH2]C[SiH2]C[SiH2]CC", "1,3,5,7-tetrasilanonane"),
        ("C[SiH2][PH]C[SiH2]C[SiH2]C", "3-phospha-2,5,7-trisilaoctane"),
    ],
)
def test_skeletal_replacement_parents_with_four_or_more_heterounits(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Oc1ccc(nc1)-c1ccccn1", "[2,2'-bipyridin]-5-ol"),
    ],
)
def test_heteroaromatic_assemblies_and_primed_locant_order(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCB(C)C", "2-(dimethylboranyl)ethan-1-ol"),
        ("C[Si](C)(C)OCCO", "2-[(trimethylsilyl)oxy]ethan-1-ol"),
        ("OCCNN", "2-hydrazinylethan-1-ol"),
    ],
)
def test_mononuclear_hydride_prefixes_and_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC1C2CC3CC(C2)CC1C3", "(adamantan-2-yl)methanol"),
        (
            "Cc1ccc(cc1)C(c1ccccc1)(c1ccccc1)c1ccccn1",
            "2-[(4-methylphenyl)di(phenyl)methyl]pyridine",
        ),
    ],
)
def test_retained_adamantane_and_one_carbon_prefix_multiplication(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCc1cnccn1", "(pyrazin-2-yl)methanol"),
    ],
)
def test_heteroaromatic_substituents_are_not_read_as_phenyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C1CCOc2ccccc12", "2,3-dihydro-4H-1-benzopyran-4-one"),
        ("O=C(O)c1ccc(=O)[nH]c1", "6-oxo-1,6-dihydropyridine-3-carboxylic acid"),
    ],
)
def test_added_hydrogen_ring_ketones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC1CCC2CCCC2C1", "5-methyloctahydro-1H-indene"),
        ("O=C1CCC2CCCCC2C1", "octahydronaphthalen-2(1H)-one"),
        ("C=C1CCC2CCCCC2C1", "2-methylidenedecahydronaphthalene"),
    ],
)
def test_saturated_fused_systems_use_hydro_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@H]1CCCN(C)C1", "(3S)-1,3-dimethylpiperidine"),
        ("C/C=C1\\C/C(=C/C)CC1", "(1Z,3E)-1,3-diethylidenecyclopentane"),
        ("C/C=C1\\CCCCC1C", "(1E)-1-ethylidene-2-methylcyclohexane"),
    ],
)
def test_ring_stereodescriptors_are_cited_once(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCCSCCCl", "1-chloro-2-[(2-chloroethyl)sulfanyl]ethane"),
    ],
)
def test_symmetric_ethers_and_sulfides_without_a_principal_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1CCC2C(C1)CCC1CCCCC21", "tetradecahydrophenanthren-2-ol"),
    ],
)
def test_saturated_tricyclic_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "OC(=O)[C@@H](N)CSSC[C@@H](N)C(=O)O",
            "(2R,2'S)-3,3'-disulfanediylbis(2-aminopropanoic acid)",
        ),
    ],
)
def test_multiplicative_names_with_dichalcogen_linkers_and_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@H](C(=O)O)[C@@H](C)C(N)=O", "(2S,3R)-4-amino-2,3-dimethyl-4-oxobutanoic acid"),
    ],
)
def test_ring_substituents_inside_chain_substituents_and_stereo_citation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCCC1[SiH2][SiH3]", "cyclohexyldisilane"),
        ("CNNNC", "1,3-dimethyltriazane"),
        ("C[SiH2]O[SiH2]O[SiH2]C", "1,5-dimethyltrisiloxane"),
        ("Cl[SiH2]O[SiH3]", "chlorodisiloxane"),
        ("c1ccccc1SOSC", "methyl(phenyl)dithioxane"),
        ("CSOSC", "dimethyldithioxane"),
        ("CSOS", "methyldithioxane"),
    ],
)
def test_heteroatom_chain_parent_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)[SiH3]", "silanecarboxylic acid"),
        ("OC(=O)[SiH2]O[SiH3]", "disiloxanecarboxylic acid"),
        ("OC(=O)[SiH2]O[SiH2]O[SiH3]", "trisiloxane-1-carboxylic acid"),
        ("C[SiH2]O[SiH2]C(=O)O", "3-methyldisiloxane-1-carboxylic acid"),
        ("OC(=O)[SiH2][SiH2]C(=O)O", "disilane-1,2-dicarboxylic acid"),
        ("NC(=O)[PH2]", "phosphanecarboxamide"),
        ("N#C[SiH3]", "silanecarbonitrile"),
        ("O=C[SiH3]", "silanecarbaldehyde"),
    ],
)
def test_carbo_suffix_on_heteroacyclic_parent_hydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)Cn1cccn1", "(1H-pyrazol-1-yl)acetic acid"),
    ],
)
def test_heteroaromatic_substituent_prefixes_cite_indicated_hydrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)CN1C=CC=CC1", "(pyridin-1(2H)-yl)acetic acid"),
        ("OC(=O)Cn1c2ccccc2c2ccccc21", "(9H-carbazol-9-yl)acetic acid"),
    ],
)
def test_added_hydrogen_substituent_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)CC1CCC2(CC1)OCCO2", "(1,4-dioxaspiro[4.5]decan-8-yl)acetic acid"),
    ],
)
def test_bridged_and_spiro_ring_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Cc1ccc2ccccc2c1-c1c(C)ccc2ccccc12", "2,2'-dimethyl-1,1'-binaphthalene"),
    ],
)
def test_assemblies_of_two_identical_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OS(=O)CCS(=O)(=O)O", "2-sulfinoethane-1-sulfonic acid"),
    ],
)
def test_sulfonic_acid_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NS(=O)(=O)CS(=O)(=O)O", "sulfamoylmethanesulfonic acid"),
    ],
)
def test_sulfonic_acid_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CNS(=O)(=O)CCS(=O)(=O)O", "2-(methylsulfamoyl)ethane-1-sulfonic acid", id="n_substituted_sulfonamide"),
        pytest.param("OCC(S(=O)(=O)N)S(=O)(=O)O", "2-hydroxy-1-sulfoethane-1-sulfonamide", id="other_heteroatom__sulfonic_acid_sulfonamide"),
        pytest.param("C=CC(S(=O)(=O)N)S(=O)(=O)O", "1-sulfoprop-2-ene-1-sulfonamide", id="unsaturated_chain__sulfonic_acid_sulfonamide"),
        pytest.param("C1CCC(S(=O)(=O)N)(CC1)S(=O)(=O)O", "1-sulfamoylcyclohexane-1-sulfonic acid", id="ring__sulfonic_acid_sulfonamide"),
    ],
)
def test_n_substituted_sulfonamide_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the -SO3H).
        ("SCS(=O)(=O)O", "sulfanylmethanesulfonic acid"),
    ],
)
def test_sulfonic_acid_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=CC(S)S(=O)(=O)O", "1-sulfanylprop-2-ene-1-sulfonic acid", id="unsaturated_chain__sulfonic_acid_thiol"),
        pytest.param("C1CCC(S)(CC1)S(=O)(=O)O", "1-sulfanylcyclohexane-1-sulfonic acid", id="ring__sulfonic_acid_thiol"),
    ],
)
def test_unsaturated_chain__and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1C(S)CCS(=O)(=O)O", "3-phenyl-3-sulfanylpropane-1-sulfonic acid"),  # CID 57312026
    ],
)
def test_phenyl_chain_sulfonic_acid_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Sc1ccccc1CS(=O)(=O)O", "(2-sulfanylphenyl)methanesulfonic acid", id="ring_with_second_substituent"),
        pytest.param("C=Cc1ccccc1CC(S)S(=O)(=O)O", "2-(2-ethenylphenyl)-1-sulfanylethane-1-sulfonic acid", id="chain_sulfonic_acid_thiol_unsaturation"),
    ],
)
def test_phenyl_ring_with_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_suffix_class_rank_has_no_duplicate_ranks_within_distinct_classes():
    assert len(SUFFIX_CLASS_RANK) == len(set(SUFFIX_CLASS_RANK.values()))


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCS", "2-aminoethane-1-thiol"),
    ],
)
def test_smiles_to_iupac_thiol_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NC1CCC(S)CC1", "4-aminocyclohexane-1-thiol", id="ring__thiol_amine"),
        pytest.param("NCC=CCS", "4-aminobut-2-ene-1-thiol", id="unsaturated_chain__thiol_amine"),
        pytest.param("C[C@H](N)CS", "(2S)-2-aminopropane-1-thiol", id="specified_stereocenter__thiol_amine"),
        pytest.param("NCC(S)CSCC", "1-amino-3-(ethylsulfanyl)propane-2-thiol", id="sulfide_coexisting"),
    ],
)
def test_ring__thiol_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("SC(=S)CCCCC(=S)S", "hexanebis(dithioic acid)"),
        ("CCCCCC(=S)S", "hexane(dithioic acid)"),
        ("OOC(=O)C(=O)OO", "ethanediperoxoic acid"),
        ("CC(=O)OS", "ethane(thioperoxoic) OS-acid"),
        ("OOC=O", "methaneperoxoic acid"),
        ("OC(=S)c1ccccc1", "benzenecarbothioic O-acid"),
        ("CS(=O)(=O)OO", "methanesulfonoperoxoic acid"),
        ("CS(=NO)(=O)O", "N-hydroxymethanesulfonimidic acid"),
        ("OC(=O)c1ccccc1S(=O)(=O)OC", "2-(methoxysulfonyl)benzoic acid"),
        ("CCCS(=O)(=O)c1ccccc1", "(propane-1-sulfonyl)benzene"),
    ],
)
def test_acid_functional_replacement_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)c1ccc(cc1)c1ccc2ccccc2c1", "4-(naphthalen-2-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)c1cccc2ccccc12", "4-(naphthalen-1-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)c1ccc2cccnc2c1", "4-(quinolin-7-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)c1ccc2[nH]ccc2c1", "4-(1H-indol-5-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)c1cccc2cc3ccccc3cc12", "4-(anthracen-1-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)C12CC3CC(CC(C3)C1)C2", "4-(adamantan-1-yl)benzoic acid"),
        ("OC(=O)C1CCC(CC1)C12CC3CC(CC(C3)C1)C2", "4-(adamantan-1-yl)cyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)C1C2CC3CC(C2)CC1C3", "4-(adamantan-2-yl)cyclohexane-1-carboxylic acid"),
        ("OC(=O)c1ccc(cc1)C1COCCOCCOCCOCCO1", "4-(1,4,7,10,13-pentaoxacyclopentadecan-2-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)C1COCCNCCOCCOCCN1", "4-(1,4,10-trioxa-7,13-diazacyclopentadecan-8-yl)benzoic acid"),
    ],
)
def test_fused_bridged_and_macrocyclic_substituents_on_a_carboxylic_acid_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)c1ccc(cc1)NOc1ccccc1", "4-(phenoxyamino)benzoic acid"),
        ("OC(=O)c1ccc(cc1)NOC", "4-(methoxyamino)benzoic acid"),
        ("OC(=O)c1ccc(cc1)N(C)OC", "4-[methoxy(methyl)amino]benzoic acid"),
        ("OC(=O)c1ccc(cc1)OSc1ccccc1", "4-[(phenylsulfanyl)oxy]benzoic acid"),
        ("OC(=O)c1ccc(cc1)SNC", "4-[(methylamino)sulfanyl]benzoic acid"),
        ("OC(=O)c1ccc(cc1)N=Nc1ccccc1", "4-(phenyldiazenyl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)N=N", "4-diazenylbenzoic acid"),
        ("OC(=O)c1ccc(cc1)NNc1ccccc1", "4-(2-phenylhydrazin-1-yl)benzoic acid"),
    ],
)
def test_heteroatom_to_heteroatom_substituent_groups_on_a_carboxylic_acid_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C1CCC(CC1)=N", "4-iminocyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)=NC", "4-(methylimino)cyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)=NO", "4-(hydroxyimino)cyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)=NOC", "4-(methoxyimino)cyclohexane-1-carboxylic acid"),
        ("OC(=O)CCC(=N)CCC(=O)O", "4-iminoheptanedioic acid"),
        ("OC(=O)c1ccc(cc1)C=N", "4-(methanimidoyl)benzoic acid"),
    ],
)
def test_imino_prefixes_under_a_senior_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "group, expected",
    [
        ("C(=N)C", "ethanimidoyl"),
        ("C(=N)", "methanimidoyl"),
        ("C(=N)CCC", "butanimidoyl"),
        ("C(=N)c1ccccc1", "benzenecarboximidoyl"),
        ("C(=N)C1CCCC1", "cyclopentanecarboximidoyl"),
        ("C(=NN)", "methanehydrazonoyl"),
        ("C(=NN)C", "ethanehydrazonoyl"),
        ("C(=NN)C1CCCCC1", "cyclohexanecarbohydrazonoyl"),
        ("C(=NO)c1ccccc1", "N-hydroxybenzenecarboximidoyl"),
        ("C(=S)CC", "propanethioyl"),
        ("C(=[Se])C", "ethaneselenoyl"),
        ("C(=S)C(=O)O", "carboxymethanethioyl"),
        ("C(=O)C(=S)O", "hydroxy(sulfanylidene)acetyl"),
        ("C(=S)C(=S)O", "hydroxy(sulfanylidene)ethanethioyl"),
        ("C(=S)C(=S)S", "sulfanyl(sulfanylidene)ethanethioyl"),
        ("C(=O)C(=O)Cl", "chloro(oxo)acetyl"),
        ("C(=O)C(N)=O", "oxamoyl"),
        ("C(=S)S", "dithiocarboxy"),
        ("C(=O)S", "sulfanylcarbonyl"),
        ("C(=S)O", "hydroxycarbonothioyl"),
        ("S(=O)(=O)c1ccccc1", "benzenesulfonyl"),
        ("S(=O)(=S)CC", "ethanesulfonothioyl"),
        ("S(=[Se])c1ccccc1", "benzenesulfinoselenoyl"),
        ("S(=N)CC", "ethanesulfinimidoyl"),
        ("[Se](=O)C", "methaneseleninyl"),
    ],
)
def test_acyl_prefixes_of_imidic_and_chalcogen_acids(group, expected):
    name = smiles_to_iupac("OC(=O)c1ccc(cc1)" + group)
    assert name.startswith("4-") and name.endswith("benzoic acid")
    assert expected in name


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=C=C1CCC(CC1)C(=O)O", "4-(oxomethylidene)cyclohexane-1-carboxylic acid"),
        ("S=C=C1CCC(CC1)C(=O)O", "4-(sulfanylidenemethylidene)cyclohexane-1-carboxylic acid"),
        ("[Se]=C=C1CCC(CC1)C(=O)O", "4-(selanylidenemethylidene)cyclohexane-1-carboxylic acid"),
        ("NN=C=C1CCC(CC1)C(=O)O", "4-(hydrazinylidenemethylidene)cyclohexane-1-carboxylic acid"),
    ],
)
def test_nonacyl_carbonic_acid_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)c1ccc(cc1)[SiH2]Cc1ccc(C)cc1", "4-{[(4-methylphenyl)methyl]silyl}benzoic acid"),
        ("OC(=O)c1ccc(cc1)C1CCC2(CC1)CCCP2", "4-(1-phosphaspiro[4.5]decan-8-yl)benzoic acid"),
        ("OC(=O)c1ccc(cc1)N1CC=NC2=NC=CN12", "4-(imidazo[1,2-b][1,2,4]triazin-1(2H)-yl)benzoic acid"),
        ("OC(=O)Cc1ccc(S(C)(=O)=O)cc1", "[4-(methanesulfonyl)phenyl]acetic acid"),
        ("OC(=O)Cc1cccc([SiH2]O[SiH3])n1", "(6-disiloxanylpyridin-2-yl)acetic acid"),
        (
            "OC(=O)C1CCC(CC1)=C1c2cccc(n2)Cc2cccc(n2)Cc2cccc(n2)Cc2cccc(n2)1",
            "4-[1,3,5,7(2,6)-tetrapyridinacyclooctaphan-2-ylidene]cyclohexane-1-carboxylic acid",
        ),
        ("OC(=O)c1ccc(cc1)[PH2]", "4-phosphanylbenzoic acid"),
        ("OC(=O)C1CCC(CC1)=[BH]", "4-boranylidenecyclohexane-1-carboxylic acid"),
    ],
)
def test_ring_and_hetero_groups_as_prefixes_use_one_ring_group_namer(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=C(Cl)c1ccc(S(=O)(=O)Cl)cc1", "4-(chlorosulfonyl)benzoyl chloride", id="acyl_halide_below_the_principal_acid_class_is_a_prefix"),
        pytest.param("CC(C)c1c(CO)nnn1-c1ccc(O)cc1", "4-[4-(hydroxymethyl)-5-(propan-2-yl)-1H-1,2,3-triazol-1-yl]phenol", id="aromatic_ring_nitrogen_cut_from_a_unit"),
        pytest.param(
            "CCNC(=O)c1cc(C(=O)NC)c2c(c1)[C@](CO)(c1ccccc1)CO2",
            "(3S)-N5-ethyl-3-(hydroxymethyl)-N7-methyl-3-phenyl-2,3-dihydro-1-benzofuran-5,7-dicarboxamide",
            id="n_locants_of_amides_on_a_fused_ring",
        ),
        pytest.param("CC(=O)OCOC(N)=O", "methylene acetate carbamate", id="carbamate_is_not_a_carboxylic_acyl_chain"),
        pytest.param(
            "Cc1cc(S(=O)(=O)N=[N+]=[N-])ccc1Oc1ccc(S(=O)(=O)N=[N+]=[N-])cc1",
            "4-[4-(azidosulfonyl)phenoxy]-3-methylbenzene-1-sulfonyl azide",
            id="sulfonyl_azide_is_not_a_sulfamoyl_group",
        ),
    ],
)
def test_polyfunctional_never_misattributes_a_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C=CCC(CN)CC(C)/C(=C\\C)CC/C=C\\C",
    ],
)
def test_group_or_stereo_element_without_a_supported_citation_raises(smiles):
    with pytest.raises(NotImplementedError):
        smiles_to_iupac(smiles)


def test_acid_carbon_bonded_to_a_chain_heteroatom_is_the_chain_end():
    assert smiles_to_iupac("OC(=O)[SiH2]C[SiH2]C[SiH2]C[SiH2]C") == "2,4,6,8-tetrasilanonan-1-oic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc2c(c1)Cc1ccc(cc12)-c1ccc(cc1)C(=O)O", "6-(4-carboxyphenyl)-9H-fluorene-2-carboxylic acid"),
        ("OC(=O)c1ccc(cc1)-c1ccncc1C(=O)O", "4-(4-carboxyphenyl)pyridine-3-carboxylic acid"),
        ("OC(=O)c1ccc(cc1)C1CCC(CC1)C(=O)O", "4-(4-carboxycyclohexyl)benzoic acid"),
        ("OC(=O)CCc1c(CC(=O)O)c(CC(=O)O)cc2ccccc12", "3-[2,3-bis(carboxymethyl)naphthalen-1-yl]propanoic acid"),
    ],
)
def test_senior_ring_is_the_parent_and_a_multiplicative_centre_keeps_all_principal_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(c1ccco1)C(O)c1ccco1", "1,2-di(furan-2-yl)-2-hydroxyethan-1-one"),
        ("n1c(-c2cccs2)cccc1-c1cccs1", "2,6-di(thiophen-2-yl)pyridine"),
        ("OC(c1ccc(C)cc1)c1ccc(C)cc1", "bis(4-methylphenyl)methanol"),
    ],
)
def test_simple_ring_groups_take_di_and_substituted_ones_take_bis(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("PPPPP", "pentaphosphane"),
        ("CO[SiH2]CC[SiH2]SC", "2-oxa-7-thia-3,6-disilaoctane"),
        ("[SiH3]OCS[SiH3]", "2-oxa-4-thia-1,5-disilapentane"),
        ("CB(C)COCCOCB(C)C", "2,9-dimethyl-4,7-dioxa-2,9-diboradecane"),
        (
            "C(C)OP(OCC)OCCOCC[N+](C)(C)C",
            "4-ethoxy-N,N,N-trimethyl-3,5,8-trioxa-4-phosphadecan-10-aminium",
        ),
        (
            "O=C(CC(=O)OC)OCCOC(CC(OCCOC(CC(=O)OC)=O)=O)=O",
            "dimethyl 3,8,10,15-tetraoxo-4,7,11,14-tetraoxaheptadecane-1,17-dioate",
        ),
        (
            "C(C)(=O)NC(CCCNC(C)=O)CC(NCCCC(CC(NCCCC(CC(NCCCC(CC(=O)OC)NC(C)=O)=O)NC(C)=O)=O)NC(C)=O)=O",
            "methyl 7,14,21,28-tetraacetamido-2,9,16,23-tetraoxo-3,10,17,24-tetraazatriacontan-30-oate",
        ),
        ("CC(=O)N(C)CCC(O)=O", "3-(N-methylacetamido)propanoic acid"),
        ("NC(=O)CCOCCOCCOCCOCCNC(=O)C", "1-acetamido-3,6,9,12-tetraoxapentadecan-15-amide"),
    ],
)
def test_skeletal_replacement_chains_and_amido_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("COCSSCCOCC[Se]C", "2,8-dioxa-4,5-dithia-11-selenadodecane", id="adjacent_sulfur_atoms_in_a_replacement_chain"),
        pytest.param("COCCSSCCOC", "1-methoxy-2-[(2-methoxyethyl)disulfanyl]ethane", id="two_ether_chains_and_a_disulfide_are_three_units"),
        pytest.param("CCOCCOCCSSCCOCCOCC", "3,6,13,16-tetraoxa-9,10-dithiaoctadecane", id="long_chain_with_a_disulfide"),
    ],
)
def test_skeletal_replacement_chains_with_adjacent_chalcogen_atoms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CNCNCNCNC", "2,4,6,8-tetraazanonane"),
        ("CCNC(=O)CCOCCOCCOCCOC", "N-ethyl-2,5,8,11-tetraoxatetradecan-14-amide"),
        ("CON(C)C(=O)ON=C(C#N)C(N)=O", "7-cyano-3-methyl-4-oxo-2,5-dioxa-3,6-diazaoct-6-en-8-amide"),
        ("NC(=N)NCCCCCNC(=N)NCCCCCNC(N)=N", "9-imino-2,8,10,16-tetraazaheptadecane-1,17-diimidamide"),
        (
            "O=C(NCCNCCNCCNCCNCCN)NCCNCCNCCNCCNCCN",
            "16-amino-N-(14-amino-3,6,9,12-tetraazatetradecan-1-yl)-2,5,8,11,14-pentaazahexadecan-1-amide",
        ),
        (
            "FC(=O)NSNCON=CC=NOCNSNC(F)=O",
            "6,11-dioxa-3,14-dithia-2,4,7,10,13,15-hexaazahexadeca-7,9-diene-1,16-dioyl difluoride",
        ),
        ("C[SiH2]C[SiH2]C[SiH2]C[SiH2]C=O", "2,4,6,8-tetrasilanonan-1-al"),
        ("N#C[SiH2]C[SiH2]C[SiH2]C[SiH2]C#N", "2,4,6,8-tetrasilanonane-1,9-dinitrile"),
        ("CNC(=O)[SiH2]C[SiH2]C[SiH2]C[SiH2]C", "N-methyl-2,4,6,8-tetrasilanonan-1-amide"),
    ],
)
def test_skeletal_replacement_parents_with_groups_on_heteroatom_bonded_carbons_and_nitrogen_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC1CCCC(C)N1NS(=O)(=O)c1ccc(NN)cc1", id="ring_nitrogen_is_not_a_hydrazine_chain_atom"),
        pytest.param("CN(C)S(=O)(=O)c1ccc(CNC(=NC)N2CCSCC2)cc1", id="ring_nitrogen_of_a_carbamimidoyl_group"),
    ],
)
def test_ring_nitrogens_are_never_chain_or_acyl_group_atoms_with_ring_copies(smiles):
    with pytest.raises(NotImplementedError):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(O)COC(=O)NCCN", "2-hydroxypropyl (2-aminoethyl)carbamate", id="blue_book_example"),
        pytest.param(
            "N#Cc1c(F)cc(C#CCCNC(=O)OCc2ccccc2)cc1F",
            "benzyl [4-(4-cyano-3,5-difluorophenyl)but-3-yn-1-yl]carbamate",
            id="nitrile_and_halogens_are_prefixes",
        ),
        pytest.param(
            "CC(C)(C)OC(=O)N(CCCNC(=O)c1cccnc1SC(F)F)Cc1cccnc1",
            "tert-butyl (3-{2-[(difluoromethyl)sulfanyl]pyridine-3-carboxamido}propyl)[(pyridin-3-yl)methyl]carbamate",
            id="two_nitrogen_substituents_and_a_junior_amide",
        ),
        pytest.param(
            "C[C@@H](NC(=O)OCc1ccccc1)C(=O)NCC",
            "benzyl [(2R)-1-(ethylamino)-1-oxopropan-2-yl]carbamate",
            id="stereocentre_in_the_nitrogen_substituent",
        ),
    ],
)
def test_carbamate_esters_are_the_parent_of_polyfunctional_molecules(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "Cc1ncoc1C(=O)Nc1cc(C#CCO)ccc1Cl",
            "N-[2-chloro-5-(3-hydroxyprop-1-yn-1-yl)phenyl]-4-methyl-1,3-oxazole-5-carboxamide",
            id="alkynyl_with_hydroxy_on_a_ring_substituent",
        ),
        pytest.param(
            "FC(F)(F)c1ccc(C=Cn2ccc3ccccc32)cc1",
            "1-{2-[4-(trifluoromethyl)phenyl]ethen-1-yl}-1H-indole",
            id="ethenyl_between_two_rings",
        ),
        pytest.param(
            "O=c1[nH]c2ccccc2nc1C=C(O)c1cccc(Br)c1",
            "3-[2-(3-bromophenyl)-2-hydroxyethen-1-yl]quinoxalin-2(1H)-one",
            id="ethenyl_with_hydroxy_and_aryl",
        ),
        pytest.param(
            "Cc1c(C)c(O)c(CCC(C)(C)O)c(C=Cc2cnco2)c1O",
            "2-(3-hydroxy-3-methylbutyl)-5,6-dimethyl-3-[2-(1,3-oxazol-5-yl)ethen-1-yl]benzene-1,4-diol",
            id="heteroaryl_ethenyl_on_a_diol",
        ),
    ],
)
def test_unsaturated_substituents_carry_heteroatoms_and_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NNc1ccccn1", "2-hydrazinylpyridine", id="blue_book_example"),
        pytest.param("NNC1=NCCN1", "2-hydrazinyl-4,5-dihydro-1H-imidazole", id="blue_book_dihydroimidazole"),
        pytest.param(
            "Cc1ncsc1C(Cc1cc(F)ccc1F)NN",
            "5-[2-(2,5-difluorophenyl)-1-hydrazinylethyl]-4-methyl-1,3-thiazole",
            id="hydrazinyl_on_a_chain_of_a_thiazole",
        ),
        pytest.param(
            "CC(C)CC1(C(Cc2ccn(C)n2)NN)CCCC1",
            "3-{2-hydrazinyl-2-[1-(2-methylpropyl)cyclopentyl]ethyl}-1-methyl-1H-pyrazole",
            id="hydrazinyl_before_methyl_in_alphanumerical_order",
        ),
        pytest.param(
            "CC(NNC)Cc1ccn(C)n1",
            "1-methyl-3-[2-(2-methylhydrazin-1-yl)propyl]-1H-pyrazole",
            id="substituted_hydrazinyl_cites_its_free_valence",
        ),
    ],
)
def test_hydrazine_beside_a_ring_with_nitrogen_is_a_hydrazinyl_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)NNC(=O)c1ccccc1C(=O)O", "2-(2-acetylhydrazine-1-carbonyl)benzoic acid", id="acyl_on_the_far_nitrogen"),
        pytest.param("CN(N)C(=O)c1ccccc1C(=O)O", "2-(1-methylhydrazine-1-carbonyl)benzoic acid", id="substituent_on_the_near_nitrogen"),
        pytest.param("CN(C)NC(=O)c1ccccc1C(=O)O", "2-(2,2-dimethylhydrazine-1-carbonyl)benzoic acid", id="two_substituents_on_the_far_nitrogen"),
        pytest.param("CNNC(=S)c1ccccc1C(=O)O", "2-(2-methylhydrazine-1-carbothioyl)benzoic acid", id="thio_analogue"),
        pytest.param(
            "COc1ccc(C(=O)NNC(=O)Cc2ccc(F)cc2)cc1S(=O)(=O)Nc1ccc(C)cc1",
            "5-{2-[(4-fluorophenyl)acetyl]hydrazine-1-carbonyl}-2-methoxy-N-(4-methylphenyl)benzene-1-sulfonamide",
            id="sulfonamide_outranks_a_diacylhydrazine",
        ),
        pytest.param(
            "Cc1ccc(NS(=O)(=O)c2c[nH]c(C(=O)NN)c2)cc1F",
            "N-(3-fluoro-4-methylphenyl)-5-(hydrazinecarbonyl)-1H-pyrrole-3-sulfonamide",
            id="sulfonamide_outranks_a_hydrazide",
        ),
    ],
)
def test_hydrazide_groups_beside_an_amide_class_parent_are_hydrazinecarbonyl_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
