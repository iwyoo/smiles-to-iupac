import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._coexisting_groups import name_via_senior_acyclic
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._seniority import SUFFIX_CLASS_RANK


def test_unsaturated_chain():
    assert smiles_to_iupac("NCC=CCO") == "4-aminobut-2-en-1-ol"


def test_specified_stereocenter():
    assert smiles_to_iupac("C[C@H](N)CO") == "(2S)-2-aminopropan-1-ol"


def test_ether_coexisting():
    assert smiles_to_iupac("NCC(O)COCC") == "1-amino-3-ethoxypropan-2-ol"


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
        ("NCC=O", "2-aminoethanal"),
    ],
)
def test_smiles_to_iupac_aldehyde_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring():
    assert smiles_to_iupac("NC1CCC(C=O)CC1") == "4-aminocyclohexane-1-carbaldehyde"


def test_unsaturated_chain__aldehyde_amine():
    assert smiles_to_iupac("NCC=CC=O") == "4-aminobut-2-enal"


def test_specified_stereocenter__aldehyde_amine():
    assert smiles_to_iupac("C[C@H](N)C=O") == "(2S)-2-aminopropanal"


def test_hydroxyl_coexisting():
    assert smiles_to_iupac("NC(CO)C=O") == "2-amino-3-hydroxypropanal"


def test_branch_substituent():
    assert smiles_to_iupac("O=CC(C)C(=O)O") == "2-methyl-3-oxopropanoic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(=O)C(Cl)C=O") == "2-chloro-3-oxobutanal"


def test_two_aldehydes_with_ketone():
    assert smiles_to_iupac("O=CC(=O)CC=O") == "2-oxobutanedial"


def test_unsaturated_chain__aldehyde_ketone():
    assert smiles_to_iupac("C=CC(=O)CC=O") == "3-oxopent-4-enal"


def test_ring__aldehyde_ketone():
    assert smiles_to_iupac("O=CC1CCC(=O)C1") == "3-oxocyclopentane-1-carbaldehyde"


def test_hydroxyl_coexistence():
    assert smiles_to_iupac("OCC(=O)CC=O") == "4-hydroxy-3-oxobutanal"


def test_phenyl_chain_aldehyde_ketone():
    assert smiles_to_iupac("c1ccccc1CC(=O)CC=O") == "3-oxo-4-phenylbutanal"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C=O") == "2-oxo-4-phenylbutanal"


def test_phenyl_substituted_benzene_ring_aldehyde_ketone():
    assert smiles_to_iupac("Cc1ccccc1CC(=O)CC=O") == "4-(2-methylphenyl)-3-oxobutanal"


def test_phenyl_chain_aldehyde_ketone_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(=O)CC=O") == "4-(2-ethenylphenyl)-3-oxobutanal"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCC(N)=O", "3-aminopropanamide"),
    ],
)
def test_smiles_to_iupac_amide_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_amines():
    assert smiles_to_iupac("NCC(N)C(N)=O") == "2,3-diaminopropanamide"


def test_ring__amide_amine():
    assert smiles_to_iupac("NC1CCCCC1C(N)=O") == "2-aminocyclohexane-1-carboxamide"


def test_unsaturated_chain__amide_amine():
    assert smiles_to_iupac("NCC=CC(N)=O") == "4-aminobut-2-enamide"


def test_specified_stereocenter__amide_amine():
    assert smiles_to_iupac("N[C@@H](C)C(N)=O") == "(2S)-2-aminopropanamide"


def test_hydroxyl_coexisting__amide_amine():
    assert smiles_to_iupac("NCC(O)C(N)=O") == "3-amino-2-hydroxypropanamide"


def test_trimethylhydrazinium_ide_amine_imide():
    # Blue Book P-74.2.1.3 worked example: 1,2,2-trimethylhydrazin-2-ium-1-ide.
    assert smiles_to_iupac("C[N-][NH+](C)C") == "1,2,2-trimethylhydrazin-2-ium-1-ide"


def test_n_substituted_amide():
    assert smiles_to_iupac("OC(=O)CC(=O)NC") == "2-(methylcarbamoyl)ethanoic acid"


def test_amide_on_ring():
    assert smiles_to_iupac("OC(=O)C1CCCCC1C(=O)N") == "2-carbamoylcyclohexane-1-carboxylic acid"


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


def test_phenyl_directly_attached_to_amine_carbon():
    assert smiles_to_iupac("c1ccccc1C(N)C(=O)O") == "2-amino-2-phenylethanoic acid"


def test_phenyl_substituted_benzene_ring_carboxylic_acid_amine():
    assert smiles_to_iupac("Cc1ccccc1CC(N)C(=O)O") == "2-amino-3-(2-methylphenyl)propanoic acid"


def test_phenyl_chain_carboxylic_acid_amine_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(N)C(=O)O") == "2-amino-3-(2-ethenylphenyl)propanoic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)C[Se](=O)O", "2-seleninoethanoic acid"),
    ],
)
def test_carboxylic_acid_seleninic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_seleninic_acid_still_routes_normally():
    assert smiles_to_iupac("C[Se](=O)O") == "methaneseleninic acid"


def test_multiple_carboxylic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC([Se](=O)O)C(=O)O")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C(O)C[Se](=O)O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C=CC[Se](=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC([Se](=O)O)CC1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)S(=O)O", "sulfinomethanoic acid"),
    ],
)
def test_carboxylic_acid_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_carboxylic_acids_not_supported__carboxylic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)CC(S(=O)O)C(=O)O")


def test_other_heteroatom_not_supported__carboxylic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C(O)CS(=O)O")


def test_unsaturated_chain_not_supported__carboxylic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C=CCS(=O)O")


def test_ring_not_supported__carboxylic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC(S(=O)O)CC1")


def test_phenyl_chain_carboxylic_acid_sulfinic_acid_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(S(=O)O)C(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)S(=O)(=O)N", "sulfamoylmethanoic acid"),
    ],
)
def test_carboxylic_acid_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_carboxylic_acids():
    assert smiles_to_iupac("OC(=O)CC(S(=O)(=O)N)C(=O)O") == "2-sulfamoylbutanedioic acid"


def test_other_heteroatom():
    assert smiles_to_iupac("OC(=O)C(O)CS(=O)(=O)N") == "2-hydroxy-3-sulfamoylpropanoic acid"


def test_unsaturated_chain__carboxylic_acid_sulfonamide():
    assert smiles_to_iupac("OC(=O)C=CCS(=O)(=O)N") == "4-sulfamoylbut-2-enoic acid"


def test_ring__carboxylic_acid_sulfonamide():
    assert (
        smiles_to_iupac("OC(=O)C1CCC(S(=O)(=O)N)CC1") == "4-sulfamoylcyclohexane-1-carboxylic acid"
    )


def test_multiple_carboxylic_acids__carboxylic_acid_sulfonic_acid():
    assert smiles_to_iupac("OC(=O)CC(S(=O)(=O)O)C(=O)O") == "2-sulfobutanedioic acid"


def test_other_heteroatom__carboxylic_acid_sulfonic_acid():
    assert smiles_to_iupac("OC(=O)C(O)CS(=O)(=O)O") == "2-hydroxy-3-sulfopropanoic acid"


def test_unsaturated_chain__carboxylic_acid_sulfonic_acid():
    assert smiles_to_iupac("OC(=O)C=CCS(=O)(=O)O") == "4-sulfobut-2-enoic acid"


def test_ring__carboxylic_acid_sulfonic_acid():
    assert smiles_to_iupac("OC(=O)C1CCC(S(=O)(=O)O)CC1") == "4-sulfocyclohexane-1-carboxylic acid"


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
            "2,2-bis[(ethanoyloxy)methyl]propane-1,3-diyl diethanoate",
        ),
        ("CC(=O)OCC=CCOC(C)=O", "but-2-ene-1,4-diyl diethanoate"),
    ],
)
def test_diester_acyloxy_generalized_backbones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC(=O)Oc1ccc(Cl)cc1OC(=O)c1ccc(Cl)cc1",
            "4-chloro-1,2-phenylene 2-(4-chlorobenzoate) 1-ethanoate",
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
        ("ClCC(=O)OCCCOC(=O)CCl", "propane-1,3-diyl bis(2-chloroethanoate)"),
    ],
)
def test_acyclic_diyl_with_generalized_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1ccc(N(=O)=O)cc1OC(C)=O", "4-nitro-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(N(C)CC)cc1OC(C)=O", "4-[ethyl(methyl)amino]-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(SC)cc1OC(C)=O", "4-(methylsulfanyl)-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(Oc2ccccc2)cc1OC(C)=O", "4-phenoxy-1,2-phenylene diethanoate"),
    ],
)
def test_functional_substituents_on_diyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(OCCOC(=O)c1ccncc1)c1ccncc1", "ethane-1,2-diyl di(pyridine-4-carboxylate)"),
        ("CC(=O)OCCOC(=O)CC#N", "ethane-1,2-diyl 2-cyanoethanoate ethanoate"),
        ("CC(=O)OC1CCCCC1OC(=O)CC(C)=O", "cyclohexane-1,2-diyl ethanoate 3-oxobutanoate"),
        ("CC(=O)Oc1ccccc1OC(=O)/C=C/C", "1,2-phenylene (2E)-but-2-enoate ethanoate"),
        (
            "CC(=O)OC1CCCC1OC(=O)[C@H](C)Cl",
            "cyclopentane-1,2-diyl (2S)-2-chloropropanoate ethanoate",
        ),
    ],
)
def test_generalized_acyl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OC1CCC(COC(C)=O)CC1", "4-[(ethanoyloxy)methyl]cyclohexyl ethanoate"),
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


def test_ring_with_extra_substituent():
    assert smiles_to_iupac("CC1CCCCC1CC1CCCCC1") == "1-(cyclohexylmethyl)-2-methylcyclohexane"


def test_amine_on_alcohol_part():
    assert smiles_to_iupac("CC(=O)OCCN") == "2-aminoethyl ethanoate"


def test_two_amines__ester_amine():
    assert smiles_to_iupac("NCC(N)C(=O)OC") == "methyl 2,3-diaminopropanoate"


def test_two_esters_names_polyester():
    assert smiles_to_iupac("NCC(=O)OCOC(=O)C") == "methylene 2-aminoethanoate ethanoate"


def test_ring__ester_amine():
    assert smiles_to_iupac("NC1CCCCC1C(=O)OC") == "methyl 2-aminocyclohexane-1-carboxylate"


def test_unsaturated_chain__ester_amine():
    assert smiles_to_iupac("NCC=CC(=O)OC") == "methyl 4-aminobut-2-enoate"


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](C)C(=O)OC")


def test_hydroxyl_coexisting__ester_amine():
    assert smiles_to_iupac("NCC(O)C(=O)OC") == "methyl 3-amino-2-hydroxypropanoate"


def test_branched_alkoxy_r_prime():
    assert smiles_to_iupac("CC(C)OCC=O") == "2-(propan-2-yloxy)ethanal"


def test_halogen_on_main_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C=O") == "2-chloro-3-methoxypropanal"


def test_specified_stereocenter__ether_aldehyde():
    assert smiles_to_iupac("COC[C@@H](C)C=O") == "(2R)-3-methoxy-2-methylpropanal"


def test_branched_alkoxy_r_prime__ether_amide():
    assert smiles_to_iupac("CC(C)OCC(N)=O") == "2-(propan-2-yloxy)ethanamide"


def test_halogen_on_main_chain_still_works__ether_amide():
    assert smiles_to_iupac("COCC(Cl)C(N)=O") == "2-chloro-3-methoxypropanamide"


def test_two_amides():
    assert smiles_to_iupac("NC(=O)C(COC)C(N)=O") == "2-(methoxymethyl)propanediamide"


def test_two_ethers():
    assert smiles_to_iupac("COCC(OC)C(N)=O") == "2,3-dimethoxypropanamide"


def test_ring__ether_amide():
    assert smiles_to_iupac("NC(=O)C1CCCCC1COC") == "2-(methoxymethyl)cyclohexane-1-carboxamide"


def test_specified_stereocenter__ether_amide():
    assert smiles_to_iupac("COC[C@@H](C)C(N)=O") == "(2R)-3-methoxy-2-methylpropanamide"


def test_branched_alkoxy_r_prime__ether_amine():
    assert smiles_to_iupac("CC(C)OCCN") == "2-(propan-2-yloxy)ethanamine"


def test_halogen_on_main_chain_still_works__ether_amine():
    assert smiles_to_iupac("COCC(Cl)CN") == "2-chloro-3-methoxypropan-1-amine"


def test_two_amines__ether_amine():
    assert smiles_to_iupac("NCC(N)COC") == "3-methoxypropane-1,2-diamine"


def test_unsaturated_chain__ether_amine():
    assert smiles_to_iupac("NCC=CCOC") == "4-methoxybut-2-en-1-amine"


def test_specified_stereocenter__ether_amine():
    assert smiles_to_iupac("N[C@@H](C)COC") == "(2S)-1-methoxypropan-2-amine"


def test_branched_alkoxy_r_prime__ether_ester():
    assert smiles_to_iupac("CC(C)OCC(=O)OC") == "methyl 2-(propan-2-yloxy)ethanoate"


def test_halogen_on_acyl_chain_still_works():
    assert smiles_to_iupac("COCC(Cl)C(=O)OC") == "methyl 2-chloro-3-methoxypropanoate"


def test_ether_on_alcohol_part():
    assert smiles_to_iupac("CC(=O)OCCOC") == "2-methoxyethyl ethanoate"


def test_ring__ether_ester():
    assert (
        smiles_to_iupac("O=C(OC)C1CCCCC1COC") == "methyl 2-(methoxymethyl)cyclohexane-1-carboxylate"
    )


def test_specified_stereocenter__ether_ester():
    assert smiles_to_iupac("COC[C@@H](C)C(=O)OC") == "methyl (2R)-3-methoxy-2-methylpropanoate"


def test_branched_alkoxy_r_prime__ether_hydroperoxide():
    assert smiles_to_iupac("CC(C)OCCOO") == "2-(propan-2-yloxy)ethaneperoxol"


def test_halogen_on_main_chain_still_works__ether_hydroperoxide():
    assert smiles_to_iupac("COCC(Cl)COO") == "2-chloro-3-methoxypropane-1-peroxol"


def test_two_ethers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COCC(OC)COO")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOCC=CCOC")


def test_specified_stereocenter_raises__ether_hydroperoxide():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OO[C@@H](C)COC")


def test_branched_alkoxy_r_prime__ether_ketone():
    assert smiles_to_iupac("CC(C)OCC(C)=O") == "1-(propan-2-yloxy)propan-2-one"


def test_two_ketones():
    assert smiles_to_iupac("O=CC(=O)COC") == "3-methoxy-2-oxopropanal"


def test_ring__ether_ketone():
    assert smiles_to_iupac("O=C1CCCCC1COC") == "2-(methoxymethyl)cyclohexan-1-one"


def test_specified_stereocenter__ether_ketone():
    assert smiles_to_iupac("COC[C@@H](C)C(C)=O") == "(3R)-4-methoxy-3-methylbutan-2-one"


def test_branched_alkoxy_r_prime__ether_thiol():
    assert smiles_to_iupac("CC(C)OCCS") == "2-(propan-2-yloxy)ethanethiol"


def test_halogen_on_main_chain_still_works__ether_thiol():
    assert smiles_to_iupac("COCC(Cl)CS") == "2-chloro-3-methoxypropane-1-thiol"


def test_unsaturated_chain__ether_thiol():
    assert smiles_to_iupac("SCC=CCOC") == "4-methoxybut-2-ene-1-thiol"


def test_specified_stereocenter__ether_thiol():
    assert smiles_to_iupac("S[C@@H](C)COC") == "(2S)-1-methoxypropane-2-thiol"


def test_chloroacetylpiperidine():
    assert smiles_to_iupac("ClCC(=O)N1CCCCC1") == "2-chloro-1-(piperidin-1-yl)ethan-1-one"


def test_branched_acyl_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)N1CCCCC1")


def test_substituted_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N1CCC(C)CC1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCOO", "2-aminoethaneperoxol"),
    ],
)
def test_smiles_to_iupac_hydroperoxide_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCCC1COO")


def test_unsaturated_chain_raises__hydroperoxide_amine():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC=CCOO")


def test_specified_stereocenter_raises__hydroperoxide_amine():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](C)COO")


def test_two_ketones__ketone_amide():
    assert smiles_to_iupac("CC(=O)CC(=O)CC(N)=O") == "3,5-dioxohexanamide"


def test_unsaturated_chain__ketone_amide():
    assert smiles_to_iupac("C=CC(=O)CC(N)=O") == "3-oxopent-4-enamide"


def test_ring__ketone_amide():
    assert smiles_to_iupac("NC(=O)C1CCC1=O") == "2-oxocyclobutane-1-carboxamide"


def test_hydroxyl_coexistence__ketone_amide():
    assert smiles_to_iupac("OCC(=O)CC(N)=O") == "4-hydroxy-3-oxobutanamide"


def test_n_substituted_amide__ketone_amide():
    assert smiles_to_iupac("CNC(=O)CC(=O)C") == "N-methyl-3-oxobutanamide"


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


def test_secondary_amine():
    assert smiles_to_iupac("CNCC(C)=O") == "1-(methylamino)propan-2-one"


def test_ring__ketone_amine():
    assert smiles_to_iupac("NC1CCC(=O)CC1") == "4-aminocyclohexan-1-one"


def test_unsaturated_chain__ketone_amine():
    assert smiles_to_iupac("NCC=CC(C)=O") == "5-aminopent-3-en-2-one"


def test_specified_stereocenter__ketone_amine():
    assert smiles_to_iupac("C[C@H](N)C(C)=O") == "(3S)-3-aminobutan-2-one"


def test_hydroxyl_coexisting__ketone_amine():
    assert smiles_to_iupac("NCC(=O)CO") == "1-amino-3-hydroxypropan-2-one"


def test_two_ketones__ketone_ester():
    assert smiles_to_iupac("CC(=O)CC(=O)CC(=O)OC") == "methyl 3,5-dioxohexanoate"


def test_unsaturated_chain__ketone_ester():
    assert smiles_to_iupac("C=CC(=O)CC(=O)OC") == "methyl 3-oxopent-4-enoate"


def test_ring__ketone_ester():
    assert smiles_to_iupac("COC(=O)C1CCC1=O") == "methyl 2-oxocyclobutane-1-carboxylate"


def test_hydroxyl_coexistence__ketone_ester():
    assert smiles_to_iupac("OCC(=O)CC(=O)OC") == "methyl 4-hydroxy-3-oxobutanoate"


def test_phenyl_chain_ketone_ester():
    assert smiles_to_iupac("c1ccccc1CC(=O)CC(=O)OC") == "methyl 3-oxo-4-phenylbutanoate"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C(=O)OC") == "methyl 2-oxo-4-phenylbutanoate"


def test_phenyl_chain_ketone_ester_unsaturation():
    assert (
        smiles_to_iupac("C=Cc1ccccc1CC(=O)CC(=O)OC") == "methyl 4-(2-ethenylphenyl)-3-oxobutanoate"
    )


def test_hydroxyl_alongside_steroid_ketone_is_named():
    assert smiles_to_iupac("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34") == "17-hydroxyestr-4-en-3-one"


@pytest.mark.parametrize(
    "smiles",
    [
        "CC.CC1CCC12CNC2",  # _spiro_heteroatom.py
        "CC.CCOCC",  # _ether.py
        "CC.CCC[N+](=O)[O-]",  # _nitro.py
        "CC.CCS(C)(=O)=O",  # _sulfone.py
        "CC.Cc1ccncc1",  # _hetero_monocyclic.py substituent path
    ],
)
def test_multi_fragment_rejected_instead_of_silently_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC(Cl)COCC(Cl)CO", "3,3'-oxybis(2-chloropropan-1-ol)"),
        ("OCCNc1ccccc1", "2-anilinoethanol"),
        ("OCC[N+](=O)[O-]", "2-nitroethanol"),
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
        ("COC(=O)c1ccccc1C(O)=O", "methyl hydrogen benzene-1,2-dicarboxylate"),
    ],
)
def test_known_compounds_through_the_fallback_engines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(Nc1ccccc1)c1ccccc1", "N-phenylbenzamide"),
        ("CC(=O)N(C)CCO", "N-(2-hydroxyethyl)-N-methylethanamide"),
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
        ("CS(=O)(=O)CCO", "2-(methanesulfonyl)ethanol"),
        ("OC(=O)CS(O)(=O)=O", "2-sulfoethanoic acid"),
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


@pytest.mark.parametrize(
    "smiles",
    [
        "S=C1CCCC=C1C",
    ],
)
def test_thioketones_are_rejected(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)(C)OC(=O)NCC(=O)O", "2-[(tert-butoxycarbonyl)amino]ethanoic acid"),
        ("O=C(O)CNC(=O)c1ccccc1", "2-(benzoylamino)ethanoic acid"),
    ],
)
def test_acyl_prefixes_with_substituents_and_alkoxycarbonylamino(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC[C@H]1CCCC[C@@H]1C", "[(1S,2S)-2-methylcyclohexyl]methanol"),
        ("OC(=O)CC[C@H]1CCCC[C@@H]1C", "3-[(1R,2S)-2-methylcyclohexyl]propanoic acid"),
    ],
)
def test_stereodescriptors_inside_substituent_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)COCCOCCOCCOCC(O)=O", "3,6,9,12-tetraoxatetradecanedioic acid"),
        ("NCCNCCNCCNCCOCCCC(O)CCCC", "1-amino-12-oxa-3,6,9-triazaicosan-16-ol"),
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
        ("OCCB(C)C", "2-(dimethylboranyl)ethanol"),
        ("C[Si](C)(C)OCCO", "2-[(trimethylsilyl)oxy]ethanol"),
        ("OCCNN", "2-hydrazinylethanol"),
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
            "OC(=O)[C@@H](N)CSSC[C@H](N)C(=O)O",
            "(2R,2'R)-3,3'-disulfanediylbis(2-aminopropanoic acid)",
        ),
    ],
)
def test_multiplicative_names_with_dichalcogen_linkers_and_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@H](C(=O)O)[C@@H](C)C(N)=O", "(2S,3R)-3-carbamoyl-2-methylbutanoic acid"),
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
    ],
)
def test_heteroatom_chain_parent_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)Cn1cccn1", "2-(1H-pyrazol-1-yl)ethanoic acid"),
    ],
)
def test_heteroaromatic_substituent_prefixes_cite_indicated_hydrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)CN1C=CC=CC1", "2-[pyridin-1(2H)-yl]ethanoic acid"),
        ("OC(=O)Cn1c2ccccc2c2ccccc21", "2-(9H-carbazol-9-yl)ethanoic acid"),
    ],
)
def test_added_hydrogen_substituent_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)CC1CCC2(CC1)OCCO2", "2-(1,4-dioxaspiro[4.5]decan-8-yl)ethanoic acid"),
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
        ("OS(=O)CCS(=O)(=O)O", "2-sulfinoethanesulfonic acid"),
    ],
)
def test_sulfonic_acid_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_other_heteroatom_not_supported__sulfonic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(S(=O)O)S(=O)(=O)O")


def test_unsaturated_chain_not_supported__sulfonic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(S(=O)O)S(=O)(=O)O")


def test_ring_not_supported__sulfonic_acid_sulfinic_acid():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(S(=O)O)(CC1)S(=O)(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NS(=O)(=O)CS(=O)(=O)O", "sulfamoylmethanesulfonic acid"),
    ],
)
def test_sulfonic_acid_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_sulfonamide():
    assert smiles_to_iupac("CNS(=O)(=O)CCS(=O)(=O)O") == "2-(methylsulfamoyl)ethanesulfonic acid"


def test_other_heteroatom__sulfonic_acid_sulfonamide():
    assert smiles_to_iupac("OCC(S(=O)(=O)N)S(=O)(=O)O") == "2-hydroxy-1-sulfoethanesulfonamide"


def test_unsaturated_chain__sulfonic_acid_sulfonamide():
    assert smiles_to_iupac("C=CC(S(=O)(=O)N)S(=O)(=O)O") == "1-sulfoprop-2-ene-1-sulfonamide"


def test_ring__sulfonic_acid_sulfonamide():
    assert (
        smiles_to_iupac("C1CCC(S(=O)(=O)N)(CC1)S(=O)(=O)O")
        == "1-sulfamoylcyclohexane-1-sulfonic acid"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the -SO3H).
        ("SCS(=O)(=O)O", "sulfanylmethanesulfonic acid"),
    ],
)
def test_sulfonic_acid_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_chain__sulfonic_acid_thiol():
    assert smiles_to_iupac("C=CC(S)S(=O)(=O)O") == "1-sulfanylprop-2-ene-1-sulfonic acid"


def test_ring__sulfonic_acid_thiol():
    assert smiles_to_iupac("C1CCC(S)(CC1)S(=O)(=O)O") == "1-sulfanylcyclohexane-1-sulfonic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1C(S)CCS(=O)(=O)O", "3-phenyl-3-sulfanylpropane-1-sulfonic acid"),  # CID 57312026
    ],
)
def test_phenyl_chain_sulfonic_acid_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_ring_with_second_substituent():
    assert smiles_to_iupac("Sc1ccccc1CS(=O)(=O)O") == "(2-sulfanylphenyl)methanesulfonic acid"


def test_phenyl_chain_sulfonic_acid_thiol_unsaturation():
    assert (
        smiles_to_iupac("C=Cc1ccccc1CC(S)S(=O)(=O)O")
        == "2-(2-ethenylphenyl)-1-sulfanylethanesulfonic acid"
    )


def test_suffix_class_rank_has_no_duplicate_ranks_within_distinct_classes():
    assert len(SUFFIX_CLASS_RANK) == len(set(SUFFIX_CLASS_RANK.values()))


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCS", "2-aminoethanethiol"),
    ],
)
def test_smiles_to_iupac_thiol_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring__thiol_amine():
    assert smiles_to_iupac("NC1CCC(S)CC1") == "4-aminocyclohexane-1-thiol"


def test_unsaturated_chain__thiol_amine():
    assert smiles_to_iupac("NCC=CCS") == "4-aminobut-2-ene-1-thiol"


def test_specified_stereocenter__thiol_amine():
    assert smiles_to_iupac("C[C@H](N)CS") == "(2S)-2-aminopropane-1-thiol"


def test_sulfide_coexisting():
    assert smiles_to_iupac("NCC(S)CSCC") == "1-amino-3-(ethylsulfanyl)propane-2-thiol"


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
        ("OC(=O)c1ccc(cc1)NNc1ccccc1", "4-(2-phenylhydrazinyl)benzoic acid"),
    ],
)
def test_heteroatom_to_heteroatom_substituent_groups_on_a_carboxylic_acid_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles", ["CNO", "OOC1CCCCC1", "COOCOC"])
def test_heteroatom_connections_are_not_prefixes_without_a_senior_group(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C1CCC(CC1)=N", "4-iminocyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)=NC", "4-(methylimino)cyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)=NO", "4-(hydroxyimino)cyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(CC1)=NOC", "4-(methoxyimino)cyclohexane-1-carboxylic acid"),
        ("OC(=O)CCC(=N)CCC(=O)O", "4-iminoheptanedioic acid"),
        ("OC(=O)c1ccc(cc1)C=N", "4-(iminomethyl)benzoic acid"),
    ],
)
def test_imino_prefixes_under_a_senior_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
