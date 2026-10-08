import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)CC[Si](C)(C)[Si](C)(C)CCC(O)=O",
            "3,3'-(1,1,2,2-tetramethyldisilane-1,2-diyl)dipropanoic acid",
        ),
    ],
)
def test_heteroatom_and_concatenated_linkers_between_chain_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCCOc1ccc(cc1)Cc1ccc(OCCO)cc1", "2,2'-[methylenebis(4,1-phenyleneoxy)]di(ethan-1-ol)"),
        (
            "OCCOCCOc1ccc(OCCOCCO)cc1",
            "2,2'-[1,4-phenylenebis(oxyethane-2,1-diyloxy)]di(ethan-1-ol)",
        ),
    ],
)
def test_ring_and_assembly_components_inside_a_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H](O)CN(C)C[C@H](C)O", "(2S,2'S)-1,1'-(methylazanediyl)di(propan-2-ol)"),
        ("S(/C1=C/CCCCCC1)/C1=C/CCCCCC/1", "(1E,1'E)-1,1'-sulfanediyldi(cyclooct-1-ene)"),
        ("C1=C/[C@@H](CC[C@@H]2/C=C/CCCCC2)CCCCC/1", "(1E,1'E,3R,3'R)-3,3'-(ethane-1,2-diyl)di(cyclooct-1-ene)"),
        ("[C@H]1(CC[C@H]3CCCC34CCCC4)CCCC12CCCC2", "(1R,1'R)-1,1'-(ethane-1,2-diyl)di(spiro[4.4]nonane)"),
        ("CC[C@@H](C)c1ccc(Sc2ccc([C@H](C)CC)cc2)cc1", "1,1'-sulfanediylbis{4-[(2R)-butan-2-yl]benzene}"),
    ],
)
def test_stereodescriptors_of_multiplied_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "C1=C\\[C@H](CC[C@@H]2/C=C/CCCCC2)CCCCC/1",
            "(1Z,3S)-3-{2-[(1R,2E)-cyclooct-2-en-1-yl]ethyl}cyclooct-1-ene",
        ),
        (
            "C1=C\\[C@@H](CC[C@H]2/C=C/CCCCC2)CCCCC/1",
            "(1Z,3R)-3-{2-[(1S,2E)-cyclooct-2-en-1-yl]ethyl}cyclooct-1-ene",
        ),
        ("[C@@H]1(CC[C@H]3CCCC34CCCC4)CCCC12CCCC2", "(1R)-1-{2-[(1S)-spiro[4.4]nonan-1-yl]ethyl}spiro[4.4]nonane"),
        ("[C@H]1(CC[C@@H]3CCCC34CCCC4)CCCC12CCCC2", "(1R)-1-{2-[(1S)-spiro[4.4]nonan-1-yl]ethyl}spiro[4.4]nonane"),
    ],
)
def test_units_of_different_configuration_take_the_senior_one_as_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(=O)NCCOCCNC(C)=O", "N,N'-[oxydi(ethane-2,1-diyl)]diacetamide"),
    ],
)
def test_amide_units_joined_through_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCN(CC)OBCCBON(CC)CC", "N,N'-[ethane-1,2-diylbis(boranediyloxy)]bis(N-ethylethanamine)"),
    ],
)
def test_amine_units_joined_through_nitrogen_to_a_heteroatom_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=C(O)CNC(=O)C(=O)NCC(=O)O", "2,2'-[oxalylbis(azanediyl)]diacetic acid"),
        ("O=C(O)CNC(=O)CCC(=O)NCC(=O)O", "2,2'-[(1,4-dioxobutane-1,4-diyl)bis(azanediyl)]diacetic acid"),
        ("OC(=O)c1ccccc1C(=O)CCC(=S)c1ccccc1C(=O)O", "2,2'-(1-oxo-4-sulfanylidenebutane-1,4-diyl)dibenzoic acid"),
        ("O=C(N)C(=O)NCC(=O)O", "(oxamoylamino)acetic acid"),
    ],
)
def test_diacyl_linking_groups_and_oxamoylamino(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("[SiH3][SiH2]C[SiH2][SiH3]", "1,1'-methylenebis(disilane)"),
        ("[SiH3]CC[SiH2]CC[SiH3]", "[silanediyldi(ethane-2,1-diyl)]bis(silane)"),
        ("C[Si](C)(C)C[Si](C)(C)C", "methylenebis(trimethylsilane)"),
        ("c1ccccc1P(c1ccccc1)CCP(c1ccccc1)c1ccccc1", "(ethane-1,2-diyl)bis(diphenylphosphane)"),
        (
            "c1ccccc1P(c1ccccc1)CP(c1ccccc1)CP(c1ccccc1)c1ccccc1",
            "[(phenylphosphanediyl)bis(methylene)]bis(diphenylphosphane)",
        ),
        ("C[Si](C)(C)C[Si](C)(C)C", "methylenebis(trimethylsilane)"),
    ],
)
def test_heteroatom_hydride_units_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_linker_with_a_substituent_that_carries_the_principal_group():
    assert (
        smiles_to_iupac("OC(=O)CN(CC(O)=O)CCN(CC(=O)O)CCN(CC(O)=O)CC(O)=O")
        == "2,2',2'',2'''-{[(carboxymethyl)azanediyl]bis(ethane-2,1-diylnitrilo)}tetraacetic acid"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CNCCNCCNCCNC", "2,5,8,11-tetraazadodecane"),
        ("NCCNCCNCCN", "N1,N2-bis(2-aminoethyl)ethane-1,2-diamine"),
        ("NCCN(CCN(C)C)C", "N1-(2-aminoethyl)-N1,N2,N2-trimethylethane-1,2-diamine"),
        ("NCCNCCNCCNCN", "N1-{2-[(2-aminoethyl)amino]ethyl}-N2-(aminomethyl)ethane-1,2-diamine"),
        (
            "CC(CN(CC(=C)C)CC(=C)C)(C)C",
            "N-(2,2-dimethylpropyl)-2-methyl-N-(2-methylprop-2-en-1-yl)prop-2-en-1-amine",
        ),
        ("C1CNCCNCCN1", "1,4,7-triazonane"),
    ],
)
def test_hetero_chain_and_macrocycle_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[O]1CC[O](C)->[Pt]<-1([Cl])[Cl]", "dichlorido(1,2-dimethoxyethane-κ2O,O')platinum"),
        (
            "C1C[O]2->[K]3456<-[O]1CC[O]->3CC[O]->4CC[O]->5CC[O]->6CC2",
            "(1,4,7,10,13,16-hexaoxacyclooctadecane-κ6O,O',O'',O''',O'''',O''''')potassium",
        ),
    ],
)
def test_polydentate_ligand_complexes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C=C1CCC(=CC(=O)O)CC1", "2,2'-(cyclohexane-1,4-diylidene)diacetic acid"),
        (
            "OC(=O)C=C1C(Cl)=C2C=CC=CC2=CC1=CC(=O)O",
            "2,2'-(1-chloronaphthalene-2,3-diylidene)diacetic acid",
        ),
    ],
)
def test_linker_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CC(c1ccccc1)c1ccccc1", "diphenylacetaldehyde"),
        ("N#CC(c1ccccc1)c1ccccc1", "diphenylacetonitrile"),
        ("SC(c1ccccc1)c1ccccc1", "diphenylmethanethiol"),
        ("NC(=O)C(c1ccccc1)c1ccccc1", "2,2-diphenylacetamide"),
        ("OS(=O)(=O)C(c1ccccc1)c1ccccc1", "diphenylmethanesulfonic acid"),
    ],
)
def test_chain_parent_cites_each_aromatic_ring_as_a_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=Cc1ccc(Cc2ccccc2)cc1", "4-benzylbenzaldehyde"),
        ("N#Cc1ccc(Cc2ccccc2)cc1", "4-benzylbenzonitrile"),
        ("Sc1ccc(Cc2ccccc2)cc1", "4-benzylbenzenethiol"),
        ("NC(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzamide"),
        ("OS(=O)(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzene-1-sulfonic acid"),
        ("NS(=O)(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzenesulfonamide"),
        ("Nc1ccc(Cc2ccccc2)cc1", "4-benzylaniline"),
        ("Oc1ccc(Cc2ccccc2)cc1", "4-benzylphenol"),
    ],
)
def test_ring_parent_cites_other_aromatic_rings_as_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_functional_group_on_a_heteroaromatic_ring_beside_another_ring_is_rejected_is_named():
    assert smiles_to_iupac("OC(=O)c1ccncc1Cc1ccccc1") == "3-benzylpyridine-4-carboxylic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1S(=O)(=O)c1ccccc1", "1,1'-sulfonyldibenzene"),
        ("c1ccccc1C#Cc1ccccc1", "1,1'-(ethyne-1,2-diyl)dibenzene"),
        ("c1ccc2ccccc2c1Cc1cccc2ccccc12", "1,1'-methylenedinaphthalene"),
        ("BrC1=CC(=CC=C1)CC1=CC(=CC=C1)Cl", "1-bromo-3-[(3-chlorophenyl)methyl]benzene"),
    ],
)
def test_unsubstituted_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCc1ccc(Cc2ccc(CCl)c(Br)c2)cc1Br", "1,1'-methylenebis[3-bromo-4-(chloromethyl)benzene]"),
    ],
)
def test_substituted_units_and_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(Oc2ccc(C(=O)O)c(Br)c2)cc1Br", "4,4'-oxybis(2-bromobenzoic acid)"),
        ("Nc1ccc(Cc2ccc(N)cc2)cc1", "4,4'-methylenedianiline"),
        (
            "OC(=O)c1ccc(Cc2ccc(C(=O)O)c(C(=O)O)c2)cc1C(=O)O",
            "4,4'-methylenedi(benzene-1,2-dicarboxylic acid)",
        ),
        ("OC(=O)c1ccc(C(=O)c2ccc(C(=O)O)cc2)cc1", "4,4'-carbonyldibenzoic acid"),
        ("OC(=O)c1ccc(C(O)c2ccc(C(=O)O)cc2)cc1", "4,4'-(hydroxymethylene)dibenzoic acid"),
        ("O=C1CCC(CC1)CC1CCC(=O)CC1", "4,4'-methylenedi(cyclohexan-1-one)"),
    ],
)
def test_functionalized_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(N(c2ccc(C(=O)O)cc2)c2ccc(C(=O)O)cc2)cc1", "4,4',4''-nitrilotribenzoic acid"),
    ],
)
def test_three_and_more_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC=CCc1ccccc1", "1,1'-(but-2-ene-1,4-diyl)dibenzene"),
        ("c1ccccc1C(C1CCCCC1)c1ccccc1", "1,1'-(cyclohexylmethylene)dibenzene"),
        (
            "c1ccccc1CC(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1",
            "1,1'-(2,2-dibenzylpropane-1,3-diyl)dibenzene",
        ),
        (
            "c1ccccc1C(SCc1ccccc1)SCc1ccccc1",
            "1,1'-[(phenylmethylene)bis(sulfanediylmethylene)]dibenzene",
        ),
    ],
)
def test_central_group_forms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Oc1ccc(OCC(C)COc2ccc(O)cc2)cc1", "4,4'-[(2-methylpropane-1,3-diyl)bis(oxy)]diphenol"),
        (
            "Oc1ccc(CCCCCCCCCCCCCCOCCCCCCCCCCCCCCc2ccc(O)cc2)cc1",
            "4,4'-[oxydi(tetradecane-14,1-diyl)]diphenol",
        ),
        (
            "C(COc1ccoc1)Oc1cocc1OCCOc1ccoc1",
            "3,3'-[furan-3,4-diylbis(oxyethane-2,1-diyloxy)]difuran",
        ),
    ],
)
def test_concatenated_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1/C=C/c1ccccc1", "1,1'-[(1E)-ethene-1,2-diyl]dibenzene"),
        ("C[C@H](c1ccccc1)[C@@H](C)c1ccccc1", "1,1'-[(2R,3S)-butane-2,3-diyl]dibenzene"),
    ],
)
def test_linker_stereodescriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccc(Cc2ccc(Cc3ccc(Cc4ccccc4)nc3)cc2)cc1", "3(2,5)-pyridina-1,7(1),5(1,4)-tribenzenaheptaphane"),
        pytest.param("c1ccc(C2CCC(C3CCC(c4ccc(C5CCC(C6CCC(C7CCCCC7)CC6)CC5)cc4)CC3)CC2)cc1",
            "1(1),4(1,4)-dibenzena-2,3,5,6(1,4),7(1)-pentacyclohexanaheptaphane", marks=pytest.mark.slow),
    ],
)
def test_linear_phane_with_any_ring_amplificants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_linear_phane_is_named_as_a_phane():
    assert (
        smiles_to_iupac("c1ccccc1Oc1cccc(Oc2cccc(Oc3ccccc3)c2)c1")
        == "2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane"
    )


def test_substituted_fused_unit():
    assert smiles_to_iupac("Clc1ccc2ccccc2c1Cc1c(Cl)ccc2ccccc12") == "1,1'-methylenebis(2-chloronaphthalene)"
    assert (
        smiles_to_iupac("C(=O)(NC1=CC=C2C=CC(=CC2=C1)S(=O)(=O)O)NC1=CC=C2C=CC(=CC2=C1)S(=O)(=O)O")
        == "7,7'-[carbonylbis(azanediyl)]di(naphthalene-2-sulfonic acid)"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccsc1Cc1ccco1", "2-[(thiophen-2-yl)methyl]furan"),
    ],
)
def test_unequal_rings_are_ranked_by_unsaturation_heteroatom_and_substituent_count(
    smiles, expected
):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCc1ccc(cc1)Oc1ccc(CCO)cc1", "2,2'-[oxybis(4,1-phenylene)]di(ethan-1-ol)"),
        ("OCCc1cc2ccccc2cc1CCO", "2,2'-(naphthalene-2,3-diyl)di(ethan-1-ol)"),
    ],
)
def test_chain_units_on_fused_and_composite_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


_ARYL_ACID = "OC(=O)c1ccc(cc1)"
_ACID_ARYL = "c1ccc(cc1)C(=O)O"


@pytest.mark.parametrize(
    "linker, expected",
    [
        ("C1CCNCC1", "piperidine-3,4-diyl"),
        ("C1CC(NCC1)", "piperidine-2,4-diyl"),
        ("C1CCOC1", "oxolane-2,3-diyl"),
        ("C1COCC1", "oxolane-3,4-diyl"),
        ("c1ccc2cc(ccc2c1)", "naphthalene-2,6-diyl"),
        ("c1cc2cc(ccc2[nH]1)", "1H-indole-2,5-diyl"),
        ("c1ccc2nc(ccc2c1)", "quinoline-2,6-diyl"),
        ("C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2,5-diyl"),
        ("C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2,7-diyl"),
        ("C1CCC2(CC1)CCCC2", "spiro[4.5]decane-1,8-diyl"),
        ("C1CC(C)CCC1", "4-methylcyclohexane-1,2-diyl"),
    ],
)
def test_ring_system_diyl_components_of_multiplicative_names(linker, expected):
    assert smiles_to_iupac(_ARYL_ACID + linker + _ACID_ARYL) == f"4,4'-({expected})dibenzoic acid"


@pytest.mark.parametrize(
    "linker, expected",
    [
        ("NN", "hydrazine-1,2-diyl"),
        ("N(C)N", "1-methylhydrazine-1,2-diyl"),
        ("N(C)N(C)", "1,2-dimethylhydrazine-1,2-diyl"),
        ("NNN", "triazane-1,3-diyl"),
        ("[PH][PH]", "diphosphane-1,2-diyl"),
        ("[PH][PH][PH]", "triphosphane-1,3-diyl"),
    ],
)
def test_homonuclear_heteroatom_hydride_linkers(linker, expected):
    assert smiles_to_iupac(_ARYL_ACID + linker + _ACID_ARYL) == f"4,4'-({expected})dibenzoic acid"


def test_diazenediyl_linker_cites_no_locants():
    assert smiles_to_iupac(_ARYL_ACID + "N=N" + _ACID_ARYL) == "4,4'-diazenediyldibenzoic acid"


def test_hydrazine_and_diazene_without_a_principal_group_are_the_parent():
    assert smiles_to_iupac("c1ccccc1NNc1ccccc1") == "1,2-diphenylhydrazine"
    assert smiles_to_iupac("c1ccccc1N=Nc1ccccc1") == "diphenyldiazene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)CC1Cc2ccc(cc2)C(CC(=O)O)c2ccc1cc2",
            "2,2'-(1,3(1,4)-dibenzenacyclopentaphane-2,4-diyl)diacetic acid",
        ),
        (
            "OC(=O)CC1Cc2ccc(cc2)CC(CC(=O)O)c2ccc1cc2",
            "2,2'-(1,4(1,4)-dibenzenacyclohexaphane-2,6-diyl)diacetic acid",
        ),
    ],
)
def test_phane_diyl_between_chain_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "linker, expected",
    [
        ("[AsH][AsH]", "diarsane-1,2-diyl"),
        ("[SbH][SbH]", "distibane-1,2-diyl"),
        ("[SnH2][SnH2]", "distannane-1,2-diyl"),
    ],
)
def test_heavier_group_14_and_15_hydride_linkers(linker, expected):
    assert smiles_to_iupac(_ARYL_ACID + linker + _ACID_ARYL) == f"4,4'-({expected})dibenzoic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)CC[AsH][AsH]CCC(=O)O", "3,3'-(diarsane-1,2-diyl)dipropanoic acid"),
        ("OC(=O)CC[PH][PH]CCC(=O)O", "3,3'-(diphosphane-1,2-diyl)dipropanoic acid"),
        ("OC(=O)CCNNCCC(=O)O", "3,3'-(hydrazine-1,2-diyl)dipropanoic acid"),
        ("OC(=O)CCN=NCCC(=O)O", "3,3'-diazenediyldipropanoic acid"),
    ],
)
def test_heteroatom_hydride_linkers_between_chain_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)c1ccc(cc1)C[SiH2]CC[SiH2]c1ccc(cc1)C(=O)O",
            "4-[(2-{[(4-carboxyphenyl)methyl]silyl}ethyl)silyl]benzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)[SiH2]C(F)C[SiH2]C(F)Cc1ccc(cc1)C(=O)O",
            "4-[(2-{[2-(4-carboxyphenyl)-1-fluoroethyl]silyl}-1-fluoroethyl)silyl]benzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)OCC(Cl)c1ccc(C)cc1C(=O)O",
            "2-[2-(4-carboxyphenoxy)-1-chloroethyl]-5-methylbenzoic acid",
        ),
    ],
)
def test_unsymmetrical_linkers_and_units_are_named_substitutively(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccccc1C1CPCC(C1)c1ccccc1", "3,5-diphenylphosphinane"),
        ("c1ccccc1C1CC[BH]C1c1ccccc1", "2,3-diphenylborolane"),
    ],
)
def test_ring_seniority_ranks_every_heteroatom_of_the_element_order(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1[Sn](C)(C)C1CCCCC1O", "2,2'-(dimethylstannanediyl)di(cyclohexan-1-ol)"),
        ("OC1CCCCC1C(O)C1CCCCC1O", "2,2'-(hydroxymethylene)di(cyclohexan-1-ol)"),
        ("OC1CCCCC1C(O)C1CCCCC1", "2-[cyclohexyl(hydroxy)methyl]cyclohexan-1-ol"),
    ],
)
def test_substituted_metal_linker_and_principal_group_count_between_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_ring_unit_with_three_different_attachments_is_not_cited_as_a_triyl():
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("CC(c1cc(CCCCS)c(CCCCS)cc1)C(C)c1cc(CCCCS)c(CCCCS)cc1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CN(C)C(=O)CCSCCC(=O)N(C)C", "3,3'-sulfanediylbis(N,N-dimethylpropanamide)"),
        ("C(N1CCCCCCCCCCC1)N1CCCCCCCCCCC1", "1,1'-methylenebis(1-azacyclododecane)"),
        ("C(N1CCOCC1)N1CCOCC1", "4,4'-methylenedimorpholine"),
        ("C(CN=CCC(=O)O)N=CCC(=O)O", "3,3'-[ethane-1,2-diylbis(azanylylidene)]dipropanoic acid"),
        (
            "OC(=O)c1ccc(CC(C(=O)OC)c2ccc(C(O)=O)cc2)cc1",
            "4,4'-(3-methoxy-3-oxopropane-1,2-diyl)dibenzoic acid",
        ),
    ],
)
def test_identical_parents_with_substituted_units_hetero_rings_imine_and_ester_linkers_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCCS(=O)CCCN", "3,3'-sulfinyldi(propan-1-amine)"),
        ("N(S(=O)CCCN)S(=O)CCCN", "3,3'-[azanediylbis(sulfinyl)]di(propan-1-amine)"),
    ],
)
def test_sulfinyl_linkers_of_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC[C@@H](C)c1ccc(Sc2ccc([C@H](C)CC)cc2)cc1",
            "1,1'-sulfanediylbis{4-[(2R)-butan-2-yl]benzene}",
        ),
        (
            "CC[C@@H](C)c1ccc(Sc2ccc([C@@H](C)CC)cc2)cc1",
            "1-[(2R)-butan-2-yl]-4-({4-[(2S)-butan-2-yl]phenyl}sulfanyl)benzene",
        ),
        (
            "C[C@@H](Br)C1([C@H](C)Br)CCCC1",
            "1-[(1R)-1-bromoethyl]-1-[(1S)-1-bromoethyl]cyclopentane",
        ),
    ],
)
def test_identical_configuration_units_are_multiplied_and_r_is_cited_before_s(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("n1nc(CO)c(Oc2c(CO)nncc2CO)c(CO)c1", "[oxydi(pyridazine-4,3,5-triyl)]tetramethanol"),
        ("n1nc(CO)c(CO)c(Oc2cnnc(CO)c2CO)c1", "[oxydi(pyridazine-5,3,4-triyl)]tetramethanol"),
        ("[SiH3]CCCCCCCCCCCCCCOCCCCCCCCCCCCCC[SiH3]", "[oxydi(tetradecane-14,1-diyl)]bis(silane)"),
        ("[SiH3]c1cc([SiH3])cc([SiH3])c1", "(benzene-1,3,5-triyl)tris(silane)"),
        ("CC(=O)N(C)[Si](C)(C=C)N(C)C(C)=O", "N,N'-[ethenyl(methyl)silanediyl]bis(N-methylacetamide)"),
    ],
)
def test_multiplicative_linker_locants_enclosure_and_multiplying_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Cl[SiH2]CC[SiH3]", "chloro(2-silylethyl)silane"),
        ("O=C(N=Nc1ccccc1)N=Nc1ccccc1", "bis(phenyldiazenyl)methanone"),
        ("S=C(N=N)N=N", "bis(diazenyl)methanethione"),
        ("[SiH3][SiH2]OCCS[SiH2][SiH3]", "3-oxa-6-thia-1,2,7,8-tetrasilaoctane"),
        ("[SiH3][SiH2]COO[SiH2][SiH3]", "[(disilanylmethyl)peroxy]disilane"),
        (
            "OC1CCC(CC1)=CCN=CCC1CCC(O)CC1",
            "4-(2-{[2-(4-hydroxycyclohexyl)ethylidene]amino}ethylidene)cyclohexan-1-ol",
        ),
        (
            "OC(=O)c1ccccc1CCOCCOCCOCCOCCc1ccccc1C(O)=O",
            "2,2'-(3,6,9,12-tetraoxatetradecane-1,14-diyl)dibenzoic acid",
        ),
    ],
)
def test_parent_choice_beside_identical_hydrides_acyl_diazenes_and_heterounits(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
