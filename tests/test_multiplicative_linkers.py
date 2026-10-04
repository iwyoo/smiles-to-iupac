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
    ],
)
def test_stereodescriptors_of_multiplied_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(=O)NCCOCCNC(C)=O", "N,N'-[oxydi(ethane-2,1-diyl)]diethanamide"),
    ],
)
def test_amide_units_joined_through_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("[SiH3][SiH2]C[SiH2][SiH3]", "1,1'-methylenebis(disilane)"),
        ("[SiH3]CC[SiH2]CC[SiH3]", "[silanediyldi(ethane-2,1-diyl)]bis(silane)"),
        ("C[Si](C)(C)C[Si](C)(C)C", "methylenebis(trimethylsilane)"),
        ("c1ccccc1P(c1ccccc1)CCP(c1ccccc1)c1ccccc1", "ethane-1,2-diylbis(diphenylphosphane)"),
        (
            "c1ccccc1P(c1ccccc1)CP(c1ccccc1)CP(c1ccccc1)c1ccccc1",
            "[(phenylphosphanediyl)bis(methylene)]bis(diphenylphosphane)",
        ),
        ("C[Si](C)(C)C[Si](C)(C)C", "methylenebis(trimethylsilane)"),
    ],
)
def test_heteroatom_hydride_units_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_linker_with_a_substituent_that_carries_the_principal_group():
    assert (
        smiles_to_iupac("OC(=O)CN(CC(O)=O)CCN(CC(=O)O)CCN(CC(O)=O)CC(O)=O")
        == "2,2',2'',2'''-{[(carboxymethyl)azanediyl]bis(ethane-2,1-diylnitrilo)}tetraacetic acid"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CNCCNCCNCCNC", "2,5,8,11-tetraazadodecane"),
        ("NCCNCCNCCN", "N1-{2-[(2-aminoethyl)amino]ethyl}ethane-1,2-diamine"),
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
        ("O=CC(c1ccccc1)c1ccccc1", "2,2-diphenylethanal"),
        ("N#CC(c1ccccc1)c1ccccc1", "2,2-diphenylethanenitrile"),
        ("SC(c1ccccc1)c1ccccc1", "diphenylmethanethiol"),
        ("NC(=O)C(c1ccccc1)c1ccccc1", "2,2-diphenylethanamide"),
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1/C=C/c1ccccc1", "1,1'-[(1E)-ethene-1,2-diyl]dibenzene"),
        ("C[C@H](c1ccccc1)[C@@H](C)c1ccccc1", "1,1'-[(2R,3S)-butane-2,3-diyl]dibenzene"),
    ],
)
def test_linker_stereodescriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_linear_phane_is_named_as_a_phane():
    assert (
        smiles_to_iupac("c1ccccc1Oc1cccc(Oc2cccc(Oc3ccccc3)c2)c1")
        == "2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane"
    )


def test_substituted_fused_unit_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Clc1ccc2ccccc2c1Cc1c(Cl)ccc2ccccc12")


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
