import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1Cc1ccccc1", "1,1'-methylenedibenzene"),
        ("c1ccccc1CCc1ccccc1", "1,1'-(ethane-1,2-diyl)dibenzene"),
        ("c1ccccc1Oc1ccccc1", "1,1'-oxydibenzene"),
        ("C1CCCCC1CCC1CCCCC1", "1,1'-(ethane-1,2-diyl)dicyclohexane"),
        ("c1ccccc1Sc1ccccc1", "1,1'-sulfanediyldibenzene"),
        ("c1ccccc1OOc1ccccc1", "1,1'-peroxydibenzene"),
        ("c1ccccc1SSc1ccccc1", "1,1'-disulfanediyldibenzene"),
        ("c1ccccc1S(=O)(=O)c1ccccc1", "1,1'-sulfonyldibenzene"),
        ("c1ccccc1C#Cc1ccccc1", "1,1'-(ethyne-1,2-diyl)dibenzene"),
        ("c1ccc2ccccc2c1Cc1cccc2ccccc12", "1,1'-methylenedinaphthalene"),
        ("c1ccncc1Cc1cccnc1", "3,3'-methylenedipyridine"),
        ("c1cc[nH]c1Cc1ccc[nH]1", "2,2'-methylenedi(1H-pyrrole)"),
    ],
)
def test_unsubstituted_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Clc1ccc(Cc2ccc(Cl)cc2)cc1", "1,1'-methylenebis(4-chlorobenzene)"),
        ("Brc1ccc(Oc2ccc(Br)cc2)cc1", "1,1'-oxybis(4-bromobenzene)"),
        ("ClCc1ccc(Cc2ccc(CCl)c(Br)c2)cc1Br", "1,1'-methylenebis[3-bromo-4-(chloromethyl)benzene]"),
        ("CC(Cl)c1cncc(Oc2cncc(C(C)Cl)c2)c1", "3,3'-oxybis[5-(1-chloroethyl)pyridine]"),
        ("Cc1ccc(Cc2ccc(C)nc2)cn1", "3,3'-methylenebis(6-methylpyridine)"),
        ("Cn1cccc1Cc1cccn1C", "2,2'-methylenebis(1-methyl-1H-pyrrole)"),
        ("CC(C)(C)c1ccc(Cc2ccc(C(C)(C)C)cc2)cc1", "1,1'-methylenebis(4-tert-butylbenzene)"),
        ("COc1ccc(Cc2ccc(OC)cc2)cc1", "1,1'-methylenebis(4-methoxybenzene)"),
    ],
)
def test_substituted_units_and_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(Oc2ccc(C(=O)O)cc2)cc1", "4,4'-oxydibenzoic acid"),
        ("OC(=O)c1ccc(Oc2ccc(C(=O)O)c(Br)c2)cc1Br", "4,4'-oxybis(2-bromobenzoic acid)"),
        ("OC(=O)C1CCC(CC1)OC1CCC(CC1)C(=O)O", "4,4'-oxydi(cyclohexane-1-carboxylic acid)"),
        ("OS(=O)(=O)c1ccc(Oc2ccc(S(=O)(=O)O)cc2)cc1", "4,4'-oxydi(benzene-1-sulfonic acid)"),
        ("Oc1ccc(Cc2ccc(O)cc2)cc1", "4,4'-methylenediphenol"),
        ("Nc1ccc(Cc2ccc(N)cc2)cc1", "4,4'-methylenedianiline"),
        ("OC1CCC(CC1)CC1CCC(O)CC1", "4,4'-methylenedi(cyclohexan-1-ol)"),
        ("OC(=O)c1ccc(Cc2ccc(C(=O)O)c(C(=O)O)c2)cc1C(=O)O", "4,4'-methylenedi(benzene-1,2-dicarboxylic acid)"),
        ("OC(=O)c1ccccc1Cc1ccccc1C(=O)O", "2,2'-methylenedibenzoic acid"),
        ("OC(=O)c1ccc(C(=O)c2ccc(C(=O)O)cc2)cc1", "4,4'-carbonyldibenzoic acid"),
        ("OC(=O)c1ccc(Nc2ccc(C(=O)O)cc2)cc1", "4,4'-azanediyldibenzoic acid"),
        ("OC(=O)c1ccc(N(C)c2ccc(C(=O)O)cc2)cc1", "4,4'-(methylazanediyl)dibenzoic acid"),
        ("OC(=O)c1ccc(C(O)c2ccc(C(=O)O)cc2)cc1", "4,4'-(hydroxymethylene)dibenzoic acid"),
        ("O=C1CCC(CC1)CC1CCC(=O)CC1", "4,4'-methylenedi(cyclohexan-1-one)"),
        ("N#Cc1ccc(Cc2ccc(C#N)cc2)cc1", "4,4'-methylenedibenzonitrile"),
    ],
)
def test_functionalized_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1C(c1ccccc1)c1ccccc1", "1,1',1''-methanetriyltribenzene"),
        ("C(c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1", "1,1',1'',1'''-methanetetrayltetrabenzene"),
        ("c1ccccc1CC(c1ccccc1)c1ccccc1", "1,1',1''-(ethane-1,1,2-triyl)tribenzene"),
        ("OC(=O)c1ccc(cc1)C(c1ccc(cc1)C(=O)O)Cc1ccc(cc1)C(=O)O", "4,4',4''-(ethane-1,1,2-triyl)tribenzoic acid"),
        ("OC(=O)c1ccc(N(c2ccc(C(=O)O)cc2)c2ccc(C(=O)O)cc2)cc1", "4,4',4''-nitrilotribenzoic acid"),
        ("NC1=CC=C(C=C1)C=C(C1=CC=C(N)C=C1)C1=CC=C(N)C=C1", "4,4',4''-(ethene-1,1,2-triyl)trianiline"),
        ("c1ccccc1Cc1cc(Cc2ccccc2)cc(Cc2ccccc2)c1", "1,1',1''-[benzene-1,3,5-triyltris(methylene)]tribenzene"),
    ],
)
def test_three_and_more_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1C(C)(C)c1ccccc1", "1,1'-(propane-2,2-diyl)dibenzene"),
        ("c1ccccc1C(C)c1ccccc1", "1,1'-(ethane-1,1-diyl)dibenzene"),
        ("c1ccccc1CC(C)c1ccccc1", "1,1'-(propane-1,2-diyl)dibenzene"),
        ("c1ccccc1C(Cl)c1ccccc1", "1,1'-(chloromethylene)dibenzene"),
        ("c1ccccc1CC(Cl)Cc1ccccc1", "1,1'-(2-chloropropane-1,3-diyl)dibenzene"),
        ("BrC1=CC=C(C=C1)C=CCc1ccc(Br)cc1", "1,1'-(prop-1-ene-1,3-diyl)bis(4-bromobenzene)"),
        ("c1ccccc1CC=CCc1ccccc1", "1,1'-(but-2-ene-1,4-diyl)dibenzene"),
        ("c1ccccc1C(C1CCCCC1)c1ccccc1", "1,1'-(cyclohexylmethylene)dibenzene"),
        ("c1ccccc1C1CC1c1ccccc1", "1,1'-(cyclopropane-1,2-diyl)dibenzene"),
        ("c1ccccc1CC(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1", "1,1'-(2,2-dibenzylpropane-1,3-diyl)dibenzene"),
        ("c1ccccc1C(SCc1ccccc1)SCc1ccccc1", "1,1'-[(phenylmethylene)bis(sulfanediylmethylene)]dibenzene"),
    ],
)
def test_central_group_forms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1COCc1ccccc1", "1,1'-[oxybis(methylene)]dibenzene"),
        ("c1ccccc1OCCOc1ccccc1", "1,1'-[ethane-1,2-diylbis(oxy)]dibenzene"),
        ("c1ccccc1OCOc1ccccc1", "1,1'-[methylenebis(oxy)]dibenzene"),
        ("c1ccccc1OCCOCCOc1ccccc1", "1,1'-[oxybis(ethane-2,1-diyloxy)]dibenzene"),
        ("c1ccccc1COCCOCc1ccccc1", "1,1'-[ethane-1,2-diylbis(oxymethylene)]dibenzene"),
        ("Oc1ccc(OCC(C)COc2ccc(O)cc2)cc1", "4,4'-[(2-methylpropane-1,3-diyl)bis(oxy)]diphenol"),
        ("OC(=O)C1CCCCC1COCCOCC1CCCCC1C(=O)O", "2,2'-[ethane-1,2-diylbis(oxymethylene)]di(cyclohexane-1-carboxylic acid)"),
        ("Oc1ccc(CCCCCCCCCCCCCCOCCCCCCCCCCCCCCc2ccc(O)cc2)cc1", "4,4'-[oxydi(tetradecane-14,1-diyl)]diphenol"),
        ("OS(=O)(=O)C1CCC(CC1)C(F)COCC(F)C1CCC(CC1)S(=O)(=O)O", "4,4'-[oxybis(1-fluoroethane-2,1-diyl)]di(cyclohexane-1-sulfonic acid)"),
        ("c1ccccc1Cc1ccc(Cc2ccccc2)cc1", "1,1'-[1,4-phenylenebis(methylene)]dibenzene"),
        ("c1ccccc1Oc1ccccc1Oc1ccccc1", "1,1'-[1,2-phenylenebis(oxy)]dibenzene"),
        ("C(COc1ccoc1)Oc1cocc1OCCOc1ccoc1", "3,3'-[furan-3,4-diylbis(oxyethane-2,1-diyloxy)]difuran"),
        ("c1ccccc1OCC(Cl)COc1ccccc1", "1,1'-[(2-chloropropane-1,3-diyl)bis(oxy)]dibenzene"),
        ("c1ccccc1Oc1cc(Oc2ccccc2)cc(Oc2ccccc2)c1", "1,1',1''-[benzene-1,3,5-triyltris(oxy)]tribenzene"),
    ],
)
def test_concatenated_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1/C=C/c1ccccc1", "1,1'-[(1E)-ethene-1,2-diyl]dibenzene"),
        ("c1ccccc1/C=C\\c1ccccc1", "1,1'-[(1Z)-ethene-1,2-diyl]dibenzene"),
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


def test_benzene_units_are_senior_to_cyclohexane_linker():
    assert smiles_to_iupac("c1ccccc1CC1CCC(CC1)Cc1ccccc1") == "1,1'-[cyclohexane-1,4-diylbis(methylene)]dibenzene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CCCCC1CC1C=CCCC1", "3,3'-methylenedi(cyclohex-1-ene)"),
        ("C1CC=CCC1CC1CCC=CC1", "4,4'-methylenedi(cyclohex-1-ene)"),
        ("C1=CCCCC1OC1C=CCCC1", "3,3'-oxydi(cyclohex-1-ene)"),
        ("C1=CCCCC1CCC1C=CCCC1", "3,3'-(ethane-1,2-diyl)di(cyclohex-1-ene)"),
    ],
)
def test_cycloalkene_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CCCCC1CC1CCCCC1", "3-(cyclohexylmethyl)cyclohex-1-ene"),
        ("C1=CCCCC1CC1=CCCCC1", "1-[(cyclohex-2-en-1-yl)methyl]cyclohex-1-ene"),
        ("c1ccsc1Cc1ccco1", "2-[(thiophen-2-yl)methyl]furan"),
        ("CC1CCCCC1CC1CCCCC1", "1-(cyclohexylmethyl)-2-methylcyclohexane"),
    ],
)
def test_unequal_rings_are_ranked_by_unsaturation_heteroatom_and_substituent_count(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    ["c1c[nH]c(-c2[nH]ccn2)n1", "C1=CCCCC1C1=CCCCC1C1=CCCCC1"],
)
def test_identical_rings_joined_directly_and_long_ring_chains_are_not_named_substitutively(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCc1ccc(cc1)Oc1ccc(CCO)cc1", "2,2'-[oxybis(4,1-phenylene)]di(ethan-1-ol)"),
        ("OCCc1ccc(cc1)Cc1ccc(CCO)cc1", "2,2'-[methylenebis(4,1-phenylene)]di(ethan-1-ol)"),
        ("OC(=O)Cc1ccc(cc1)Oc1ccc(CC(O)=O)cc1", "2,2'-[oxybis(4,1-phenylene)]diethanoic acid"),
        ("OC(=O)Cc1ccc(Cc2ccc(CC(=O)O)cc2)cc1", "2,2'-[methylenebis(4,1-phenylene)]diethanoic acid"),
        ("NCCc1ccc(cc1)Sc1ccc(CCN)cc1", "2,2'-[sulfanediylbis(4,1-phenylene)]di(ethan-1-amine)"),
        ("OCCc1cccc(c1)Oc1cccc(CCO)c1", "2,2'-[oxybis(3,1-phenylene)]di(ethan-1-ol)"),
        ("OCCc1cc2ccccc2cc1CCO", "2,2'-(naphthalene-2,3-diyl)di(ethan-1-ol)"),
        ("OCCc1ccc2cc(CCO)ccc2c1", "2,2'-(naphthalene-2,6-diyl)di(ethan-1-ol)"),
        ("OCCc1ccc2[nH]c(CCO)cc2c1", "2,2'-(1H-indole-2,5-diyl)di(ethan-1-ol)"),
        ("OCCN1CCN(CCO)CC1", "2,2'-(piperazine-1,4-diyl)di(ethan-1-ol)"),
    ],
)
def test_chain_units_on_fused_and_composite_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
