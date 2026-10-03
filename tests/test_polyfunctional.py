import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCOc1ccccc1", "2-phenoxyethanol"),
        ("OCCOCCO", "2,2'-oxydi(ethan-1-ol)"),
        ("OC(=O)COCC(O)=O", "2,2'-oxydiethanoic acid"),
        ("CC(=O)CCOCCC(C)=O", "4,4'-oxydi(butan-2-one)"),
        ("OCCN(CCO)CCO", "2,2',2''-nitrilotri(ethan-1-ol)"),
        ("OCC(Cl)COCC(Cl)CO", "3,3'-oxybis(2-chloropropan-1-ol)"),
        ("OCCOCCCO", "3-(2-hydroxyethoxy)propan-1-ol"),
        ("OCCSC", "2-(methylsulfanyl)ethanol"),
        ("OCCN(C)C", "2-(dimethylamino)ethanol"),
        ("OCCNc1ccccc1", "2-anilinoethanol"),
        ("CC(C)OCCO", "2-(propan-2-yloxy)ethanol"),
        ("OCCOCCOC", "2-(2-methoxyethoxy)ethanol"),
        ("CC(C)(C)OCCO", "2-tert-butoxyethanol"),
        ("OCC[N+](=O)[O-]", "2-nitroethanol"),
        ("OCCC#N", "3-hydroxypropanenitrile"),
        ("NCCOCCN", "2,2'-oxydi(ethan-1-amine)"),
        ("COCC(=O)O", "2-methoxyethanoic acid"),
        ("OC(=O)CCC(=O)C", "4-oxopentanoic acid"),
        ("OCCCC(N)C(O)=O", "2-amino-5-hydroxypentanoic acid"),
        ("CC(O)C(=O)C", "3-hydroxybutan-2-one"),
        ("OCC(Cl)CBr", "3-bromo-2-chloropropan-1-ol"),
        ("O=CCCO", "3-hydroxypropanal"),
        ("NC(CCc1ncccc1Cl)C(=O)O", "2-amino-4-(3-chloropyridin-2-yl)butanoic acid"),
    ],
)
def test_chain_parent_with_heteroatom_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C[13CH2][18OH]",
        "C1CCCCC1OC#N",
        "N#CSCCSC#N",
        "CC=C=O",
        "CC#[N+][N-]C",
        "CC(N)N(C)C",
        "C=O",
        "OCCSSCCO",
    ],
)
def test_shapes_the_chain_engine_cannot_name_are_rejected_not_misnamed(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COc1ccc(O)cc1", "4-methoxyphenol"),
        ("Oc1ccc(cc1)[N+](=O)[O-]", "4-nitrophenol"),
        ("Nc1ccc(cc1)C(O)=O", "4-aminobenzoic acid"),
        ("Oc1ccccc1C(O)=O", "2-hydroxybenzoic acid"),
        ("COc1ccc(cc1)C=O", "4-methoxybenzaldehyde"),
        ("Nc1ccc(N)cc1", "benzene-1,4-diamine"),
        ("OC1CCC(CC1)N", "4-aminocyclohexan-1-ol"),
        ("OC(=O)C1CCC(CC1)N", "4-aminocyclohexane-1-carboxylic acid"),
        ("OC(=O)c1cccnc1", "pyridine-3-carboxylic acid"),
        ("Oc1cccnc1", "pyridin-3-ol"),
        ("CC(C)c1ccc(O)cc1C", "3-methyl-4-(propan-2-yl)phenol"),
        ("OC1CCCCC1c1ccccc1", "2-phenylcyclohexan-1-ol"),
        ("CC(C)Cc1ccc(cc1)C(C)C(O)=O", "2-[4-(2-methylpropyl)phenyl]propanoic acid"),
    ],
)
def test_ring_parent_with_heteroatom_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1(CCCC1)C1CCCC1", "[1,1'-bi(cyclopentan)]-1-ol"),
        ("Oc1ccc(cc1)-c1ccccc1", "[1,1'-biphenyl]-4-ol"),
        ("OC(=O)c1ccc(cc1)-c1ccccc1", "[1,1'-biphenyl]-4-carboxylic acid"),
        ("Nc1ccc(cc1)-c1ccccc1", "[1,1'-biphenyl]-4-amine"),
        ("Oc1ccc(cc1)-c1ccc(O)cc1", "[1,1'-biphenyl]-4,4'-diol"),
        ("Clc1ccc(cc1)-c1ccc(O)cc1", "4'-chloro-[1,1'-biphenyl]-4-ol"),
        ("OC1CCC(CC1)C1CCCCC1", "[1,1'-bi(cyclohexan)]-4-ol"),
        ("COc1ccc(cc1)-c1ccccc1", "4-methoxy-1,1'-biphenyl"),
        ("OC(=O)c1ccc(cc1)-c1ccc(cc1)C(O)=O", "[1,1'-biphenyl]-4,4'-dicarboxylic acid"),
    ],
)
def test_identical_rings_joined_directly_form_a_ring_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Clc1ccc(cc1)[N+](=O)[O-]", "1-chloro-4-nitrobenzene"),
        ("COCCOC", "1,2-dimethoxyethane"),
        ("COc1ccccc1OC", "1,2-dimethoxybenzene"),
        ("C1CCC([N+](=O)[O-])CC1", "nitrocyclohexane"),
        ("COC1CCCCC1", "methoxycyclohexane"),
        ("COC=C", "methoxyethene"),
        ("CSc1ccc(Cl)cc1", "1-chloro-4-(methylsulfanyl)benzene"),
        ("CN(C)C1CCCCC1", "N,N-dimethylcyclohexanamine"),
        ("c1ccccc1CN(C)C1CCCCC1", "N-benzyl-N-methylcyclohexanamine"),
        ("CN(CC)c1ccc(Cl)cc1", "4-chloro-N-ethyl-N-methylaniline"),
        ("CC(C)NC(C)C", "N-(propan-2-yl)propan-2-amine"),
    ],
)
def test_parents_without_a_principal_group_and_n_substituted_amines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COc1cc(C=O)ccc1O", "4-hydroxy-3-methoxybenzaldehyde"),
        ("C=CCc1ccc(O)c(OC)c1", "2-methoxy-4-(prop-2-en-1-yl)phenol"),
        ("NC(=O)c1ccccc1O", "2-hydroxybenzamide"),
        ("Oc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]", "2,4-dinitrophenol"),
        ("CC(C)c1ccc(C)cc1O", "5-methyl-2-(propan-2-yl)phenol"),
        ("CC(C)C1CCC(C)CC1O", "5-methyl-2-(propan-2-yl)cyclohexan-1-ol"),
        ("OCc1ccccc1", "phenylmethanol"),
        ("OCCNCCO", "2,2'-azanediyldi(ethan-1-ol)"),
        ("COC(=O)CC(=O)O", "methyl hydrogen propanedioate"),
        ("COC(=O)c1ccccc1C(O)=O", "methyl hydrogen benzene-1,2-dicarboxylate"),
    ],
)
def test_known_compounds_through_the_fallback_engines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCC(=O)NC", "3-hydroxy-N-methylpropanamide"),
        ("CC(=O)Oc1ccccc1C(O)=O", "2-(ethanoyloxy)benzoic acid"),
        ("OC(=O)CC(=O)NC", "2-(methylcarbamoyl)ethanoic acid"),
    ],
)
def test_prefixes_and_n_substituents_are_cited_alphabetically(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Nc1ccc(O)cc1", "N-(4-hydroxyphenyl)ethanamide"),
        ("O=C(Nc1ccccc1)c1ccccc1", "N-phenylbenzamide"),
        ("CC(=O)N(C)CCO", "N-(2-hydroxyethyl)-N-methylethanamide"),
        ("CNC(=O)c1ccccc1", "N-methylbenzamide"),
        ("CC(=O)NC(C)C(O)=O", "2-(ethanoylamino)propanoic acid"),
    ],
)
def test_n_substituted_amide_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@@H](N)CO", "(2R)-2-aminopropan-1-ol"),
        ("O[C@H]1CCCC[C@@H]1N", "(1S,2S)-2-aminocyclohexan-1-ol"),
        ("COC[C@@H](C)C(C)=O", "(3R)-4-methoxy-3-methylbutan-2-one"),
        ("C[C@@H](O)[C@H](N)C(=O)O", "(2S,3R)-2-amino-3-hydroxybutanoic acid"),
        ("C[C@@H]1C[C@H](C)CCC1", "(1R,3S)-1,3-dimethylcyclohexane"),
    ],
)
def test_stereodescriptors_on_the_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereodescriptors_on_a_substituent_are_rejected_not_dropped_is_named():
    assert smiles_to_iupac("Oc1ccc(cc1)[C@H](C)Cl") == "4-[(1S)-1-chloroethyl]phenol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)CC(O)(CC(O)=O)C(O)=O", "2-hydroxypropane-1,2,3-tricarboxylic acid"),
        ("N#CC(CC#N)CC#N", "propane-1,2,3-tricarbonitrile"),
        ("OC(=O)CC(CC(O)=O)C(O)=O", "propane-1,2,3-tricarboxylic acid"),
        ("OC(=O)CC(=CC(O)=O)C(O)=O", "prop-1-ene-1,2,3-tricarboxylic acid"),
        ("NC(=O)CC(CC(N)=O)C(N)=O", "propane-1,2,3-tricarboxamide"),
    ],
)
def test_carbo_suffix_when_the_groups_do_not_fit_one_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCc1ccc(CCO)cc1", "2,2'-(1,4-phenylene)di(ethan-1-ol)"),
        ("OC(=O)Cc1ccc(CC(O)=O)cc1", "2,2'-(1,4-phenylene)diethanoic acid"),
        ("OC(=O)CC1CCC(CC(O)=O)CC1", "2,2'-(cyclohexane-1,4-diyl)diethanoic acid"),
        ("OCC(Cl)c1ccc(C(Cl)CO)cc1", "2,2'-(1,4-phenylene)bis(2-chloroethan-1-ol)"),
        ("OCc1ccccc1CO", "(1,2-phenylene)dimethanol"),
        ("NCc1ccc(CN)cc1", "(1,4-phenylene)dimethanamine"),
        ("OCOCO", "oxydimethanol"),
    ],
)
def test_chain_units_on_a_ring_linker_are_named_multiplicatively(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCOCCCO", "3-(2-hydroxyethoxy)propan-1-ol"),
        ("OCC(CO)(CO)CO", "2,2-bis(hydroxymethyl)propane-1,3-diol"),
        ("OC(=O)CCCC(O)=O", "pentanedioic acid"),
    ],
)
def test_chains_that_hold_every_group_stay_substitutive(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Oc1ccc2ccccc2c1", "naphthalen-2-ol"),
        ("Nc1cccc2ccccc12", "naphthalen-1-amine"),
        ("Nc1ccc2ccc(N)cc2c1", "naphthalene-2,7-diamine"),
        ("Oc1ccc2ccccc2c1O", "naphthalene-1,2-diol"),
        ("O=C(O)c1ccc2ccccc2c1", "naphthalene-2-carboxylic acid"),
        ("OC(=O)c1ccc2ccccc2c1O", "1-hydroxynaphthalene-2-carboxylic acid"),
        ("NC(=O)c1cccc2ccccc12", "naphthalene-1-carboxamide"),
        ("N#Cc1ccc2ccccc2c1", "naphthalene-2-carbonitrile"),
        ("O=Cc1ccc2ccccc2c1", "naphthalene-2-carbaldehyde"),
        ("OC(=O)c1c[nH]c2ccccc12", "1H-indole-3-carboxylic acid"),
        ("OC(=O)c1cc2ccccc2[nH]1", "1H-indole-2-carboxylic acid"),
        ("Oc1cccc2cccnc12", "quinolin-8-ol"),
        ("Oc1ccc2cc[nH]c2c1", "1H-indol-6-ol"),
        ("Oc1c2ccccc2cc2ccccc12", "anthracen-9-ol"),
        ("Oc1cc2ccccc2c2ccccc12", "phenanthren-9-ol"),
        ("COc1ccc2ccccc2c1", "2-methoxynaphthalene"),
    ],
)
def test_fused_aromatic_parents_with_functional_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Nc1ccc(cc1)S(N)(=O)=O", "4-aminobenzene-1-sulfonamide"),
        ("OS(=O)(=O)c1ccc(N)cc1", "4-aminobenzene-1-sulfonic acid"),
        ("OS(=O)(=O)CCN", "2-aminoethanesulfonic acid"),
        ("CN(C)S(=O)(=O)c1ccccc1N", "2-amino-N,N-dimethylbenzene-1-sulfonamide"),
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
        ("OCCN1CCCC1", "2-(pyrrolidin-1-yl)ethanol"),
        ("OCCN1CCOCC1", "2-(morpholin-4-yl)ethanol"),
        ("OCCN1CCN(C)CC1", "2-(4-methylpiperazin-1-yl)ethanol"),
        ("OCC1CCNCC1", "(piperidin-4-yl)methanol"),
        ("OCCc1ccc2ccccc2c1", "2-(naphthalen-2-yl)ethanol"),
        ("NCCc1c[nH]c2ccccc12", "2-(1H-indol-3-yl)ethanamine"),
        ("Oc1ccc(cc1)N1CCOCC1", "4-(morpholin-4-yl)phenol"),
        ("NC(CCc1ncccc1Cl)C(=O)O", "2-amino-4-(3-chloropyridin-2-yl)butanoic acid"),
    ],
)
def test_heterocyclic_and_fused_ring_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1CCNCC1", "piperidin-4-ol"),
        ("OC(=O)C1CCNCC1", "piperidine-4-carboxylic acid"),
        ("O=C1CCCN1", "pyrrolidin-2-one"),
        ("OC1C=CCCC1", "cyclohex-2-en-1-ol"),
        ("OCCN1CCN(CCO)CC1", "2,2'-(piperazine-1,4-diyl)di(ethan-1-ol)"),
        ("Nc1ccc(cc1)N1CCCC1", "4-(pyrrolidin-1-yl)aniline"),
        ("Clc1ccc(cc1)C1CCNCC1", "4-(4-chlorophenyl)piperidine"),
        ("Cc1ccc(cc1)N1CCOCC1", "4-(4-methylphenyl)morpholine"),
        ("c1ccc(cc1)-c1ccccn1", "2-phenylpyridine"),
    ],
)
def test_heterocyclic_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "S=C1CCCC=C1C",
        "OCCC1CCC2CCCCC2C1",
    ],
)
def test_ring_assembly_substituents_thioketones_and_saturated_fused_rings_are_rejected(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCc1ccc(cc1)-c1ccccc1", "2-([1,1'-biphenyl]-4-yl)ethanol"),
        ("OC(=O)Cc1ccc(cc1)-c1ccccc1", "2-([1,1'-biphenyl]-4-yl)ethanoic acid"),
        ("Clc1ccc(cc1)-c1ccc(CCO)cc1", "2-(4'-chloro-[1,1'-biphenyl]-4-yl)ethanol"),
        ("OCCC1CCC(CC1)C1CCCCC1", "2-([1,1'-bi(cyclohexan)]-4-yl)ethanol"),
        ("NCc1ccccc1-c1ccccc1", "([1,1'-biphenyl]-2-yl)methanamine"),
        ("c1ccc(-c2ccc(-c3ccccn3)cc2)cc1", "2-([1,1'-biphenyl]-4-yl)pyridine"),
    ],
)
def test_ring_assembly_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected




@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)(C)OC(=O)NCC(=O)O", "2-[(tert-butoxycarbonyl)amino]ethanoic acid"),
        ("NCC(=O)NCC(=O)O", "2-[(2-aminoethanoyl)amino]ethanoic acid"),
        ("O=C(O)CNC(=O)c1ccccc1", "2-(benzoylamino)ethanoic acid"),
        ("OC(=O)CNC(=O)OCc1ccccc1", "2-[(benzyloxycarbonyl)amino]ethanoic acid"),
        ("NCC(=O)Oc1ccc(cc1)C(O)=O", "4-[(2-aminoethanoyl)oxy]benzoic acid"),
        ("CC(=O)Nc1ccc(cc1)C(O)=O", "4-(ethanoylamino)benzoic acid"),
    ],
)
def test_acyl_prefixes_with_substituents_and_alkoxycarbonylamino(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Oc1ccc(cc1)[C@H](C)Cl", "4-[(1S)-1-chloroethyl]phenol"),
        ("OC[C@H]1CCCC[C@@H]1C", "[(1S,2S)-2-methylcyclohexyl]methanol"),
        ("OC(=O)CC[C@H]1CCCC[C@@H]1C", "3-[(1R,2S)-2-methylcyclohexyl]propanoic acid"),
        ("NC[C@H](Cl)c1ccc(O)cc1", "4-[(1R)-2-amino-1-chloroethyl]phenol"),
        ("c1ccccc1O[C@H](C)CC", "{[(2R)-butan-2-yl]oxy}benzene"),
        ("OCCO[C@H](C)CC", "2-{[(2R)-butan-2-yl]oxy}ethanol"),
    ],
)
def test_stereodescriptors_inside_substituent_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)COCCOCCOCCOCC(O)=O", "3,6,9,12-tetraoxatetradecanedioic acid"),
        ("OCCOCCOCCOCCOCCO", "3,6,9,12-tetraoxatetradecane-1,14-diol"),
        ("NCCOCCOCCOCCOCCN", "3,6,9,12-tetraoxatetradecane-1,14-diamine"),
        ("NCCNCCNCCNCCNCCN", "3,6,9,12-tetraazatetradecane-1,14-diamine"),
        ("ClCCOCCOCCOCCOCCCl", "1,14-dichloro-3,6,9,12-tetraoxatetradecane"),
        ("NCCNCCNCCNCCOCCCC(O)CCCC", "1-amino-12-oxa-3,6,9-triazaicosan-16-ol"),
        ("COCCSCCOCCNCC", "2,8-dioxa-5-thia-11-azatridecane"),
        ("OCCOCCO", "2,2'-oxydi(ethan-1-ol)"),
    ],
)
def test_skeletal_replacement_parents_with_four_or_more_heterounits(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCc1ccc(nc1)-c1ccccn1", "([2,2'-bipyridin]-5-yl)methanol"),
        ("OCc1ccc(nc1)-c1ccc(Cl)cn1", "(5'-chloro-[2,2'-bipyridin]-5-yl)methanol"),
        ("OCc1ccc(s1)-c1cccs1", "([2,2'-bithiophen]-5-yl)methanol"),
        ("c1ccc(-c2ccc(-c3ccco3)s2)s1", "2-([2,2'-bithiophen]-5-yl)furan"),
        ("Oc1ccc(nc1)-c1ccccn1", "[2,2'-bipyridin]-5-ol"),
        ("OC(=O)c1ccc(nc1)-c1ccccn1", "[2,2'-bipyridine]-5-carboxylic acid"),
        ("Clc1ccc(nc1)-c1ccccn1", "5-chloro-2,2'-bipyridine"),
        ("Oc1cc(O)c(c(c1)O)-c1ccc(O)cc1", "[1,1'-biphenyl]-2,4,4',6-tetrol"),
        ("Cc1ccc(cc1)-c1ccc(C)cc1O", "4,4'-dimethyl-[1,1'-biphenyl]-2-ol"),
    ],
)
def test_heteroaromatic_assemblies_and_primed_locant_order(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC[Si](C)(C)C", "2-(trimethylsilyl)ethanol"),
        ("OCC[Si](C)(C)c1ccccc1", "2-[dimethyl(phenyl)silyl]ethanol"),
        ("OCCP(c1ccccc1)c1ccccc1", "2-(diphenylphosphanyl)ethanol"),
        ("OCCB(C)C", "2-(dimethylboranyl)ethanol"),
        ("C[Si](C)(C)OCCO", "2-[(trimethylsilyl)oxy]ethanol"),
        ("OCCN[Si](C)(C)C", "2-[(trimethylsilyl)amino]ethanol"),
        ("OCCNN", "2-hydrazinylethanol"),
        ("OC(=O)CC[Si](C)(C)C", "3-(trimethylsilyl)propanoic acid"),
        ("c1ccccc1[Si](C)(C)C", "trimethyl(phenyl)silane"),
        ("Clc1ccc(cc1)P(c1ccc(F)cc1)", "(4-chlorophenyl)(4-fluorophenyl)phosphane"),
    ],
)
def test_mononuclear_hydride_prefixes_and_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCC12CC3CC(CC(C3)C1)C2", "(adamantan-1-yl)methanol"),
        ("OCC1C2CC3CC(C2)CC1C3", "(adamantan-2-yl)methanol"),
        ("OC(=O)CC12CC3CC(CC(C3)C1)C2", "2-(adamantan-1-yl)ethanoic acid"),
        ("OC12CC3CC(CC(C3)C1)C2", "adamantan-1-ol"),
        ("C1C2CC3CC1CC(C2)C3", "adamantane"),
        ("Cc1ccc(cc1)C(c1ccccc1)(c1ccccc1)c1ccccn1", "2-[(4-methylphenyl)di(phenyl)methyl]pyridine"),
        ("c1ccc(cc1)C(c1ccccc1)(c1ccccc1)c1ccccn1", "2-(triphenylmethyl)pyridine"),
    ],
)
def test_retained_adamantane_and_one_carbon_prefix_multiplication(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCc1ccncn1", "(pyrimidin-4-yl)methanol"),
        ("OCc1cnccn1", "(pyrazin-2-yl)methanol"),
        ("OCc1ncccn1", "(pyrimidin-2-yl)methanol"),
        ("OCc1nnccc1", "(pyridazin-3-yl)methanol"),
        ("OCc1ncncn1", "(1,3,5-triazin-2-yl)methanol"),
        ("OCc1ccc(Cl)cn1", "(5-chloropyridin-2-yl)methanol"),
    ],
)
def test_heteroaromatic_substituents_are_not_read_as_phenyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CCCCCC(CCC(C)C)CCCCC", "[6-(3-methylbutyl)undecyl]benzene"),
        ("CC(C)CCC(CCCCCC1CCCCC1)CCCCC", "[6-(3-methylbutyl)undecyl]cyclohexane"),
        ("c1ccncc1CCCCCC(CCC(C)C)CCCCC", "3-[6-(3-methylbutyl)undecyl]pyridine"),
    ],
)
def test_nested_enclosing_marks_escalate_in_ring_parent_modules(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
