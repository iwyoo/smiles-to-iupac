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
        "N#CC(CC#N)CC#N",
        "OC(=O)CC(O)(CC(O)=O)C(O)=O",
        "CC(N)N(C)C",
        "C=O",
        "OCC[Si](C)(C)C",
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


def test_identical_rings_joined_directly_are_not_named_as_a_substituted_ring():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1(CCCC1)C1CCCC1")


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
