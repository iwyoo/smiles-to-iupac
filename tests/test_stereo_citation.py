import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H]1CC[C@@H](C)CC1", "(1s,4s)-1,4-dimethylcyclohexane"),
        ("C[C@@H]1CC[C@@H](C)CC1", "(1r,4r)-1,4-dimethylcyclohexane"),
        ("O[C@H]1CC[C@H](O)CC1", "(1r,4r)-cyclohexane-1,4-diol"),
        ("C1CC[C@H]2CCCC[C@@H]2C1", "(4ar,8ar)-decahydronaphthalene"),
        ("C1CC[C@H]2CCCC[C@H]2C1", "(4as,8as)-decahydronaphthalene"),
    ],
)
def test_pseudoasymmetric_centres_use_lower_case_descriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H](Cl)C(Cl)C(=O)O", "(3S)-2,3-dichlorobutanoic acid"),
        ("N[C@H]1CCCCC1Cl", "(1S)-2-chlorocyclohexan-1-amine"),
        ("C1CC[C@H](NC1)C(=O)O", "(2S)-piperidine-2-carboxylic acid"),
        ("CC(O)/C=C/C", "(3E)-pent-3-en-2-ol"),
        ("C[C@H]1CCC2CCCCC2C1", "(2S)-2-methyldecahydronaphthalene"),
    ],
)
def test_only_the_specified_elements_are_cited(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O[C@H]1CCC[C@@H](O)C1", "(1R,3S)-cyclohexane-1,3-diol"),
        ("O[C@@H]1CCC[C@H](O)C1", "(1R,3S)-cyclohexane-1,3-diol"),
        ("N[C@H]1CCC[C@@H](N)C1", "(1R,3S)-cyclohexane-1,3-diamine"),
        ("O[C@H]1CC[C@@H](O)C1", "(1R,3S)-cyclopentane-1,3-diol"),
        ("OC(=O)[C@H]1CCC[C@@H](C(O)=O)C1", "(1R,3S)-cyclohexane-1,3-dicarboxylic acid"),
    ],
)
def test_r_is_cited_at_the_lower_locant_of_a_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O[C@H]1CCCC[C@@H]1C1CCC(CC1)C1CCCCC1", "(11R,12S)-[11,21:24,31-tercyclohexan]-12-ol"),
        ("C[C@H]1CCCC[C@@H]1C1CCCCC1C1CCCCC1", "(11S,12S)-12-methyl-11,21:22,31-tercyclohexane"),
    ],
)
def test_stereodescriptors_in_an_assembly_of_three_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
