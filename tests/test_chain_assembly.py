import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Oc1ccc(cc1)-c1ccc(cc1)-c1ccccc1", "[11,21:24,31-terphenyl]-14-ol"),
        ("Oc1ccc(cc1)-c1cccc(c1)-c1ccccc1", "[11,21:23,31-terphenyl]-14-ol"),
        ("Nc1ccc(-c2ccc(-c3ccc(N)cc3)cc2)cc1", "[11,21:24,31-terphenyl]-14,34-diamine"),
        ("OC(=O)c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "[11,21:24,31-terphenyl]-14-carboxylic acid"),
        ("Cc1ccc(-c2ccc(-c3ccc(C)cc3)cc2)cc1", "14,34-dimethyl-11,21:24,31-terphenyl"),
        ("OCc1ccccc1-c1ccccc1-c1ccccc1", "([11,21:22,31-terphenyl]-12-yl)methanol"),
        ("CC(O)c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "1-([11,21:24,31-terphenyl]-14-yl)ethanol"),
        ("OC1CCC(CC1)C1CCC(CC1)C1CCCCC1", "[11,21:24,31-tercyclohexan]-14-ol"),
        ("Oc1ccc(nc1)-c1ccc(nc1)-c1ccccn1", "[12,22:25,32-terpyridin]-35-ol"),
        ("CCOc1ccc(-c2ccc(-c3ccc(CC(=O)O)cc3)cc2)cc1", "2-(34-ethoxy-[11,21:24,31-terphenyl]-14-yl)ethanoic acid"),
    ],
)
def test_unbranched_assemblies_with_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccc(cc1)-c1cc(cc(c1)-c1ccccc1)-c1ccccc1", "25-phenyl-11,21:23,31-terphenyl"),
        ("Oc1cc(-c2ccccc2)cc(-c2ccccc2)c1", "[11,21:23,31-terphenyl]-25-ol"),
        ("Oc1ccc(-c2ccccc2)cc1-c1ccccc1", "[11,21:23,31-terphenyl]-24-ol"),
    ],
)
def test_branched_assemblies_use_the_longest_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_ring_chain_still_is_named():
    assert smiles_to_iupac("c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccccc1C") == '12-methyl-1,7(1),2,3,4,5,6(1,4)-heptabenzenaheptaphane'


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCC(CC1)=C1CCCCC1", "1,1'-bi(cyclohexylidene)"),
        ("OC1CCCCC1=C1CCCCC1", "[1,1'-bi(cyclohexylidene)]-2-ol"),
        ("C1CCC(CC1)=C1CCCCC1C", "2-methyl-1,1'-bi(cyclohexylidene)"),
        ("C1=CCCCC1C1CCCC=C1", "1,1'-bi(cyclohex-2-ene)"),
    ],
)
def test_double_bond_junction_and_unsaturated_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_rings_that_number_differently_are_not_an_assembly():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CC=CCC1C1CC=CCC1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O[C@@H]1CCCC[C@H]1C1CCCCC1", "(1S,2R)-[1,1'-bi(cyclohexan)]-2-ol"),
        ("O[C@@H]1CCCC[C@@H]1C1CCCCC1", "(1R,2R)-[1,1'-bi(cyclohexan)]-2-ol"),
    ],
)
def test_stereodescriptors_in_a_ring_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
