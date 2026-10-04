import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._fusion_numbering_hex import _periphery, hex_numberings


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccc2c(c1)ccc1c3ccccc3ccc21", "1,2,3,4,4a,4b,5,6,6a,7,8,9,10,10a,10b,11,12,12a"),
        ("c1ccc2cc3cc4ccccc4cc3cc2c1", "1,2,3,4,4a,5,5a,6,6a,7,8,9,10,10a,11,11a,12,12a"),
        ("c1ccc2cc3c(ccc4ccccc43)cc2c1", "1,2,3,4,4a,5,6,6a,7,7a,8,9,10,11,11a,12,12a,12b"),
        ("c1ccc2c(c1)c1ccccc1c1ccccc21", "1,2,3,4,4a,4b,5,6,7,8,8a,8b,9,10,11,12,12a,12b"),
        ("c1ccc2c(c1)ccc1ccc3ccccc3c12", "1,2,3,4,4a,5,6,6a,7,8,8a,9,10,11,12,12a,12b,12c"),
        ("c1cc2ccc3cccc4ccc(c1)c2c34", "1,2,3,3a,4,5,5a,6,7,8,8a,9,10,10a"),
        ("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "1,2,3,3a,4,5,6,6a,6b,7,8,9,9a,10,11,12,12a,12b"),
    ],
)
def test_peripheral_numbering_of_hexagonal_systems(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    numbering = hex_numberings(mol)[0]
    cycle, _ = _periphery(mol)
    walk = [numbering[a] for a in cycle]
    start = walk.index("1")
    assert ",".join(walk[start:] + walk[:start]) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Oc1ccc2ccc3cccc4ccc1c2c34", "pyren-1-ol"),
        ("Oc1cc2ccc3cccc4ccc(c1)c2c34", "pyren-2-ol"),
        ("Cc1cc2ccccc2c2ccc3ccccc3c12", "5-methylchrysene"),
        ("Oc1cc2c(ccc3ccccc32)c2ccccc12", "chrysen-6-ol"),
        ("Cc1c2ccccc2c(C)c2c1ccc1ccccc12", "7,12-dimethylbenzo[a]anthracene"),
        ("Oc1ccc2c(c1)c1ccccc1c1ccccc21", "triphenylen-2-ol"),
        ("Oc1ccc2cc3cc4ccccc4cc3cc2c1", "tetracen-2-ol"),
        ("Oc1ccc2c(c1)c1cccc3cccc2c13", "fluoranthen-8-ol"),
        ("O=C1Cc2cccc3cccc1c23", "acenaphthylen-1(2H)-one"),
        ("Oc1cc2c3c(N)cccc3cc3ccc4cccc1c4c32", "10-aminobenzo[a]pyren-12-ol"),
    ],
)
def test_substituted_larger_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccc2c(c1)C=Cc1ccccc1N2", "5H-dibenzo[b,f]azepine"),
        ("NC(=O)N1c2ccccc2C=Cc2ccccc21", "5H-dibenzo[b,f]azepine-5-carboxamide"),
        ("CN(C)CCCN1c2ccccc2CCc2ccccc21", "3-(10,11-dihydro-5H-dibenzo[b,f]azepin-5-yl)-N,N-dimethylpropan-1-amine"),
        ("CN(C)CCC=C1c2ccccc2CCc2ccccc12", "3-(10,11-dihydro-5H-dibenzo[a,d][7]annulen-5-ylidene)-N,N-dimethylpropan-1-amine"),
        ("O=C1c2ccccc2CCc2ccccc12", "10,11-dihydro-5H-dibenzo[a,d][7]annulen-5-one"),
        ("c1ccc2c(c1)C=Cc1ccccc1O2", "dibenzo[b,f]oxepine"),
    ],
)
def test_seven_membered_fused_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("NC(=O)N1CCCC1", "pyrrolidine-1-carboxamide"),
        ("OC(=O)N1CCOCC1", "morpholine-4-carboxylic acid"),
        ("NC(=O)n1c2ccccc2cc1", "1H-indole-1-carboxamide"),
        ("NC(=O)N1c2ccccc2Sc2ccccc21", "10H-phenothiazine-10-carboxamide"),
    ],
)
def test_acyl_groups_on_a_ring_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
