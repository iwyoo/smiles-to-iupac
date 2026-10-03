import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure, adjacency, halogen_substituents
from smiles_to_iupac._substituents import name_branch


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC=C1CCC(C)CC1", "1-ethylidene-4-methylcyclohexane"),
        ("CC1CCCCC1=C", "1-methyl-2-methylidenecyclohexane"),
        ("C=C1CCCCC1=C", "1,2-dimethylidenecyclohexane"),
        ("C=C1CCCCC1C#C", "1-ethynyl-2-methylidenecyclohexane"),
        ("C=C(C)C1CCC1=C", "1-methylidene-2-(prop-1-en-2-yl)cyclobutane"),
        ("ClC1CCC(=C)CC1", "1-chloro-4-methylidenecyclohexane"),
        ("ClC(Cl)=C1CCCC1", "(dichloromethylidene)cyclopentane"),
        ("CCC(CC)=C1CCCC1", "(pentan-3-ylidene)cyclopentane"),
        ("CC(C)(C)C=C1CCCC1", "(2,2-dimethylpropylidene)cyclopentane"),
        ("C=C=C1CCCCC1", "ethenylidenecyclohexane"),
        ("C=CC=C1CCCCC1", "(prop-2-en-1-ylidene)cyclohexane"),
    ],
)
def test_ring_parent_with_ylidene_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=CC1CCCCC1", "ethenylcyclohexane"),
        ("C#CC1CCCCC1", "ethynylcyclohexane"),
        ("CC=CC1CCCCC1", "(prop-1-en-1-yl)cyclohexane"),
        ("CC(=C)C1CCCCC1", "(prop-1-en-2-yl)cyclohexane"),
        ("C=CCC1CCCCC1", "(prop-2-en-1-yl)cyclohexane"),
        ("CC(C)=CC1CCCCC1", "(2-methylprop-1-en-1-yl)cyclohexane"),
        ("CCCC(=CC)C1CCCCC1", "(hex-2-en-3-yl)cyclohexane"),
        ("C=CC=CC1CCCCC1", "(buta-1,3-dien-1-yl)cyclohexane"),
        ("CC#CC#CC1CCCCC1", "(penta-1,3-diyn-1-yl)cyclohexane"),
        ("ClC=CC1CCCC1", "(2-chloroethen-1-yl)cyclopentane"),
        ("C=C(C1CCCCC1)C1CCCCC1", "1,1'-(ethene-1,1-diyl)dicyclohexane"),
    ],
)
def test_ring_parent_with_enyl_and_ynyl_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CCCCC1=CC", "3-ethylidenecyclohexene"),
        ("C=C1C=CC=C1", "5-methylidenecyclopenta-1,3-diene"),
        ("CC(=C)C1=CCCCC1", "1-(prop-1-en-2-yl)cyclohexene"),
        ("C=C1CC=CC1", "4-methylidenecyclopentene"),
    ],
)
def test_unsaturated_ring_parent_with_exocyclic_multiple_bonds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=C(CC)CC", "3-methylidenepentane"),
        ("CCCC(=C)CCC", "4-methylideneheptane"),
        ("C=C(CC)C(C)C", "2-methyl-3-methylidenepentane"),
        ("C=C(C=C)C=C", "3-methylidenepenta-1,4-diene"),
        ("CC=CC(=CC)C=CC", "4-ethylidenehepta-2,5-diene"),
        ("C=CC(C=C)CC=C", "3-ethenylhexa-1,5-diene"),
        ("CCC(=CC)CC", "3-ethylpent-2-ene"),
        ("C#CC(CC)C(=C)C", "3-ethyl-2-methylpent-1-en-4-yne"),
        ("CCCCCC(C=C)CCCCCC", "6-ethenyldodecane"),
        ("C=CCCCC(C=C)CCCCCC", "6-ethenyldodec-1-ene"),
    ],
)
def test_chain_parent_keeps_longest_chain_and_cites_the_rest_as_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(C)=Cc1ccccc1", "(2-methylprop-1-en-1-yl)benzene"),
        ("C=C(C)c1ccccc1", "(prop-1-en-2-yl)benzene"),
        ("C=CC=Cc1ccccc1", "(buta-1,3-dien-1-yl)benzene"),
        ("C(=Cc1ccccc1)c1ccccc1", "1,1'-(ethene-1,2-diyl)dibenzene"),
        ("C1CCCCC1=Cc1ccccc1", "(cyclohexylidenemethyl)benzene"),
        ("C1CCCCC1=C(C)c1ccccc1", "(1-cyclohexylideneethyl)benzene"),
        ("C1CCCC1=CC1CCCCC1", "(cyclopentylidenemethyl)cyclohexane"),
        ("C1CCCCCCC1C=Cc1ccccc1", "(2-phenylethen-1-yl)cyclooctane"),
        ("C1CCC1=C1CCCC1", "cyclobutylidenecyclopentane"),
        ("C1CCCCCCC1=Cc1ccccc1", "benzylidenecyclooctane"),
    ],
)
def test_aromatic_parent_and_ring_pairs_with_unsaturated_bridge(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC1CCCCC1=CC", "2-ethylidenecyclohexan-1-ol"),
        ("OC1CCC(=C)CC1", "4-methylidenecyclohexan-1-ol"),
        ("OC1CCCCC1C=C", "2-ethenylcyclohexan-1-ol"),
        ("OCC(=C)CC", "2-methylidenebutan-1-ol"),
        ("OCCC(=CC)CCC", "3-ethylidenehexan-1-ol"),
        ("OCCC(=C(C)C)CC", "3-ethyl-4-methylpent-3-en-1-ol"),
        ("O=C1CCCCC1=C", "2-methylidenecyclohexan-1-one"),
        ("CC(=O)CC(=C)CCC", "4-methylideneheptan-2-one"),
        ("NC1CCCCC1=C", "2-methylidenecyclohexan-1-amine"),
        ("NCCC(=CC)CCC", "3-ethylidenehexan-1-amine"),
        ("SC1CCCCC1=C", "2-methylidenecyclohexane-1-thiol"),
        ("OC(=O)C(=C)CC", "2-methylidenebutanoic acid"),
        ("OC(=O)C1CCCCC1=C", "2-methylidenecyclohexane-1-carboxylic acid"),
        ("O=CC(=C)CC", "2-methylidenebutanal"),
        ("O=CC1CCCCC1=C", "2-methylidenecyclohexane-1-carbaldehyde"),
        ("N#CC(=C)CC", "2-methylidenebutanenitrile"),
        ("N#CC1CCCCC1=C", "2-methylidenecyclohexane-1-carbonitrile"),
        ("NC(=O)C1CCCCC1=C", "2-methylidenecyclohexane-1-carboxamide"),
        ("COC(=O)C(=C)CC", "methyl 2-methylidenebutanoate"),
        ("COC(=O)C1CCCCC1=C", "methyl 2-methylidenecyclohexane-1-carboxylate"),
        ("OC(=O)c1ccccc1C=C", "2-ethenylbenzoic acid"),
        ("Oc1ccccc1C=C(C)C", "2-(2-methylprop-1-en-1-yl)phenol"),
        ("Nc1ccccc1C=C", "2-ethenylaniline"),
        ("C=CC(C)NCCCCC", "N-(but-3-en-2-yl)pentan-1-amine"),
        ("CC(=C)CNCC(=C)C", "2-methyl-N-(2-methylprop-2-en-1-yl)prop-2-en-1-amine"),
    ],
)
def test_ylidene_and_enyl_prefixes_alongside_principal_characteristic_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_prop_2_en_1_yl_ester_alcohol_part():
    assert smiles_to_iupac("CC(=O)OCC(=C)C").startswith("2-methylprop-2-en-1-yl ")


@pytest.mark.parametrize(
    "smiles,root,parent,expected",
    [
        ("[SiH]#CCC", 1, 0, "propylidyne"),
        ("[SiH]#CC", 1, 0, "ethylidyne"),
        ("[SiH]#CC=C", 1, 0, "prop-2-en-1-ylidyne"),
        ("[SiH2]=CCC", 1, 0, "propylidene"),
        ("[SiH2]=C(C)C", 1, 0, "propan-2-ylidene"),
        ("[SiH2]=C(C)C(C)C", 1, 0, "3-methylbutan-2-ylidene"),
        ("[SiH2]=C=C", 1, 0, "ethenylidene"),
        ("[SiH2]=C(Cl)C", 1, 0, "1-chloroethylidene"),
        ("[SiH2]=C(CCC)CC=C", 1, 0, "hept-1-en-4-ylidene"),
        ("[SiH2]=C(CC=C)CC#CC", 1, 0, "oct-1-en-6-yn-4-ylidene"),
        ("[SiH3]C=C", 1, 0, "ethenyl"),
        ("[SiH3]C#C", 1, 0, "ethynyl"),
        ("[SiH3]C(=C)CC", 1, 0, "but-1-en-2-yl"),
        ("[SiH3]C(C=C)C(CC)=CC", 1, 0, "4-ethylhexa-1,4-dien-3-yl"),
        ("[SiH3]C(C=C)CC(CC=C)C=CC", 1, 0, "5-(prop-2-en-1-yl)octa-1,6-dien-3-yl"),
        ("[SiH3]C(CCCCC)CC=C", 1, 0, "non-1-en-4-yl"),
        ("[SiH3]C(Cl)=CC", 1, 0, "1-chloroprop-1-en-1-yl"),
        ("[SiH3]C(C)(C)C", 1, 0, "tert-butyl"),
        ("[SiH2]=Cc1ccccc1", 1, 0, "benzylidene"),
        ("[SiH]#Cc1ccccc1", 1, 0, "benzylidyne"),
    ],
)
def test_name_branch_free_valence_kinds(smiles, root, parent, expected):
    mol = Chem.MolFromSmiles(smiles)
    name, _ = name_branch(adjacency(mol), root, parent, halogen_substituents(mol), mol=mol)
    assert name == expected


def test_name_branch_without_mol_still_reads_every_bond_as_single():
    mol = Chem.MolFromSmiles("[SiH3]CC=C")
    assert name_branch(adjacency(mol), 1, 0, {}) == ("propyl", False)


@pytest.mark.parametrize(
    "smiles",
    [
        "C/C=C1\\CCCCC1C",
        "C/C=C/C=C",
    ],
)
def test_specified_double_bond_geometry_on_a_prefix_is_not_silently_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_double_bond_geometry_on_a_prefix_is_cited_inside_the_prefix():
    assert smiles_to_iupac("C/C=C/C1CCCCC1") == "[(1E)-prop-1-en-1-yl]cyclohexane"
