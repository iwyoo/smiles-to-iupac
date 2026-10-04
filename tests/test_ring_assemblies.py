import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)Cc1ccc(cc1)-c1cc(-c3ccccc3)c(CC(=O)O)cc1",
            "2,2'-([11,21:23,31-terphenyl]-14,24-diyl)diethanoic acid",
        ),
        (
            "c1(-c3ccc4ccccc4c3CC(=O)O)cccc2c(CC(=O)O)cccc12",
            "2,2'-([1,2'-binaphthalene]-1',5-diyl)diethanoic acid",
        ),
        (
            "OC(=O)C=C2C=CC(c3ccc4ccccc4c3)c3ccccc23",
            "2-[[1,2'-binaphthalen]-4(1H)-ylidene]ethanoic acid",
        ),
    ],
)
def test_assembly_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "[11,21:24,31-terphenyl]-14-carboxylic acid"),
        ("Oc1ccc(nc1)-c1ccc(nc1)-c1ccccn1", "[12,22:25,32-terpyridin]-35-ol"),
    ],
)
def test_unbranched_assemblies_with_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_ring_chain_still_is_named():
    assert (
        smiles_to_iupac(
            "c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccccc1C"
        )
        == "12-methyl-1,7(1),2,3,4,5,6(1,4)-heptabenzenaheptaphane"
    )


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1=C1CCCCC1", "[1,1'-bi(cyclohexylidene)]-2-ol"),
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
    ],
)
def test_stereodescriptors_in_a_ring_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1cc2ccccc2cc1-c1ccc2ccccc2c1-c1ccc2ccccc2c1", "12,21:22,32-ternaphthalene"),
        ("c1ccc2[nH]c(cc2c1)-c1cc2ccccc2n1-c1cc2ccccc2[nH]1", "11H,31H-12,21:22,32-terindole"),
    ],
)
def test_assemblies_of_three_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCc2ccccc2C1C1CCc2ccccc2C1", "1,1',2,2',3,3',4,4'-octahydro-1,2'-binaphthalene"),
        ("C1CCC2CCCCC2C1-C1CCC2CCCCC2C1", "icosahydro-1,2'-binaphthalene"),
        (
            "C1Cc2ccccc2C1C1Cc2ccccc2C1C1Cc2ccccc2C1",
            "12,13,22,23,32,33-hexahydro-11H,21H,31H-11,22:21,32-terindene",
        ),
    ],
)
def test_hydro_and_indicated_hydrogen_of_fused_assembly_components(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1(=C2c3ccccc3-c3ccccc23)c2ccccc2-c2ccccc12", "9,9'-bifluorenylidene"),
        (
            "O=C1c2ccccc2C(=C2c3ccccc3C(=O)c3ccccc23)c2ccccc12",
            "[9,9'-bianthracenylidene]-10,10'-dione",
        ),
    ],
)
def test_double_bond_junction_between_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_4_chlorobiphenyl():
    assert smiles_to_iupac("Clc1ccc(cc1)-c1ccccc1") == "4-chloro-1,1'-biphenyl"


def test_bicyclopropane():
    assert smiles_to_iupac("C1CC1C1CC1") == "1,1'-bi(cyclopropane)"


def test_bipyrrole_carbon_attached_indicated_hydrogen():
    assert smiles_to_iupac("C1=CNC(=C1)C2=CC=CN2") == "1H,1'H-2,2'-bipyrrole"


def test_bipyrazole_nitrogen_attached():
    assert smiles_to_iupac("c1cnn(-n2cccn2)c1") == "1,1'-bipyrazole"


def test_branched_three_ring_assembly_is_named_on_its_longest_chain():
    # P-28.6: 1,3,5-triphenylbenzene.
    assert (
        smiles_to_iupac("c1ccc(-c2cc(-c3ccccc3)cc(-c4ccccc4)c2)cc1")
        == "25-phenyl-11,21:23,31-terphenyl"
    )


def test_halogen_substituent():
    assert smiles_to_iupac("C1C(Cl)C1C1CC1C1CC1") == "12-chloro-11,21:22,31-tercyclopropane"


def test_unsaturated_ring_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCCC1C1=CCCCC1C1=CCCCC1")


def test_mixed_thiophene_and_furan():
    assert smiles_to_iupac("c1ccc(-c2ccc(-c3ccco3)s2)s1") == "2-([2,2'-bithiophen]-5-yl)furan"


def test_terpyridine():
    assert smiles_to_iupac("C1=CC=NC(=C1)C2=NC(=CC=C2)C3=CC=CC=N3") == "12,22:26,32-terpyridine"


def test_terpyrrole_middle_ring_indicated_hydrogen():
    assert smiles_to_iupac("c1ccn(-c2cc(-n3cccc3)c[nH]2)c1") == "21H-11,22:24,31-terpyrrole"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC(=C2CCCC2)C1", "1,1'-bi(cyclopentylidene)"),  # CID 549125
    ],
)
def test_ring_assembly_ylidene_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituent__ring_assembly_ylidene():
    assert smiles_to_iupac("C1CCC(Cl)C1=C1CCCC1") == "2-chloro-1,1'-bi(cyclopentylidene)"


def test_bicyclic_ring_assembly_ylidene_bridge_carbon_attachment():
    # Attached at the one-carbon bridge (locant 7) instead of locant 2.
    assert smiles_to_iupac("C1CC2CCC1C2=C1C2CCC1CC2") == "7,7'-bi(bicyclo[2.2.1]heptanylidene)"


def test_bicyclic_ylidene_between_different_ring_systems_is_a_substituent():
    # Different ring systems are not "identical cyclic systems" (P-28.1): one is the parent.
    assert (
        smiles_to_iupac("C1CCC2C(=C3CC4CCC3C4)CCC2C1")
        == "1-(bicyclo[2.2.1]heptan-2-ylidene)octahydro-1H-indene"
    )


def test_bicyclic_ring_assembly_ylidene_halogen_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CC2CCC1C2=C1CC2CCC1C2")
