import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)Cc1ccc(cc1)-c1cc(-c3ccccc3)c(CC(=O)O)cc1",
            "2,2'-([1¹,2¹:2³,3¹-terphenyl]-1⁴,2⁴-diyl)diacetic acid",
        ),
        (
            "c1(-c3ccc4ccccc4c3CC(=O)O)cccc2c(CC(=O)O)cccc12",
            "2,2'-([1,2'-binaphthalene]-1',5-diyl)diacetic acid",
        ),
        (
            "OC(=O)C=C2C=CC(c3ccc4ccccc4c3)c3ccccc23",
            "[[1,2'-binaphthalen]-4(1H)-ylidene]acetic acid",
        ),
    ],
)
def test_assembly_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "[1¹,2¹:2⁴,3¹-terphenyl]-1⁴-carboxylic acid"),
        ("Oc1ccc(nc1)-c1ccc(nc1)-c1ccccn1", "[1²,2²:2⁵,3²-terpyridin]-3⁵-ol"),
    ],
)
def test_unbranched_assemblies_with_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_ring_chain_still_is_named():
    assert (
        smiles_to_iupac(
            "c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccc(cc1)-c1ccccc1C"
        )
        == "1²-methyl-1,7(1),2,3,4,5,6(1,4)-heptabenzenaheptaphane"
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
        ("c1cc2ccccc2cc1-c1ccc2ccccc2c1-c1ccc2ccccc2c1", "1²,2¹:2²,3²-ternaphthalene"),
        ("c1ccc2[nH]c(cc2c1)-c1cc2ccccc2n1-c1cc2ccccc2[nH]1", "1¹H,3¹H-1²,2¹:2²,3²-terindole"),
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
            "1²,1³,2²,2³,3²,3³-hexahydro-1¹H,2¹H,3¹H-1¹,2²:2¹,3²-terindene",
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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Clc1ccc(cc1)-c1ccccc1", "4-chloro-1,1'-biphenyl", id="4_chlorobiphenyl"),
        pytest.param("C1CC1C1CC1", "1,1'-bi(cyclopropane)", id="bicyclopropane"),
        pytest.param("C1=CNC(=C1)C2=CC=CN2", "1H,1'H-2,2'-bipyrrole", id="bipyrrole_carbon_attached_indicated_hydrogen"),
        pytest.param("c1cnn(-n2cccn2)c1", "1,1'-bipyrazole", id="bipyrazole_nitrogen_attached"),
        pytest.param("c1ccc(-c2cc(-c3ccccc3)cc(-c4ccccc4)c2)cc1", "2⁵-phenyl-1¹,2¹:2³,3¹-terphenyl", id="branched_three_ring_assembly_is_named_on_its_longest_chain"),
        pytest.param("C1C(Cl)C1C1CC1C1CC1", "1²-chloro-1¹,2¹:2²,3¹-tercyclopropane", id="halogen_substituent"),
        pytest.param("C1CCNC(C1)c1ccccn1", "1,2,3,4,5,6-hexahydro-2,2'-bipyridine", id="two_rings_mancude_and_saturated_take_hydro_prefixes"),
        pytest.param("C1CCOC1c1ccco1", "2,3,4,5-tetrahydro-2,2'-bifuran", id="five_membered_saturated_ring_leaves_the_chalcogen_unhydrogenated"),
        pytest.param("C1CCC(CC1)c1ccccc1", "cyclohexylbenzene", id="benzene_and_cyclohexane_pair_is_the_stated_exception"),
        pytest.param("C1CCC(CC1)c1ccc(cc1)c1ccccc1", "1¹,1²,1³,1⁴,1⁵,1⁶-hexahydro-1¹,2¹:2⁴,3¹-terphenyl", id="three_rings_mixing_benzene_and_cyclohexane_take_hydro_prefixes"),
        pytest.param("C1CCSC1c1ccc(s1)c1cccs1", "1²,1³,1⁴,1⁵-tetrahydro-1²,2²:2⁵,3²-terthiophene", id="three_rings_five_membered_saturated_ring"),
    ],
)
def test_4_chlorobiphenyl_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCCC1C1=CCCCC1C1=CCCCC1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccc(-c2ccc(-c3ccco3)s2)s1", "2-([2,2'-bithiophen]-5-yl)furan", id="mixed_thiophene_and_furan"),
        pytest.param("C1=CC=NC(=C1)C2=NC(=CC=C2)C3=CC=CC=N3", "1²,2²:2⁶,3²-terpyridine", id="terpyridine"),
        pytest.param("c1ccn(-c2cc(-n3cccc3)c[nH]2)c1", "21H-1¹,2²:2⁴,3¹-terpyrrole", id="terpyrrole_middle_ring_indicated_hydrogen"),
    ],
)
def test_mixed_thiophene_and_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC(=C2CCCC2)C1", "1,1'-bi(cyclopentylidene)"),  # CID 549125
    ],
)
def test_ring_assembly_ylidene_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1CCC(Cl)C1=C1CCCC1", "2-chloro-1,1'-bi(cyclopentylidene)", id="halogen_substituent__ring_assembly_ylidene"),
        pytest.param("C1CC2CCC1C2=C1C2CCC1CC2", "7,7'-bi(bicyclo[2.2.1]heptanylidene)", id="bicyclic_ring_assembly_ylidene_bridge_carbon_attachment"),
        pytest.param("C1CCC2C(=C3CC4CCC3C4)CCC2C1", "1-(bicyclo[2.2.1]heptan-2-ylidene)octahydro-1H-indene", id="bicyclic_ylidene_between_different_ring_systems_is_a_substituent"),
    ],
)
def test_halogen_substituent__and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bicyclic_ring_assembly_ylidene_halogen_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CC2CCC1C2=C1CC2CCC1C2")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1C1CCC(CC1)=C1CCCCC1", "4'-cyclohexylidene[1,1'-bi(cyclohexan)]-2-ol"),
        ("C1CCCCC1C1CCC(CC1)=C1CCCCC1", "4-cyclohexyl-1,1'-bi(cyclohexylidene)"),
        ("OC1CCCCC1=C1CCCCC1C1CCCCC1", "2'-cyclohexyl[1,1'-bi(cyclohexylidene)]-2-ol"),
    ],
)
def test_three_rings_with_a_double_bond_junction_name_a_two_ring_assembly_with_ring_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
