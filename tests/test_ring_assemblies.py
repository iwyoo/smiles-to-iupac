import pytest
from rdkit import Chem
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._ring_system_seniority import ring_seniority_key


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
            "([1,2'-binaphthalen]-4(1H)-ylidene)acetic acid",
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
        ("C1=CCCCC1C1CCCC=C1", "[1,1'-bi(cyclohexane)]-2,2'-diene"),
        ("C1=CCCCC1C1=CCCCC1", "[1,1'-bi(cyclohexane)]-1,2'-diene"),
        ("C1=CCCC(C1)C1CCCCC1", "[1,1'-bi(cyclohexane)]-3-ene"),
        ("C1=CC=CC(C1)C1C=CC=CC1", "[1,1'-bi(cyclohexane)]-2,2',4,4'-tetraene"),
        ("C1CC=CC1C1CC=CC1", "[1,1'-bi(cyclopentane)]-2,3'-diene"),
        ("C1CCC#CCCC1C1CCCCCCC1", "[1,1'-bi(cyclooctane)]-4-yne"),
        ("OC1CCC(C=C1)C1CCCCC1", "[1,1'-bi(cyclohexane)]-2-en-4-ol"),
        ("OC1CC=CCC1C1CC=CCC1", "[1,1'-bi(cyclohexane)]-3',4-dien-2-ol"),
        ("C1=CCCCC1C1=CCCCC1C", "6-methyl[1,1'-bi(cyclohexane)]-1,2'-diene"),
        ("ClC1=CCCCC1C1CCCCC1Cl", "2,2'-dichloro[1,1'-bi(cyclohexane)]-2-ene"),
    ],
)
def test_double_bond_junction_and_unsaturated_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1=CCCCC1C1=CCCCC1C1=CCCCC1", "[1¹,2¹:2²,3¹-tercyclohexane]-1¹,2²,3²-triene"),
        ("C1CC=CCC1C1CCC(CC1)C1CCCCC1", "[1¹,2¹:2⁴,3¹-tercyclohexane]-1³-ene"),
        ("C1CCCCC1C1CC#CCC1C1CCCCC1", "[1¹,2¹:2²,3¹-tercyclohexane]-2⁴-yne"),
        ("C1=CCCCC1C1CC#CCC1C1CCCCC1", "[1¹,2¹:2²,3¹-tercyclohexane]-1²-en-2⁴-yne"),
    ],
)
def test_ring_chain_unsaturation_is_cited_after_the_bracket(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C12(C=CC=C1C=CC=C2)C12C=CC=C1C=CC=C2", "3a,3'a-biindene"),
        ("N1C=CC(=CC=C1)C=1CC=NC=CC1", "1H,3'H-4,4'-biazepine"),
        ("N1(CC=CC=C1)C1=NC=CC=C1", "2H-1,2'-bipyridine"),
        ("C1(SCCCCCCCCCCCC1)C1COCCCCCCCCCCC1", "3'-oxa-2-thia-1,1'-bi(cyclotetradecane)"),
        ("C1(SCCCCCCCCCC1)=C1SCCCCCCCCCC1", "2,2'-dithia-1,1'-bi(cyclododecylidene)"),
        ("O1C(C=CC1=O)=C1C(OC=C1)=O", "2'H,5H-[2,3'-bifuranylidene]-2',5-dione"),
    ],
)
def test_assembly_indicated_hydrogen_primes_after_the_number_and_replacement_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CC2CCC1CC2C1CC2CCC1C=C2", "[2,2'-bi(bicyclo[2.2.2]octane)]-5,5'-diene"),
        ("C12C(C3CC4PCC3CC4)CC(NC1)CC2", "5-aza-5'-phospha-2,2'-bi(bicyclo[2.2.2]octane)"),
        ("ClC1CC2CCC1CC2C1CC2CCC1CC2Cl", "5,5'-dichloro-2,2'-bi(bicyclo[2.2.2]octane)"),
        ("C1C2C3CC(C3)C2C1C1C2C3CC(C3)C2C1", "3,3'-bi(tricyclo[4.1.1.0^2,5]octane)"),
        (
            "C12(C34CSC(C=C3)(C56CSC(C=C5)CC6)CC4)SCC(C=C1)CC2",
            "1²,2³,3³-trithia[1¹,2¹:2⁴,3¹-terbicyclo[2.2.2]octane]-1⁵,2⁵,3⁵-triene",
        ),
    ],
)
def test_assembly_of_von_baeyer_components(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def _senior_assembly(senior, junior):
    a, b = Chem.MolFromSmiles(senior), Chem.MolFromSmiles(junior)
    return ring_seniority_key(a, range(a.GetNumAtoms()), True) < ring_seniority_key(b, range(b.GetNumAtoms()), True)


@pytest.mark.parametrize(
    "senior,junior",
    [
        pytest.param("C1=CC=PC(=C1)C2=CC=CC=P2", "c1ccc(cc1)-c1ccccc1-c1ccccc1", id="heterocycle_before_carbocycle"),
        pytest.param("c1ccnc(c1)-c1ccccn1", "C1(C=CC=CO1)C1C=CC=CO1", id="nitrogen_before_oxygen"),
        pytest.param("C1(C=CC=CO1)C1C=CC=CO1", "C1(C=CC=CS1)C1C=CC=CS1", id="earlier_heteroatom"),
        pytest.param(
            "C1=CC=C2C(=C1)C=CC(=N2)C3=NC4=CC=CC=C4C=C3", "c1ccnc(c1)-c1cccc(n1)-c1ccccn1", id="more_rings"
        ),
        pytest.param("C1=CC=C(NC=C1)C1=CC=CC=CN1", "c1ccnc(c1)-c1ccccn1", id="more_atoms"),
        pytest.param("C1=CC(=NN=C1)C2=NN=CC=C2", "C1=CC=NC(=C1)C2=CN=CC=C2", id="more_heteroatoms"),
        pytest.param("O1C=COC(=C1)C1=COC=CO1", "O1C(=CSC=C1)C1=CSC=CO1", id="more_earlier_heteroatoms"),
    ],
)
def test_senior_ring_assembly(senior, junior):
    assert _senior_assembly(senior, junior)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "c1ccc(cc1)-c1ccccc1Cc1cccc2ccccc12",
            "2-[(naphthalen-1-yl)methyl]-1,1'-biphenyl",
            id="assembly_with_more_atoms_than_the_fused_system",
        ),
        pytest.param(
            "c1ccc(cc1)-c1ccccc1CC1CCCCC1", "2-(cyclohexylmethyl)-1,1'-biphenyl", id="assembly_beside_a_monocycle"
        ),
        pytest.param(
            "C1=CC=C(Cc2ccccc2-c2ccccc2)c2ccccc2C=C1",
            "5-[([1,1'-biphenyl]-2-yl)methyl]benzo[8]annulene",
            id="fused_system_before_an_assembly_of_equal_size",
        ),
    ],
)
def test_ring_assembly_beside_other_ring_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
