import pytest
from rdkit import Chem
from smiles_to_iupac import smiles_to_iupac
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
        pytest.param("OC1=CCCC=C1C1=CCCCC1", "[1,1'-bi(cyclohexane)]-1(6),1',2-trien-2-ol", id="ring_bond_closing_the_numbering_beside_a_suffix"),
        pytest.param("NC1=CCCC=C1C1CCCCC1", "[1,1'-bi(cyclohexane)]-1(6),2-dien-2-amine", id="compound_locant_beside_an_amine_suffix"),
        pytest.param("NC1=CCC=C1C1CCCC1", "[1,1'-bi(cyclopentane)]-1(5),2-dien-2-amine", id="compound_locant_in_five_membered_rings"),
        pytest.param(
            "OC1=CCCC=C1C1=CCCC=C1",
            "[1,1'-bi(cyclohexane)]-1(6),1',2,5'-tetraen-2-ol",
            id="second_ring_numbered_without_a_compound_locant",
        ),
        pytest.param("C1=CCCC=C1C1CCCCC1", "[1,1'-bi(cyclohexane)]-1,5-diene", id="numbering_without_a_compound_locant_when_free"),
    ],
)
def test_double_bond_junction_and_unsaturated_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O[C@@H]1CCCC[C@H]1C1CCCCC1", "(1S,2R)-[1,1'-bi(cyclohexan)]-2-ol"),
        pytest.param(
            "OC[C@H]1C[C@@H]1C1CC1", "[(1R,2S)-[1,1'-bi(cyclopropan)]-2-yl]methanol", id="substituent_of_a_chain_parent"
        ),
        pytest.param(
            "OC[C@@H]1C[C@H]1C1CC1", "[(1S,2R)-[1,1'-bi(cyclopropan)]-2-yl]methanol", id="substituent_of_the_enantiomer"
        ),
        pytest.param(
            "NC[C@H]1C[C@@H]1C1CC1", "1-[(1R,2S)-[1,1'-bi(cyclopropan)]-2-yl]methanamine", id="substituent_of_an_amine_chain"
        ),
        pytest.param(
            "OC[C@H]1CCCC[C@@H]1C1CCCCC1", "[(1R,2S)-[1,1'-bi(cyclohexan)]-2-yl]methanol", id="substituent_of_six_membered_rings"
        ),
        pytest.param(
            "OC(C)[C@H]1C[C@@H]1C1CC1", "1-[(1R,2S)-[1,1'-bi(cyclopropan)]-2-yl]ethan-1-ol", id="substituent_of_a_secondary_alcohol"
        ),
        pytest.param(
            "C[C@H]1C[C@@H]1[C@H]1C[C@@H]1CO",
            "[(1R,1'S,2S,2'S)-2'-methyl[1,1'-bi(cyclopropan)]-2-yl]methanol",
            id="substituent_with_stereo_in_both_rings",
        ),
        pytest.param(
            "OC[C@H]1C[C@@H]1[C@H]1C[C@@H]1CO",
            "[(1R,1'R,2S,2'S)-[1,1'-bi(cyclopropan)]-2,2'-diyl]dimethanol",
            id="multiplied_chain_parent_on_the_assembly",
        ),
        pytest.param(
            "OC[C@H]1C[C@@H]1[C@@H]1C[C@H]1CO",
            "[(1R,1'S,2S,2'R)-[1,1'-bi(cyclopropan)]-2,2'-diyl]dimethanol",
            id="multiplied_parent_numbering_gives_the_lowest_locant_to_r",
        ),
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
        ("C1C(=CCc2ccccc12)C1C=C2C=CC=CC2=CC1", "1,2',3',4-tetrahydro-2,2'-binaphthalene"),
        (
            "C1(C3=NC4=C(C=C3)C=CCC4)C=CC2=C(C1)N=C(C5=CC6=NC=CC=C6C=C5)C=C2",
            "1⁷,1⁸,2⁷,2⁸-tetrahydro-1²,2⁷:2²,3⁷-terquinoline",
        ),
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
        pytest.param("C1CNN(C1)N1CCCN1", "1,1'-bipyrazolidine", id="saturated_diazolidines_are_not_the_azoles"),
        pytest.param("c1ccc(-c2cc(-c3ccccc3)cc(-c4ccccc4)c2)cc1", "2⁵-phenyl-1¹,2¹:2³,3¹-terphenyl", id="branched_three_ring_assembly_is_named_on_its_longest_chain"),
        pytest.param("C1C(Cl)C1C1CC1C1CC1", "1²-chloro-1¹,2¹:2²,3¹-tercyclopropane", id="halogen_substituent"),
        pytest.param("C1CCNC(C1)c1ccccn1", "1,2,3,4,5,6-hexahydro-2,2'-bipyridine", id="two_rings_mancude_and_saturated_take_hydro_prefixes"),
        pytest.param("C1CCOC1c1ccco1", "2,3,4,5-tetrahydro-2,2'-bifuran", id="five_membered_saturated_ring_leaves_the_chalcogen_unhydrogenated"),
        pytest.param("c1ccc([C@@H]2CCCNC2)nc1", "(3'R)-1',2',3',4',5',6'-hexahydro-2,3'-bipyridine", id="stereocentre_on_the_primed_component"),
        pytest.param(
            "C[C@@H]1CCC[C@H](C)N1c1ccccn1", "(2R,6S)-2,6-dimethyl-1,2,3,4,5,6-hexahydro-1,2'-bipyridine", id="r_takes_the_lower_locant_on_a_tie"
        ),
        pytest.param("C1CCC(CC1)c1ccccc1", "cyclohexylbenzene", id="benzene_and_cyclohexane_pair_is_the_stated_exception"),
        pytest.param("C1CCC(CC1)c1ccc(cc1)c1ccccc1", "1¹,1²,1³,1⁴,1⁵,1⁶-hexahydro-1¹,2¹:2⁴,3¹-terphenyl", id="three_rings_mixing_benzene_and_cyclohexane_take_hydro_prefixes"),
        pytest.param("C1CCSC1c1ccc(s1)c1cccs1", "1²,1³,1⁴,1⁵-tetrahydro-1²,2²:2⁵,3²-terthiophene", id="three_rings_five_membered_saturated_ring"),
        pytest.param("C1(=CC=CCC1)C1=CC=CC=C1", "2,3-dihydro-1,1'-biphenyl", id="two_rings_benzene_and_cyclohexadiene_take_hydro_prefixes"),
        pytest.param("C1=CCCCC1c1ccccc1", "1,2,3,4-tetrahydro-1,1'-biphenyl", id="benzene_and_cyclohexene_take_hydro_prefixes"),
        pytest.param("C1CC=CC=C1C1=CC=C(C=C1)C1=CC=CC=C1", "1²,1³-dihydro-1¹,2¹:2⁴,3¹-terphenyl", id="three_rings_benzene_and_cyclohexadiene_take_hydro_prefixes"),
    ],
)
def test_4_chlorobiphenyl_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1cc(-c2cc[nH]n2)n[nH]1", "1H,1'H-3,3'-bipyrazole", id="nh_azoles_cite_the_indicated_hydrogen"),
        pytest.param("Cc1c[nH]nc1-c1n[nH]cc1C", "4,4'-dimethyl-1H,1'H-3,3'-bipyrazole", id="carbon_substituent_precedes_indicated_hydrogen"),
        pytest.param("Cn1ccnc1-c1ncc[nH]1", "1-methyl-1H,1'H-2,2'-biimidazole", id="substituted_nitrogen_keeps_its_indicated_hydrogen"),
        pytest.param("c1c[nH]c(n1)-c1cnc[nH]1", "1H,3'H-2,4'-biimidazole", id="lowest_attachment_locant_before_indicated_hydrogen"),
        pytest.param("Cn1nccc1-c1ccnn1C", "2,2'-dimethyl-2H,2'H-3,3'-bipyrazole", id="substituted_nitrogen_beside_the_junction"),
    ],
)
def test_azole_assembly_cites_the_indicated_hydrogen_of_both_rings(smiles, expected):
    # P-28.2.1, P-28.2.3, P-58.2.1
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("O=c1[nH]c(-c2ccc(=O)[nH]c2)ccc1", "[2,3'-bipyridine]-6,6'(1H,1'H)-dione", id="added_hydrogen_after_the_suffix_locants"),
        pytest.param("S=c1[nH]nc(-c2n[nH]c(=S)s2)s1", "[2,2'-bi(1,3,4-thiadiazole)]-5,5'(4H,4'H)-dithione", id="sulfanylidene_suffixes"),
        pytest.param("O=c1ccn(-c2ccc(=O)[nH]c2)cc1", "[1,3'-bipyridine]-4,6'(1H,1'H)-dione", id="junction_on_the_ring_nitrogen_keeps_its_added_hydrogen"),
        pytest.param("O=c1ccc(-c2ccc(=O)oc2)co1", "6H,6'H-[3,3'-bipyran]-6,6'-dione", id="indicated_hydrogen_accommodates_the_suffix"),
        pytest.param(
            "O=c1cc(-c2cc(=O)[nH][nH]2)[nH][nH]1",
            "1,1',2,2'-tetrahydro-5H,5'H-[3,3'-bipyrazole]-5,5'-dione",
            id="hydro_prefixes_beside_indicated_hydrogen",
        ),
    ],
)
def test_assembly_of_ring_ketones_cites_the_hydrogen_of_each_ring(smiles, expected):
    # P-28.2.3, P-31.1.4.3.4, P-58.2.1
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1=CCCCC1C1=CCCCC1C1=CCCCC1", "[1¹,2¹:2²,3¹-tercyclohexane]-1¹,2²,3²-triene"),
        ("C1CC=CCC1C1CCC(CC1)C1CCCCC1", "[1¹,2¹:2⁴,3¹-tercyclohexan]-1³-ene"),
        ("C1CCCCC1=C1CCC(CC1)C1CCCCC1", "[1¹,2¹:2⁴,3¹-tercyclohexan]-1¹(2¹)-ene"),
        ("C1CCCCC1=C1CCC(CC1)=C1CCCCC1", "[1¹,2¹:2⁴,3¹-tercyclohexane]-1¹(2¹),2⁴(3¹)-diene"),
        ("C1CCCCC1=C1C=CC(CC1)C1C=CCCC1", "[1¹,2¹:2⁴,3¹-tercyclohexane]-1¹(2¹),2²,3²-triene"),
        ("C1CCCCC1C1CC#CCC1C1CCCCC1", "[1¹,2¹:2²,3¹-tercyclohexan]-2⁴-yne"),
        ("C1=CCCCC1C1CC#CCC1C1CCCCC1", "[1¹,2¹:2²,3¹-tercyclohexan]-1²-en-2⁴-yne"),
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
        pytest.param("ClC1CC2CCC1C2=C1CC2CCC1C2", "7-(bicyclo[2.2.1]heptan-2-ylidene)-2-chlorobicyclo[2.2.1]heptane", id="identical_ring_systems_the_one_with_more_substituents_is_the_parent"),
    ],
)
def test_halogen_substituent__and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1C1CCC(CC1)=C1CCCCC1", "4'-cyclohexylidene[1,1'-bi(cyclohexan)]-2-ol"),
        ("OC1CCCCC1=C1CCCCC1C1CCCCC1", "2'-cyclohexyl[1,1'-bi(cyclohexylidene)]-2-ol"),
    ],
)
def test_three_rings_with_a_double_bond_junction_name_a_two_ring_assembly_with_ring_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_double_bond_between_assembly_components_is_cited_with_both_locants():
    assert smiles_to_iupac("CC[C@H]1CC/C(CC1)=C1\\CC[C@@H](CC)CC1") == "[1(1')E,4s,4's]-4,4'-diethyl-1,1'-bi(cyclohexylidene)"


def test_assembly_holding_the_principal_group_is_the_parent_beside_a_fused_ring_substituent():
    assert smiles_to_iupac("OC1CCCCN1C1CCCCN1c1ccc2ccccc2c1") == "1'-(naphthalen-2-yl)[1,2'-bipiperidin]-2-ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C12(C=CC=C1C=CC=C2)C12C=CC=C1C=CC=C2", "3a,3'a-biindene"),
        ("N1C=CC(=CC=C1)C=1CC=NC=CC1", "1H,3'H-4,4'-biazepine"),
        ("N1(CC=CC=C1)C1=NC=CC=C1", "2H-1,2'-bipyridine"),
        ("C1(SCCCCCCCCCCCC1)C1COCCCCCCCCCCC1", "3'-oxa-2-thia-1,1'-bi(cyclotetradecane)"),
        ("C1(SCCCCCCCCCC1)=C1SCCCCCCCCCC1", "2,2'-dithia-1,1'-bi(cyclododecylidene)"),
        ("O1C(C=CC1=O)=C1C(OC=C1)=O", "2'H,5H-[2,3'-bifuranylidene]-2',5-dione"),
        # P-31.1.4.2.4, P-58.2: only the indicated hydrogens the mancude ring needs; the other saturated position is added hydrogen
        ("O=C1C=CC(=C2CC=COC2=O)OC1", "2'H,4'H-[2,3'-bipyranylidene]-2',5(6H)-dione"),
        ("O=C1CC=CC(=C2C(=O)C=COC2)O1", "2'H,4'H-[2,3'-bipyranylidene]-4',6(5H)-dione"),
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
        ("OC1(CC2CCC1C2)C1CC2CCC1C2", "[2,2'-bi(bicyclo[2.2.1]heptan)]-2-ol"),
        ("OC1CC2CCC1CC2C1CC2CCC1CN2", "5-aza[2,2'-bi(bicyclo[2.2.2]octan)]-5'-ol"),
        ("OC(=O)CC1CC2CCC1C2C1CC2CCC1C2", "{[2,7'-bi(bicyclo[2.2.1]heptan)]-2'-yl}acetic acid"),
        (
            "CC(=O)NC1CC2CCC1CC2C1CC2CCC1CC2NC(C)=O",
            "N,N'-{[2,2'-bi(bicyclo[2.2.2]octane)]-5,5'-diyl}diacetamide",
        ),
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "C(c1cc(nnc1)-c1cccnn1)c1ccnc(c1)-c1cccnc1",
            "5-[([2,3'-bipyridin]-4-yl)methyl]-3,3'-bipyridazine",
            id="assembly_of_rings_without_a_ring_specification",
        ),
        pytest.param(
            "C(c1cc2ccccc2nc1-c1ccc2ccccc2n1)c1ccc(nc1)-c1ncc(cc1)-c1ccccn1",
            "3-[([1²,2²:2⁵,3²-terpyridin]-1⁵-yl)methyl]-2,2'-biquinoline",
            id="fused_pair_before_a_chain_of_three_rings",
        ),
        pytest.param(
            "C(c1ccc(nc1)-c1ccccn1)c1ccc(cc1)-c1ccccc1",
            "5-[([1,1'-biphenyl]-4-yl)methyl]-2,2'-bipyridine",
            id="heterocyclic_assembly_before_a_carbocyclic_one",
        ),
    ],
)
def test_senior_ring_assembly_is_the_parent_beside_another_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param(
            "NC(=O)c1cccc(CN2CCCCC2)c1-c1ccccc1",
            "6-[(piperidin-1-yl)methyl][1,1'-biphenyl]-2-carboxamide",
            id="piperidine_substituent_beside_a_biphenyl_amide",
        ),
        pytest.param(
            "CN1CCCCC1CCc1ccccc1-c1ccc(F)cc1O",
            "4-fluoro-2'-[2-(1-methylpiperidin-2-yl)ethyl][1,1'-biphenyl]-2-ol",
            id="piperidine_substituent_beside_a_biphenyl_ol",
        ),
        pytest.param(
            "CCCCCC1CCC(c2cc(C=O)ccc2-c2ccccc2)CC1",
            "2-(4-pentylcyclohexyl)[1,1'-biphenyl]-4-carbaldehyde",
            id="cyclohexyl_substituent_beside_a_biphenyl_aldehyde",
        ),
    ],
)
def test_assembly_with_a_principal_group_keeps_other_rings_as_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[Sb]1(CCCCC1)[Sb]1CCCCC1", "1,1'-bistibinane", id="stibinane_joined_through_the_group_15_atoms"),
        pytest.param("P1(CCCCC1)P1CCCCC1", "1,1'-biphosphinane", id="phosphinane_joined_through_the_ring_atoms"),
        pytest.param("N1(CCCCC1)N1CCCCC1", "1,1'-bipiperidine", id="retained_ring_name_is_kept"),
        pytest.param("[SiH]1(CCCC1)[SiH]1CCCC1", "1,1'-bisilolane", id="silicon_ring_assembly"),
        pytest.param("C1COC(O1)C1OCCO1", "2,2'-bi(1,3-dioxolane)", id="ring_name_with_heteroatom_locants_is_enclosed"),
        pytest.param("C1CCNCC1C1CCCCN1", "2,3'-bipiperidine", id="carbon_junction_with_lowest_locant_on_the_unprimed_ring"),
        pytest.param(
            "Cl[Sb]1CCCCC1[Sb]1CCCCC1C",
            "1'-chloro-2-methyl-1,2'-bistibinane",
            id="substituents_follow_the_junction_locants",
        ),
        pytest.param("[Sb]1(CCCCC1)[Sb]1CCCC1", "1-(stibolan-1-yl)stibinane", id="different_rings_are_not_an_assembly"),
    ],
)
def test_assembly_of_saturated_heteromonocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-44.1.2.2, P-44.2.1, P-29.2: two different rings joined through their nitrogens are parent and substituent, the senior ring being the parent
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1CCN(CC1)N1CCCC1", "1-(pyrrolidin-1-yl)piperidine", id="more_skeletal_atoms"),
        pytest.param("C1CCN(CC1)N1CCOCC1", "4-(piperidin-1-yl)morpholine", id="oxygen_heteroatom_beside_nitrogen"),
        pytest.param("C1CCN(CC1)N1CCNCC1", "1-(piperidin-1-yl)piperazine", id="more_heteroatoms"),
        pytest.param("CN1CCN(CC1)N1CCCCC1", "1-methyl-4-(piperidin-1-yl)piperazine", id="substituent_on_the_parent"),
        pytest.param("OC1CCN(CC1)N1CCCC1", "1-(pyrrolidin-1-yl)piperidin-4-ol", id="suffix_on_the_saturated_ring"),
        pytest.param("c1ccn(N2CCOCC2)c1", "4-(1H-pyrrol-1-yl)morpholine", id="azole_substituent"),
        pytest.param("c1ccn(N2CCCCC2)c1", "1-(1H-pyrrol-1-yl)piperidine", id="azole_beside_a_six_membered_ring"),
        pytest.param("c1ccn(N2CCNCC2)c1", "1-(1H-pyrrol-1-yl)piperazine", id="azole_beside_a_diazine"),
        pytest.param("c1cn(N2CCOCC2)cn1", "4-(1H-imidazol-1-yl)morpholine", id="imidazole_substituent"),
        pytest.param("Cc1ccn(N2CCOCC2)c1", "4-(3-methyl-1H-pyrrol-1-yl)morpholine", id="substituted_azole_group"),
        pytest.param("c1ccc2c(c1)ccn2N1CCOCC1", "1-(morpholin-4-yl)-1H-indole", id="two_rings_outrank_one"),
        pytest.param("c1cnn(N2CCCC2)c1", "1-(pyrrolidin-1-yl)-1H-pyrazole", id="more_heteroatoms_in_an_azole"),
        pytest.param("OC(=O)c1ccn(N2CCCC2)c1", "1-(pyrrolidin-1-yl)-1H-pyrrole-3-carboxylic acid", id="acid_on_the_azole"),
    ],
)
def test_different_rings_joined_nitrogen_to_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-62.2.2, P-62.2.4.1.2, P-28.2.1: the substituents on the nitrogen of an amine of a ring assembly take the locant N, with the number and primes of the group for a diamine
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CNc1ccc(-c2ccccc2)cc1", "N-methyl[1,1'-biphenyl]-4-amine", id="secondary_amine"),
        pytest.param("CN(C)c1ccc(-c2ccccc2)cc1", "N,N-dimethyl[1,1'-biphenyl]-4-amine", id="tertiary_amine"),
        pytest.param("CCNc1ccc(-c2ccccc2)cc1", "N-ethyl[1,1'-biphenyl]-4-amine", id="ethyl_on_nitrogen"),
        pytest.param("CNc1cccc(-c2ccccc2)c1", "N-methyl[1,1'-biphenyl]-3-amine", id="amine_at_position_3"),
        pytest.param("CNc1ccccc1-c1ccccc1", "N-methyl[1,1'-biphenyl]-2-amine", id="amine_at_position_2"),
        pytest.param("CNc1ccc(-c2ccccc2)cc1C", "N,3-dimethyl[1,1'-biphenyl]-4-amine", id="ring_substituent_beside_the_n_prefix"),
        pytest.param("CNc1ccc(-c2ccccc2)cc1Cl", "3-chloro-N-methyl[1,1'-biphenyl]-4-amine", id="alphanumerical_order_with_n"),
        pytest.param("c1ccc(Nc2ccc(-c3ccccc3)cc2)cc1", "N-phenyl[1,1'-biphenyl]-4-amine", id="assembly_outranks_a_single_ring"),
        pytest.param("c1ccc(Nc2ccccc2-c2ccccc2)cc1", "N-phenyl[1,1'-biphenyl]-2-amine", id="phenyl_on_the_nitrogen_at_position_2"),
        pytest.param("CN(c1ccccc1)c1ccc(-c2ccccc2)cc1", "N-methyl-N-phenyl[1,1'-biphenyl]-4-amine", id="two_groups_on_the_nitrogen"),
        pytest.param("CNc1ccc(-c2ccncc2)cc1", "N-methyl-4-(pyridin-4-yl)aniline", id="different_rings_are_not_an_assembly"),
        pytest.param("CNC1(C2CC2)CC1", "N-methyl[1,1'-bi(cyclopropan)]-1-amine", id="cycloalkane_assembly"),
        pytest.param("FCCNC1(C2CC2)CC1", "N-(2-fluoroethyl)[1,1'-bi(cyclopropan)]-1-amine", id="substituted_group_on_the_nitrogen"),
        pytest.param("CNC1(C2CCCC2)CCCC1", "N-methyl[1,1'-bi(cyclopentan)]-1-amine", id="five_membered_rings"),
        pytest.param("CNc1ccc(-c2ccc(NC)cc2)cc1", "N4,N4'-dimethyl[1,1'-biphenyl]-4,4'-diamine", id="diamine_one_group_on_each_nitrogen"),
        pytest.param(
            "CN(C)c1ccc(-c2ccc(N(C)C)cc2)cc1", "N4,N4,N4',N4'-tetramethyl[1,1'-biphenyl]-4,4'-diamine", id="diamine_two_groups_on_each_nitrogen"
        ),
        pytest.param("CCNc1ccc(-c2ccc(NC)cc2)cc1", "N4'-ethyl-N4-methyl[1,1'-biphenyl]-4,4'-diamine", id="diamine_different_groups"),
        pytest.param("CNc1ccc(-c2ccc(N)cc2)cc1", "N4'-methyl[1,1'-biphenyl]-4,4'-diamine", id="diamine_with_one_primary_group"),
        pytest.param("CNc1ccc(-c2ccc(O)cc2)cc1", "4'-(methylamino)[1,1'-biphenyl]-4-ol", id="alcohol_outranks_the_amine"),
        pytest.param("CNc1ccc2ccccc2c1-c1ccc2ccccc2c1", "N-methyl[1,2'-binaphthalen]-2-amine", id="fused_system_assembly"),
        pytest.param(
            "CNc1ccc2ccccc2c1-c1c(NC)ccc2ccccc12", "N2,N2'-dimethyl[1,1'-binaphthalene]-2,2'-diamine", id="fused_system_assembly_diamine"
        ),
    ],
)
def test_substituted_amines_of_a_ring_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC1CCCCN1C1CCCCN1", "[1,2'-bipiperidin]-2-ol", id="suffix_on_the_ring_joined_through_its_heteroatom"),
        pytest.param("OC1COC(CO1)C1COCCO1", "[2,2'-bi(1,4-dioxan)]-5-ol", id="ring_name_with_locants_is_enclosed_before_a_suffix"),
        pytest.param("O=C1CCCCN1C1CCCCN1", "[1,2'-bipiperidin]-2-one", id="oxo_on_one_ring_only"),
        pytest.param("OC(=O)CC1CCN(CC1)C1CCNCC1", "([1,4'-bipiperidin]-4-yl)acetic acid", id="assembly_group"),
        pytest.param(
            "CC(=O)NC1CCNCC1C1CCNCC1NC(C)=O", "N,N'-([3,4'-bipiperidin]-3',4-diyl)diacetamide", id="assembly_linker"
        ),
        pytest.param("OC1CCN(CC1)C1CCN(CC1)C1CCNCC1", "[1¹,2⁴:2¹,3⁴-terpiperidin]-1⁴-ol", id="chain_of_three_rings"),
        pytest.param(
            "OC1COCCCCCCCCCCC1C1CNCCCCCCCCCCC1",
            "4'-oxa-3-aza[1,1'-bi(cyclotetradecan)]-2'-ol",
            id="replacement_prefixes_of_a_large_ring_assembly_before_the_bracket",
        ),
        pytest.param(
            "OC(=O)CC1COCCCCCCCCCCC1C1COCCCCCCCCCCC1",
            "(3,4'-dioxa[1,1'-bi(cyclotetradecan)]-2'-yl)acetic acid",
            id="replacement_prefixes_inside_an_assembly_group",
        ),
    ],
)
def test_assembly_of_saturated_heteromonocycles_as_parent_group_and_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
