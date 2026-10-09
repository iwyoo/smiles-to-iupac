import pytest
from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.slow
def test_morphine():
    smiles = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"
    assert (
        smiles_to_iupac(smiles)
        == "17-methyl-7,8-didehydro-4,5α-epoxymorphinan-3,6α-diol"
    )


@pytest.mark.slow
def test_codeine():
    smiles = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@H]3[C@H](C=C4)O"
    assert (
        smiles_to_iupac(smiles)
        == "3-methoxy-17-methyl-7,8-didehydro-4,5α-epoxymorphinan-6α-ol"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("B", "borane"),
        ("ClB", "chloroborane"),
    ],
)
def test_smiles_to_iupac_simple_borane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituted_alkyl_chain():
    assert smiles_to_iupac("ClCCB") == "(2-chloroethyl)borane"


# P-67.1.2.6.2 amides and hydrazides of the boron acids and silicic acid
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("NB(N)N", "boranetriamine", id="boron_triamide"),
        pytest.param("N[Si](N)(N)N", "silanetetramine", id="tetra_elides_its_final_a"),
        pytest.param("N[SiH2]N", "silanediamine", id="silicon_diamide"),
        pytest.param("NNB(NN)NN", "1,1′,1′′-boranetriyltrihydrazine", id="boron_trihydrazide"),
        pytest.param("NN[BH2]", "boranylhydrazine", id="boron_monohydrazide"),
        pytest.param("[Si](O)(O)(O)OC#N", "cyanatosilicic acid", id="cyanato_prefix_on_silicic_acid"),
        pytest.param("[Si](O)(O)(OC#N)OC#N", "dicyanatosilicic acid", id="multiplied_cyanato_prefix"),
    ],
)
def test_amides_and_cyanates_of_boron_and_silicon_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("BB", "diborane(4)", id="parent_hydride_cites_its_hydrogen_count"),
        pytest.param("CBB", "1-methyldiborane(4)", id="substituted_diborane"),
        pytest.param("CB(C)B(C)B(C)C", "1,1,2,3,3-pentamethyltriborane(5)", id="triborane_chain"),
        pytest.param("ClBB(Cl)Cl", "1,1,2-trichlorodiborane(4)", id="lowest_locants_for_prefixes"),
        pytest.param("B=BB", "triborene(5)", id="double_bond_keeps_the_saturated_hydrogen_count"),
        pytest.param("BB(B)B", "2-boranyltriborane(5)", id="branched_skeleton_cites_the_boron_branch"),
        pytest.param("OC(=O)CBBB", "[triboran(5)-1-yl]acetic acid", id="senior_group_makes_the_chain_a_prefix"),
        pytest.param("CB(C)OB(C)C", "tetramethyldiboroxane", id="full_substitution_omits_locants"),
    ],
)
def test_borane_chains(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[N+](C)(C)[B-](Cl)(Cl)Cl", "N,N-dimethylmethanamine—trichloroborane (1/1)", id="charge_separated_adduct"),
        pytest.param("N->B", "ammonia—borane (1/1)", id="dative_bond_adduct"),
        pytest.param("CC[S+](CC)[BH2-]C", "(ethylsulfanyl)ethane—methylborane (1/1)", id="sulfur_donor"),
        pytest.param("CO[NH2+][BH3-]", "O-methylhydroxylamine(N—B)borane (1/1)", id="attachment_cited_for_several_donors"),
        pytest.param("[BH3-][NH2+]C(=O)Nc1ccccc1", "N-phenylurea(N'—B)borane (1/1)", id="urea_unsubstituted_nitrogen"),
        pytest.param("[BH3-][NH+](C)C(=O)Nc1ccccc1", "N-methyl-N'-phenylurea(N—B)borane (1/1)", id="urea_substituted_nitrogen"),
        pytest.param("CN(C)C(=[NH+][GaH3-])N(C)C", "N,N,N',N'-tetramethylguanidine(N''—Ga)gallane (1/1)", id="guanidine_imino_nitrogen"),
        pytest.param("C[NH+]([BH3-])C(=N)NC", "N,N'-dimethylguanidine(N—B)borane (1/1)", id="guanidine_amino_nitrogen"),
    ],
)
def test_lewis_adducts_of_boranes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C=CB", "ethenylborane", id="unsaturated_substituent"),
        pytest.param("CB(c1ccccc1)C", "dimethyl(phenyl)borane", id="mixed_alkyl_and_phenyl_substituents"),
        pytest.param("Cc1ccccc1B", "(2-methylphenyl)borane", id="substituted_phenyl"),
        pytest.param("C1CCCCC1B", "cyclohexylborane", id="non_aromatic_ring"),
        pytest.param("CB(c1ccc(Cl)cc1)", "(4-chlorophenyl)(methyl)borane", id="halogenated_phenyl_mixed_with_alkyl"),
    ],
)
def test_unsaturated_substituent_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_asymmetric_substituents():
    assert smiles_to_iupac("CCB(C)O") == "ethyl(methyl)borinic acid"
    assert smiles_to_iupac("CC(C)B(C)O") == "methyl(propan-2-yl)borinic acid"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OB(C)CCB(C)O", id="second_borinic_acid_group"),
        pytest.param("NCB(C)O", id="unrecognized_heteroatom"),
    ],
)
def test_rejects_second_borinic_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_benzene_ring():
    assert smiles_to_iupac("c1ccccc1B(O)O") == "phenylboronic acid"


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OB(O)CCB(O)O", id="second_boronic_acid_group"),
        pytest.param("NCB(O)O", id="unrecognized_heteroatom__boronic_acid"),
    ],
)
def test_rejects_second_boronic_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Re]([OH2])Cl", "aquachlorido(methyl)rhenium"),
    ],
)
def test_coordination_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C[Ti](Cl)(Cl)Cl.[Na+]",
        "C[Hg]c1ccc(S(=O)(=O)O)cc1",
    ],
)
def test_coordination_out_of_scope_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[SbH2][SnH3]", "stannylstibane"),
    ],
)
def test_metal_pair_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metal_pair_linked_by_carbon_chain():
    assert (
        smiles_to_iupac("c1ccccc1[Bi](c1ccccc1)CCC[Pb](CC)(CC)CC")
        == "diphenyl[3-(triethylplumbyl)propyl]bismuthane"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1cc(C(=O)O)ccc1[Hg]C", "(4-carboxyphenyl)(methyl)mercury"),
        ("CC[Sn](CC)(CC)c1ccc(cc1)[Ge](C)(C)C", "trimethyl[4-(triethylstannyl)phenyl]germane"),
    ],
)
def test_substituted_aryl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "[Mo+](C#[O+])(C#[O+])(C#[O+])C.C1=C[CH-]C=C1",
            "tricarbonyl(η5-cyclopenta-2,4-dien-1-yl)(methyl)molybdenum",
        ),
    ],
)
def test_dinuclear_and_hapto_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1[Hg][Sb](c1ccccc1)c1ccccc1", "diphenylstibanyl(phenyl)mercury"),
    ],
)
def test_class1_metal_with_class2_metal_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Ti](=O)(Cl)Cl", "dichlorido(methyl)oxidotitanium"),
        ("[Ti](O)(O)(O)O", "tetrahydroxidotitanium"),
        ("N#C[Au]C#N", "dicyanidogold"),
        ("[Ti](OC)(OC)(OC)OC", "tetramethanolatotitanium"),
        ("CC(=O)[Pt](C)(P(CC)(CC)CC)P(CC)(CC)CC", "acetyl(methyl)bis(triethylphosphane)platinum"),
        ("[Fe](N=O)(C#[O+])C", "carbonyl(methyl)nitrosyliron"),
    ],
)
def test_anionic_and_acyl_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "Cl[Pd-]1(Cl)Cl[Pd-](Cl)(Cl)Cl1.[Na+].[Na+]",
            "disodium di-\u03bc-chlorido-tetrachlorido-1\u03ba2Cl,2\u03ba2Cl-dipalladate(2-)",
        ),
    ],
)
def test_mu_bridged_dinuclear(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Pb]([SnH3])([SnH3])([SnH3])[SnH3]", "plumbanetetrayltetrakis(stannane)"),
        ("[GeH3][GeH2][GeH2][GeH3]", "tetragermane"),
        ("CC[Sn](Cl)(CC)[Sn](CC)(CC)Cl", "1,2-dichloro-1,1,2,2-tetraethyldistannane"),
    ],
)
def test_metal_chain_and_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[ReH2]c1ccc2ccccc2c1", "dihydrido(naphthalen-2-yl)rhenium"),
    ],
)
def test_naphthyl_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CP(C)CCP(CC)CC", "[2-(dimethylphosphanyl)ethyl]di(ethyl)phosphane"),
    ],
)
def test_polyphosphane_and_chelates(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


K = "\u03ba"
MU = "\u03bc"
DASH = "\u2014"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "O=[C]1[Fe]23([C]#[O+])([C]#[O+])([C]#[O+])[C](=O)[Fe]12([C]#[O+])([C]#[O+])([C]#[O+])[Fe]3([C]#[O+])([C]#[O+])([C]#[O+])[C]#[O+]",
            f"di-{MU}-carbonyl-decacarbonyl-1{K}3C,2{K}3C,3{K}4C-triangulo-triiron(3 Fe{DASH}Fe)",
        ),
        (
            "[O+]#[C][Ru]1([C]#[O+])([C]#[O+])([C]#[O+])[Ru]([C]#[O+])([C]#[O+])([C]#[O+])([C]#[O+])[Ru]1([C]#[O+])([C]#[O+])([C]#[O+])[C]#[O+]",
            f"dodecacarbonyl-1{K}4C,2{K}4C,3{K}4C-triangulo-triruthenium(3 Ru{DASH}Ru)",
        ),
        (
            "[O+]#[C][Mn]([C]#[O+])([C]#[O+])([C]#[O+])([C]#[O+])[Re]([C]#[O+])([C]#[O+])([C]#[O+])([C]#[O+])[C]#[O+]",
            f"decacarbonyl-1{K}5C,2{K}5C-rheniummanganese(Re{DASH}Mn)",
        ),
        (
            "C[O]1->[Cu]([Cl])[O](C)->[Cu]1[Cl]",
            f"dichlorido-1{K}Cl,2{K}Cl-di-{MU}-methanolato-dicopper",
        ),
        (
            "C(#[O+])[Ru]1([H]->[Ru]1(C#[O+])(C#[O+]))(C#[O+])C#[O+]",
            f"pentacarbonyl-1{K}3C,2{K}2C-{MU}-hydrido-diruthenium(Ru{DASH}Ru)",
        ),
    ],
)
def test_polynuclear_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCC(=C(C)C)CC1", "(propan-2-ylidene)cyclohexane"),
    ],
)
def test_exocyclic_ylidene_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Clc1ccc(cc1)[GaH2]", "(4-chlorophenyl)gallane"),
    ],
)
def test_group13_hydride_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[Tl](C)O", "dimethylthallanol", id="group13_hydroxy_suffix"),
        pytest.param("C[Al](C)[O-].[Na+]", "sodium dimethylalumanolate", id="group13_olate"),
        pytest.param("[Ga](SSCC)(SSCC)SSCC", "tris(ethyldisulfanyl)gallane", id="disulfanyl_on_group13"),
        pytest.param("C[Sn](C)(C)SC", "trimethyl(methylsulfanyl)stannane", id="sulfanyl_on_group14"),
        pytest.param(
            "CCCCCCCCCCCCCCCCCC(=O)O[Al](OC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC",
            "alumanetriyl tri(octadecanoate)",
            id="pseudoester_of_one_hydride_atom",
        ),
        pytest.param(
            "CCCC[Sn](CCCC)(OC(=O)CCC(=O)OC)OC(=O)CCC(=O)OC",
            "dimethyl dibutylstannanediyl dibutanedioate",
            id="pseudoester_beside_alkyl_esters",
        ),
        pytest.param("OC(=O)c1ccc(B(S)O)cc1", "4-(thioborono)benzoic acid", id="chalcogen_analogue_of_borono"),
        pytest.param("Oc1ccc(BOB)cc1", "4-diboroxanylphenol", id="diboroxanyl_prefix"),
        pytest.param(
            "OC1CCCCC1C[BH]C[SiH2]C[SiH2]C[SiH2]C",
            "2-(2,4,6-trisila-8-boranonan-9-yl)cyclohexan-1-ol",
            id="skeletal_replacement_prefix_with_boron",
        ),
        pytest.param("C[Al]1CCCC2CCCCC12", "1-methyldecahydro-1-benzaluminine", id="hydro_prefixes_on_an_aluminium_ring"),
        pytest.param(
            "CCOB(c1ccccc1)OSC", "O-ethyl OS-methyl phenylborono(thioperoxoate)", id="boron_acid_ester_with_a_peroxy_group"
        ),
    ],
)
def test_group13_and_boron_acid_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bare_metal_hydride_name():
    assert smiles_to_iupac("[AlH3]") == "alumane"


def test_two_metal_atoms_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Al](C)C.C[Ga](C)C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Cl[Sn](c1ccccc1)(c1ccccc1)c1ccccc1", "chlorotri(phenyl)stannane"),
    ],
)
def test_group14_hydride_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_group13_hydride_unaffected():
    assert smiles_to_iupac("CC[Al](CC)CC") == "triethylalumane"


def test_unsaturated_substituent_across_groups():
    assert smiles_to_iupac("C=C[As](C=C)C=C") == "tri(ethenyl)arsane"
    assert smiles_to_iupac("CCCC[Sn](CCCC)(CCCC)C=C") == "tributyl(ethenyl)stannane"


def test_multiple_bond_directly_to_metal_is_a_ylidene_prefix():
    assert smiles_to_iupac("C=[Sb]CC") == "ethyl(methylidene)stibane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Ca]C", "dimethylcalcium"),
    ],
)
def test_group2_organometallic_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Mg](Cl)c1ccccc1", "phenylmagnesium chloride"),
    ],
)
def test_grignard_rmx_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_inorganic_dihalide_of_group2_metal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cl[Mg]Cl")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Mg]CC", "ethyl(methyl)magnesium"),
    ],
)
def test_group2_additive_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_metal_atoms_raises__group1_2_organometallic():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Li]C.[Na]CC")


def test_group1_metal_with_substituted_alkyl():
    assert smiles_to_iupac("[Li]CO") == "hydroxymethyllithium"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Li][CH3][Li]", "μ-methyl-dilithium"),
        ("Cl[Hg]c1ccc([Hg]Cl)s1", "dichlorido-1κCl,2κCl-μ-thiophene-2,5-diyl-dimercury"),
        (
            "C[As](C)c1ccncc1[Hg]c1ccc([Hg]O)s1",
            "[4-(dimethylarsanyl)pyridin-3-yl]-1κC-hydroxido-2κO-μ-thiophene-2,5-diyl-dimercury",
        ),
        ("C[Sb](C)c1ccc(cc1)[Hg]C", "[4-(dimethylstibanyl)phenyl](methyl)mercury"),
        ("CB(C)c1ccc(cc1)[Hg]C", "[4-(dimethylboranyl)phenyl](methyl)mercury"),
    ],
)
def test_carbon_bridged_dilithium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


H = "\u03b7"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[cH]12->[Cr]3456<-[cH]1[cH]->3[cH]->4[cH]->5[cH]->62", f"({H}6-benzene)chromium"),
        (
            "[O+]#[C][Mn]1234([C]#[O+])([C]#[O+])[CH]5[CH]->1=[CH]->2[CH]->3=[CH]->45",
            f"tricarbonyl({H}5-cyclopenta-2,4-dien-1-yl)manganese",
        ),
        ("[CH2]1[CH]2=[CH2]->[Cr]<-21", f"({H}3-allyl)chromium"),
        (
            "[CH2]1[CH]2=[CH2]->[Cr]<-213456([CH2][CH]->3=[CH2]->4)[CH2][CH]->5=[CH2]->6",
            f"tris({H}3-allyl)chromium",
        ),
        (
            "C1C[CH]2->[Rh]34<-[CH]1=[CH]->3CC[CH]->4=2",
            f"[(1,2,5,6-{H})-cycloocta-1,5-diene]rhodium",
        ),
        (
            "[O+]#[C][Fe]123([C]#[O+])([C]#[O+])<-[CH]4=[CH]->1C1CC4[CH]->2=[CH]->31",
            f"[(2,3,5,6-{H})-bicyclo[2.2.1]hepta-2,5-diene]tricarbonyliron",
        ),
        (
            "[O+]#[C][Mo+]123456([C]#[O+])([C]#[O+])[CH]7[CH]->1=[CH]->2[CH]->3=[CH]->4[CH]->5=[CH]->67",
            f"tricarbonyl({H}7-cyclohepta-2,4,6-trien-1-yl)molybdenum(1+)",
        ),
        (
            "[O+]#[C][Fe]123([C]#[O+])([C]#[O+])<-[CH2]=[CH]->1[CH]->2=[CH2]->3",
            f"({H}4-buta-1,3-diene)tricarbonyliron",
        ),
    ],
)
def test_hapto_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "NCC[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->41)[cH]->5[cH]->62",
            f"tricarbonyl[2-({H}6-phenyl)ethan-1-amine]chromium",
        ),
        (
            "C[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->4[c]->51CCN)[cH]->62",
            f"tricarbonyl[2-(2-methyl-{H}6-phenyl)ethan-1-amine]chromium",
        ),
        (
            "CC(N(C)C)[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->4[c]->51P(c1ccccc1)c1ccccc1)[cH]->62",
            f"tricarbonyl{{1-[2-(diphenylphosphanyl)-{H}6-phenyl]-N,N-dimethylethan-1-amine}}chromium",
        ),
        (
            "C[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->4[cH]->51)[cH]->62",
            f"tricarbonyl({H}6-methylbenzene)chromium",
        ),
        (
            "c1ccc([B-](c2ccccc2)(c2ccccc2)[c]23->[Rh+]456789%10(<-[CH]%11=[CH]->4CC[CH]->5=[CH]->6CC%11)"
            "<-[cH]([cH]->7[cH]->82)[cH]->9[cH]->%103)cc1",
            f"[(1,2,5,6-{H})-cycloocta-1,5-diene][triphenyl({H}6-phenyl)borato]rhodium",
        ),
    ],
)
def test_hapto_ring_inside_larger_ligand(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "[O+]#[C][Mo]123456([C]#[O+])([CH]7C=CC=C[CH]->1=[CH]->27)[CH]1[CH]->3=[CH]->4[CH]->5=[CH]->61",
            "dicarbonyl[(1–3-η)-cyclohepta-2,4,6-trien-1-yl](η5-cyclopenta-2,4-dien-1-yl)molybdenum",
        ),
        (
            "C[C]12->[Fe]345([C]#[O+])([C]#[O+])[CH]1[c]->31cccc[c]->41[CH]->5=2",
            "dicarbonyl[2-methyl-(1–3,3a,7a-η)-1H-inden-1-yl]iron",
        ),
        (
            "COC(=O)[C]12[Mn]345([C]#[O+])([C]#[O+])([C]#[O+])<-[CH](=[CH]->31)[CH]->4=[CH]->52",
            "tricarbonyl[1-(methoxycarbonyl)-η5-cyclopenta-2,4-dien-1-yl]manganese",
        ),
        (
            "C[Si](C)(C)[C]12[Mn]345([C]#[O+])([C]#[O+])([C]#[O+])<-[CH](=[CH]->31)[CH]->4=[CH]->52",
            "tricarbonyl[1-(trimethylsilyl)-η5-cyclopenta-2,4-dien-1-yl]manganese",
        ),
        (
            "C[C]12->[Fe]345([C]#[O+])([C]#[O+])<-[CH](=[CH]->3[N]->41)[CH]->5=2",
            "dicarbonyl(2-methyl-η5-1H-pyrrol-1-yl)iron",
        ),
        (
            "[O+]#[C][Fe]1234([C]#[O+])<-[CH]5=[CH]->1[N]->2[CH]->3=[CH]->45",
            "dicarbonyl(η5-1H-pyrrol-1-yl)iron",
        ),
        (
            "[O+]#[C][Fe]1234([C]#[O+])([C]#[O+])[CH2][CH]->1=[CH]->2[CH]->3=[CH2]->4",
            "tricarbonyl(η5-penta-2,4-dien-1-yl)iron",
        ),
        ("[Cl][Pt]1([Cl])<-[CH]#[CH]->1", "(η2-acetylene)dichloridoplatinum"),
        (
            "[O+]#[C][Fe]123([C]#[O+])([C]#[O+])<-[cH]4cc[c]56->[Fe]1789([C]#[O+])([C]#[O+])<-[cH]([cH]->75)[cH]->8[c]->9-6[cH]->2[cH]->34",
            "{μ-[2(1–3,3a,8a-η):1(4–6-η)]azulene}-(pentacarbonyl-1κ3C,2κ2C)diiron(Fe—Fe)",
        ),
        (
            "[O+]#[C][Fe]12([C]#[O+])([C]#[O+])<-[CH2]=[CH]->1[CH]1=[CH2]->[Fe]<-12([C]#[O+])([C]#[O+])[C]#[O+]",
            "{μ-[2(1,2-η):1(3,4-η)]buta-1,3-diene}-(hexacarbonyl-1κ3C,2κ3C)diiron(Fe—Fe)",
        ),
    ],
)
def test_extended_hapto_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metallacycle_with_unsaturated_ligand():
    assert smiles_to_iupac("[Ni]1(C=C)CCCC1") == "1-ethenyl-1-nickelacyclopentane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1[Pt](Cl)(Cl)O[Pt]1(Cl)Cl", "1,1,3,3-tetrachloro-2-oxa-1,3-diplatinacyclobutane"),
    ],
)
def test_metallacycle_with_several_metals(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "[Si]1(Cl)(Cl)[Fe](C#[O+])(C#[O+])(C#[O+])(C#[O+])CCC1",
            "2,2,2,2-tetracarbonyl-1,1-dichloro-1-sila-2-ferracyclopentane",
        ),
        ("C[Ir]1(C#[O+])=CC(C)=CC(C)=C1", "1-carbonyl-1,3,5-trimethyl-1-iridabenzene"),
    ],
)
def test_hetero_and_metallabenzene_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C1CC2C(C1)C[Ti]2.C1=C[CH]C=C1.C1=C[CH]C=C1",
            "6,6-di(η5-cyclopenta-2,4-dien-1-yl)-6-titanabicyclo[3.2.0]heptane",
        ),
    ],
)
def test_bicyclic_metallacycle_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CC2C=CC1[Pt]2", "7-platinabicyclo[2.2.1]hepta-2,5-diene"),
    ],
)
def test_unsaturated_bicyclic_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C1c2ccccc2[Pt]2(P(C)(C)CP2(C)C)c2ccccc12",
            "9,9-[methylenebis(dimethylphosphane)]-10H-9-platinaanthracene", marks=pytest.mark.slow),
        (
            "[Pt]1(P(c2ccccc2)(c2ccccc2)c2ccccc2)(P(c2ccccc2)(c2ccccc2)c2ccccc2)C2CC1C(OC)CCC2OC",
            "2,5-dimethoxy-7,7-bis(triphenylphosphane)-7-platinabicyclo[4.1.1]octane",
        ),
        (
            "CC1=C(C)[Pt](P(c2ccccc2)(c2ccccc2)c2ccccc2)(P(c2ccccc2)(c2ccccc2)c2ccccc2)C(C)=C1C",
            "2,3,4,5-tetramethyl-1,1-bis(triphenylphosphane)-1-platinacyclopenta-2,4-diene",
        ),
    ],
)
def test_metallacycles_with_ring_and_chelating_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1c2ccccc2[Pt](Cl)(Cl)C=C1", "1,1-dichloro-4H-1-platinanaphthalene"),
        ("C1=Cc2ccccc2C[Pt]1(Cl)Cl", "2,2-dichloro-1H-2-platinanaphthalene"),
        ("C1c2cccc(C)c2[Pt](Cl)(Cl)C=C1", "1,1-dichloro-8-methyl-4H-1-platinanaphthalene"),
        ("C1c2ccccc2-c2ccccc2[Pt]1(Cl)Cl", "9,9-dichloro-10H-9-platinaphenanthrene"),
        ("C1=Cc2ccccc2[Pt]1(Cl)Cl", "1,1-dichloro-1-platinaindene"),
        (
            "C1c2ccccc2-c2ccccc2[Pt]1(Cl)Cl".replace("C1c2ccccc2-c2ccccc2", "C1c2ccccc2-c2ccccc2"),
            "9,9-dichloro-10H-9-platinaphenanthrene",
        ),
    ],
)
def test_fused_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC2(C1)C[Ti]C2", "2-titanaspiro[3.3]heptane"),
        ("C1C2CC3C1C[Ti]3C2", "3-titanatricyclo[3.2.1.0^3,6]octane"),
    ],
)
def test_spiro_and_polycyclic_metallacycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[C@H](C)C1CC[Pt](Cl)(Cl)C1", "3-[(2S)-butan-2-yl]-1,1-dichloro-1-platinacyclopentane"),
    ],
)
def test_metallacycle_substituent_stereocentre(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_stereocentre_on_a_metal_ligand_is_not_dropped():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](Cl)[Pt]1(Cl)CCCC1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "OC(=O)CC1CC2CC[Pt](Cl)(Cl)C12",
            "{2,2-dichloro-2-platinabicyclo[3.2.0]heptan-7-yl}acetic acid",
        ),
        (
            "OC(=O)CC1C[Pt](Cl)(Cl)c2ccccc12",
            "(1,1-dichloro-2,3-dihydro-1-platinainden-3-yl)acetic acid",
        ),
    ],
)
def test_polycyclic_metallacycle_cited_as_a_substituent_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC2CC[Pt]1(Cl)(Cl)C2c1ccccc1", "1,1-dichloro-7-phenyl-1-platinabicyclo[2.2.1]heptane"),
        (
            "C1CC2CC[Pt]1(Cl)(Cl)C2N(C)C",
            "1,1-dichloro-N,N-dimethyl-1-platinabicyclo[2.2.1]heptan-7-amine",
        ),
        ("C1CC(O)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentan-3-ol"),
        ("C1CC(C(=O)O)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentane-3-carboxylic acid"),
        ("C1CC(=O)C[Pt]1(Cl)Cl", "1,1-dichloro-1-platinacyclopentan-3-one"),
    ],
)
def test_metallacycle_ring_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC[Pt+]1(Cl)Cl", "1,1-dichloro-1-platinacyclobutan-1-ium"),
    ],
)
def test_metallacycle_ylidene_hydrido_ionic_and_metals(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@H]1CC[Pt]1(Cl)Cl", "(2S)-1,1-dichloro-2-methyl-1-platinacyclobutane"),
    ],
)
def test_metallacycle_ring_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC[Pt]1=C", "1-methylidene-1-platinacyclobutane"),
    ],
)
def test_heteroatom_rings_and_ylidene_on_ring_metal(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=C[CH]C=C1.C1=C[CH]C=C1.[Os]", "osmocene"),  # neutral biradical form, CID 102601604
    ],
)
def test_metallocene_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


PRIME = "′"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Os+2].OCC[c-]1cccc1.[cH-]1cccc1", "2-(osmocen-1-yl)ethan-1-ol"),
        (
            "[Fe+2].CC(=O)[c-]1cccc1.CC(=O)[c-]1cccc1",
            f"1,1{PRIME}-(ferrocene-1,1{PRIME}-diyl)di(ethan-1-one)",
        ),
        ("[Fe+2].C[c-]1cccc1.[cH-]1cccc1", "1-methylferrocene"),
        ("[Fe+2].OC(=O)[c-]1cccc1.[cH-]1cccc1", "ferrocene-1-carboxylic acid"),
        ("[Fe+2].OC(=O)[c-]1cccc1.OC(=O)[c-]1cccc1", f"ferrocene-1,1{PRIME}-dicarboxylic acid"),
        ("[Fe+2].N[c-]1cccc1.[cH-]1cccc1", "ferrocen-1-amine"),
        ("[Fe+2].O[c-]1cccc1.[cH-]1cccc1", "ferrocen-1-ol"),
        ("[Fe+2].N#C[c-]1cccc1.[cH-]1cccc1", "ferrocene-1-carbonitrile"),
        ("[Fe+2].OC(=O)[c-]1cccc1.CC(=O)[c-]1cccc1", f"1{PRIME}-acetylferrocene-1-carboxylic acid"),
        ("[Fe+2].OC(=O)CC[c-]1cccc1.[cH-]1cccc1", "3-(ferrocen-1-yl)propanoic acid"),
        ("[V+2].CN(C)C(C)[c-]1cccc1.[cH-]1cccc1", "N,N-dimethyl-1-(vanadocen-1-yl)ethan-1-amine"),
    ],
)
def test_substituted_metallocenes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ansa_metallocene_cited_as_divalent_group():
    name = smiles_to_iupac("OC(=O)CC([c-]1cccc1)CC[c-]1cccc1.[Fe+2]")
    assert name == f"3,5-(ferrocene-1,1{PRIME}-diyl)pentanoic acid"


def test_two_ruthenocenes_joined_by_a_chain():
    smiles = (
        "C[C]12->[Ru]3456789(<-[CH](=[CH]->3[CH2]->41)[CH]->5=2)<-[CH]1=[CH]->6[CH2]->7[C]->8(CC[C]23->"
        "[Ru]45678%10%11(<-[CH]%12=[CH]->4[CH2]->5[C]->6(C)=[CH]->7%12)<-[CH](=[CH]->8[CH2]->%102)[CH]->%11=3)=[CH]->91"
    )
    assert (
        smiles_to_iupac(smiles)
        == f"1,1{PRIME}{PRIME}-(ethane-1,2-diyl)bis(1{PRIME}-methylruthenocene)"
    )


def test_benzoferrocene():
    smiles = "c1cc[c]23->[Fe]456789%10(<-[CH]%11=[CH]->4[CH2]->5[CH]->6=[CH]->7%11)<-[CH](=[CH]->8[c]->92c1)[CH2]->%103"
    assert smiles_to_iupac(smiles) == "benzoferrocene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC1[C]23->[Fe]456789%10(<-[CH]%11=[CH]->4[CH2]->5[C]->6(C[C]45->[Fe]6%12%13%14%15%16%17(<-[CH](=[CH]->6[CH2]->%124)[CH]->%13=5)<-[CH]4=[CH]->%14[CH]->%15(C)[C]->%161=[CH]->%174)=[CH]->7%11)<-[CH](=[CH]->8[C]->9=2C)[CH2]->%103",
            "1^2,2,3^2-trimethyl-1,3(1,1′)-diferrocenacyclotetraphane",
        ),
    ],
)
def test_substituted_phanes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_single_bridged_ferrocene_warns_that_no_pin_is_defined():
    from smiles_to_iupac import NonPreferredNameWarning

    with pytest.warns(NonPreferredNameWarning, match="phane"):
        name = smiles_to_iupac("C(CC[c-]1cccc1)[c-]1cccc1.[Fe+2]")
    assert name == f"1,3-(ferrocene-1,1{PRIME}-diyl)propane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[SiH4]", "silane"),
    ],
)
def test_smiles_to_iupac_silane_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_silane():
    assert smiles_to_iupac("[SiH2]1[SiH2][SiH2][SiH2][SiH2]1") == "pentasilolane"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CCCCC1=C(CC)C2CCCCC2", "1,1'-(propan-1-yl-1-ylidene)dicyclohexane"),
        ("C1CCCCC1=CC=C2CCCCC2", "1,1'-(ethane-1,2-diylidene)dicyclohexane"),
        (
            "OC(=O)C1CCC(CC1)=CCC2CCC(C(=O)O)CC2",
            "4,4'-(ethan-1-yl-2-ylidene)di(cyclohexane-1-carboxylic acid)",
        ),
    ],
)
def test_yl_ylidene_linker_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C=C2CCCCC2CC(=O)O", "2,2'-(cyclohexan-1-yl-2-ylidene)diacetic acid"),
        (
            "OC(=O)CC1CCC(CC1)C1CCC(CC1)=CC(=O)O",
            "2,2'-([1,1'-bi(cyclohexan)]-4-yl-4'-ylidene)diacetic acid",
        ),
    ],
)
def test_ring_yl_ylidene_linker_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_trimethylammonium_methylide():
    assert smiles_to_iupac("[CH2-][N+](C)(C)C") == "(N,N-dimethylmethanaminiumyl)methanide"


def test_tertiary_ammonium_ylide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2-][N+](C)C")


def test_mixed_substituent_sulfonium_ylide():
    assert smiles_to_iupac("[CH2-][S+](C)CC") == "[ethyl(methyl)sulfaniumyl]methanide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[CH-][N+](C)(C)C", "1-(N,N-dimethylmethanaminiumyl)ethan-1-ide"),
        ("C[CH-][P+](C)(C)C", "1-(trimethylphosphaniumyl)ethan-1-ide"),
        ("C[CH-][O+](C)C", "1-(dimethyloxidaniumyl)ethan-1-ide"),
        ("C[CH-][S+](C)C", "1-(dimethylsulfaniumyl)ethan-1-ide"),
    ],
)
def test_ylide_with_a_substituted_anion_carbon(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C#CC(CC)C(=C)C", "3-ethyl-2-methylpent-1-en-4-yne"),
    ],
)
def test_chain_parent_keeps_longest_chain_and_cites_the_rest_as_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCCCC1=C(C)c1ccccc1", "(1-cyclohexylideneethyl)benzene"),
    ],
)
def test_aromatic_parent_and_ring_pairs_with_unsaturated_bridge(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CC1CCCCC1=C", "2-methylidenecyclohexane-1-carbaldehyde"),
        ("N#CC1CCCCC1=C", "2-methylidenecyclohexane-1-carbonitrile"),
    ],
)
def test_ylidene_and_enyl_prefixes_alongside_principal_characteristic_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "OC(=O)CC1CCCC\\C1=C/C",
    ],
)
def test_specified_double_bond_geometry_on_a_prefix_is_not_silently_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_double_bond_geometry_on_a_prefix_is_cited_inside_the_prefix():
    assert smiles_to_iupac("C/C=C/C1CCCCC1") == "[(1E)-prop-1-en-1-yl]cyclohexane"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)c1ccc(cc1)C=NOCON=Cc1ccc(cc1)C(=O)O",
            "4,4'-(3,5-dioxa-2,6-diazahepta-1,6-diene-1,7-diyl)dibenzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)C(Cl)=NCCN=C(Cl)c1ccc(cc1)C(=O)O",
            "4,4'-{ethane-1,2-diylbis[azanylylidene(chloromethanylylidene)]}dibenzoic acid",
        ),
    ],
)
def test_ylylidene_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_imine_linker_without_senior_unit_group_is_not_multiplicative():
    with pytest.raises(Exception):
        smiles_to_iupac("c1ccccc1C=NCCN=Cc1ccccc1")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C=[SiH2]", "methylidenesilane"),
        ("C=[GeH2]", "methylidenegermane"),
        ("CC=[SnH2]", "ethylidenestannane"),
        ("CC#[SiH]", "ethylidynesilane"),
        ("CCC#[SiH]", "propylidynesilane"),
        ("C=[Si](C)C", "dimethyl(methylidene)silane"),
        ("C#[Si][Si]#CC", "ethylidyne(methylidyne)disilane"),
        ("C=[SiH][SiH3]", "methylidenedisilane"),
        ("C[Si](C)(C)[Si](C)(C)C", "hexamethyldisilane"),
    ],
)
def test_ylidene_and_ylidyne_prefixes_on_group_14_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccccc1[SiH2][SiH2]c1ccccc1", "1,2-diphenyldisilane"),
        ("C1CCCCC1[SiH2]C1CCCCC1", "dicyclohexylsilane"),
        ("[SiH3][Si]([SiH3])([SiH3])[SiH3]", "2,2-disilyltrisilane"),
        ("[SiH3][SiH]([SiH3])[SiH2][SiH3]", "2-silyltetrasilane"),
        ("[SiH3][SiH2][SiH]([SiH3])[SiH2][SiH]([SiH3])[SiH2][SiH3]", "3,5-disilylheptasilane"),
    ],
)
def test_group_14_hydride_outranks_two_rings_and_branched_silicon_chain_is_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1[SiH]=C", "2-(methylidenesilyl)cyclohexan-1-ol"),
        ("OC1CCCCC1P=C", "2-(methylidenephosphanyl)cyclohexan-1-ol"),
        ("OC1CCCCC1[SiH]=[SiH2]", "2-(silylidenesilyl)cyclohexan-1-ol"),
        ("OC1CCCCC1[SiH2][SiH]=C", "2-(2-methylidenedisilanyl)cyclohexan-1-ol"),
        ("OC1CCCCC1N=C", "2-(methylideneamino)cyclohexan-1-ol"),
    ],
)
def test_ylidene_substituent_on_heteroatom_attaching_a_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[Si](C)(C)O", "trimethylsilanol"),
        ("C[Si](C)(O)O", "dimethylsilanediol"),
        ("C[Si](C)(C)OC", "methoxytri(methyl)silane"),
        ("Cl[Si](Cl)(Cl)C", "trichloro(methyl)silane"),
        ("C[Si](C)(C)N", "trimethylsilanamine"),
        ("C[Si](N)(N)N", "methylsilanetriamine"),
        ("CN[Si](O)(O)O", "(methylamino)silanetriol"),
        ("C[Si](C)(C)C([Si](C)(C)C)[SnH](O)C([Si](C)(C)C)[Si](C)(C)C", "bis[bis(trimethylsilyl)methyl]stannanol"),
        ("C[Si](C)(C)N(C)C", "N,N,1,1,1-pentamethylsilanamine"),
        ("C[Si](C)(C)N[Si](C)(C)C", "1,1,1-trimethyl-N-(trimethylsilyl)silanamine"),
        ("C[Sn](C)(C)O", "trimethylstannanol"),
    ],
)
def test_group_14_hydride_with_hydroxy_amino_alkoxy_and_halogen_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles", ["[Pb]=[Pb]", "[Sn]=[Sn]"])
def test_non_single_bond_in_a_metal_chain_raises(smiles):
    with pytest.raises(NotImplementedError):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O1[SiH2]O[SiH2]O[SiH2]1", "1,3,5,2,4,6-trioxatrisilinane"),
        ("C1N=PN=PN1", "1,6-dihydro-1,3,5,2,4-triazadiphosphinine"),
        ("[CH3][Bi]1[O][Bi]([CH3])[O]1", "2,4-dimethyl-1,3,2,4-dioxadibismetane"),
    ],
)
def test_alternating_heteroatom_ring_cites_every_heteroatom_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O", "water"),
        ("N", "ammonia"),
        ("Cl", "hydrogen chloride"),
        ("[SeH2]", "hydrogen selenide"),
    ],
)
def test_common_hydride_names_carry_no_pin(smiles, expected):
    with pytest.warns(NonPreferredNameWarning, match="P-21.1.1.2"):
        assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("N1NNNN1", "pentazolidine"),
        ("N1NNNNN1", "hexazinane"),
        ("N1N=NN=N1", "1H-pentazole"),
        ("O1OCOCO1", "1,2,3,5-tetroxane"),
        ("[SiH2]1[SiH2][SiH2][SiH2][SiH2][SiH2][SiH2][SiH2][SiH2][SiH2][SiH2][SiH2]1", "dodecasilacyclododecane"),
        ("S1SSSSSSSS1", "nonathionane"),
        ("O1OOOOOOOO1", "nonoxonane"),
        ("S1SSSSSSSSSSSS1", "tridecathiacyclotridecane"),
    ],
)
def test_homogeneous_heteromonocycle_elides_multiplier_vowel_and_omits_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[GeH2]S", "methylgermanethiol"),
        ("[SiH3][SeH]", "silaneselenol"),
        ("[SiH2](S)S", "silanedithiol"),
        ("[SiH2](O)S", "sulfanylsilanol"),
    ],
)
def test_group_14_hydrides_take_chalcogenol_suffixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C[P@](CCC)C1=CC=CC=C1", "(R)-methyl(phenyl)(propyl)phosphane", id="stereogenic_phosphorus_parent"),
        pytest.param("C[P@@](CCC)C1=CC=CC=C1", "(S)-methyl(phenyl)(propyl)phosphane", id="stereogenic_phosphorus_parent_enantiomer"),
        pytest.param("C[Si@@H](O)CCC", "(R)-methyl(propyl)silanol", id="stereogenic_silicon_with_hydroxy_suffix"),
        pytest.param("C[Si@H](O)CCC", "(S)-methyl(propyl)silanol", id="stereogenic_silicon_enantiomer"),
    ],
)
def test_stereodescriptor_of_a_mononuclear_hydride_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[TlH3]", "thallane", id="thallane"),
        pytest.param("FB(F)[Si](C)(C)C", "(difluoroboranyl)tri(methyl)silane", id="boranyl_on_the_senior_silane"),
        pytest.param("B(O)(O)C1=CC(=C(C(=O)O)C=C1)[N+](=O)[O-]", "4-borono-2-nitrobenzoic acid", id="borono_prefix_beside_a_carboxylic_acid"),
        pytest.param("B(O)(O)C1=CC=C(C(=O)O)C=C1", "4-boronobenzoic acid", id="borono_on_benzoic_acid"),
        pytest.param("CBNBC", "1-methyl-N-(methylboranyl)boranamine", id="boranamine_with_a_boranyl_group"),
        pytest.param("[GeH3][Se][GeH2][Se][GeH3]", "trigermaselenane", id="alternating_germanium_selenium_chain"),
    ],
)
def test_noncarbon_hydride_parents_and_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("B(SC)(SC)SC", "trimethyl borotrithioate"),
        ("B(OC)(O)SCC", "S-ethyl O-methyl hydrogen borothioate"),
        ("CB(O)[O-]", "hydrogen methylboronate"),
        ("CB([O-])[O-].[Na+].[Na+]", "disodium methylboronate"),
        ("CB(C)[O-].[Na+]", "sodium dimethylborinate"),
        ("C(CCC)BNBNBCCCC", "N,N'-bis(butylboranyl)boranediamine"),
        ("C1(=CC=CC=C1)B1OC2=C(N=CN2)O1", "2-phenyl-2H,4H-[1,3,2]dioxaborolo[4,5-d]imidazole"),
    ],
)
def test_boron_acid_esters_anions_boranediamines_and_fused_indicated_hydrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O[Si](O)(O)S", "thiosilicic acid"),
        ("CCO[Si](SC)(SC)SC", "O-ethyl S,S,S-trimethyl trithiosilicate"),
        ("OB(O)B(O)O", "hypodiboric acid"),
        ("C1(O[Pb]2(OC(C3=C1C=CC=C3)=O)OC(C3=C(C(O2)=O)C=CC=C3)=O)=O", "3,3′-spirobi[[2,4,3]benzodioxaplumbepine]-1,1′,5,5′-tetrone"),
        ("C1O[Sn](C)(C)OCc2ccccc12", "3,3-dimethyl-1,5-dihydro-3H-2,4,3-benzodioxastannepine"),
    ],
)
def test_silicic_boric_acids_and_group_14_fused_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1(=CC=CC=C1)B(NCCNB(C1=CC=CC=C1)C1=CC=CC=C1)C1=CC=CC=C1", "N1,N2-bis(diphenylboranyl)ethane-1,2-diamine"),
        ("NCCNB(C)C", "N1-(dimethylboranyl)ethane-1,2-diamine"),
    ],
)
def test_diamines_whose_nitrogens_carry_boranyl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CCC(CO)B(C(CC)CO)OC", "methyl bis(1-hydroxybutan-2-yl)borinate"),
        ("OCCB(CCO)O", "bis(2-hydroxyethyl)borinic acid"),
        ("OCCB(O)O", "(2-hydroxyethyl)boronic acid"),
    ],
)
def test_boron_acids_outrank_hydroxy_groups_of_their_organyl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C[Si](C)(C)O[Si](C)(C)C", "hexamethyldisiloxane", id="one_arrangement_needs_no_locants"),
        pytest.param("C[Si](C)(C)O[Si](C)(C)c1ccccc1", "pentamethyl(phenyl)disiloxane", id="moving_the_phenyl_gives_the_same_compound"),
        pytest.param("C[SiH2]O[SiH2]C", "1,3-dimethyldisiloxane", id="two_placements_need_locants"),
        pytest.param(
            "C[Si](O[Si](C)(C)C1CC2C=CC1C2)(c1ccccc1)c1ccccc1",
            "1-(bicyclo[2.2.1]hept-5-en-2-yl)-1,1,3-trimethyl-3,3-diphenyldisiloxane",
            id="exchanging_groups_between_the_silicons_gives_isomers",
        ),
    ],
)
def test_locants_of_heteroatom_chains_follow_the_arrangement_rule(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_chalcogen_chain_with_a_double_bond_between_chalcogens_is_not_named_as_saturated():
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("CS(C)(C)(C)=S")
