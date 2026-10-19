import pytest
from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)CC[Si](C)(C)[Si](C)(C)CCC(O)=O",
            "3,3'-(1,1,2,2-tetramethyldisilane-1,2-diyl)dipropanoic acid",
        ),
    ],
)
def test_heteroatom_and_concatenated_linkers_between_chain_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCCOc1ccc(cc1)Cc1ccc(OCCO)cc1", "2,2'-[methylenebis(4,1-phenyleneoxy)]di(ethan-1-ol)"),
        (
            "OCCOCCOc1ccc(OCCOCCO)cc1",
            "2,2'-[1,4-phenylenebis(oxyethane-2,1-diyloxy)]di(ethan-1-ol)",
        ),
        (
            "CC(=O)[C@H]1CC=C(CC)C[C@@H]1C(C)=O",
            "1,1'-[(1S,2S)-4-ethylcyclohex-4-ene-1,2-diyl]di(ethan-1-one)",
        ),
    ],
)
def test_ring_and_assembly_components_inside_a_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [("CCn1c2ccccc2c2cc(Cc3ccc4c(c3)c3ccccc3n4CC)ccc21", "3,3'-methylenebis(9-ethyl-9H-carbazole)")],
)
def test_identical_n_substituted_fused_units_joined_by_a_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(=O)c1ccc(C(C)=O)cc1", "1,1'-(1,4-phenylene)di(ethan-1-one)"),
        ("CNCc1ccc(CNC)cc1", "1,1'-(1,4-phenylene)bis(N-methylmethanamine)"),
        ("CC(=O)c1cc(C(C)=O)cc(C(C)=O)c1", "1,1',1''-(benzene-1,3,5-triyl)tri(ethan-1-one)"),
        ("CCC(=O)N1CCN(C(=O)CC)CC1", "1,1'-(piperazine-1,4-diyl)di(propan-1-one)"),
        ("CC(=O)c1ccc(Oc2ccc(C(C)=O)cc2)cc1", "1,1'-[oxybis(4,1-phenylene)]di(ethan-1-one)"),
        ("O=C(c1ccccc1)c1ccc(C(=O)c2ccccc2)cc1", "(1,4-phenylene)bis(phenylmethanone)"),
        ("O=C(CC(=O)O)c1ccc(C(=O)CC(=O)O)cc1", "3,3'-(1,4-phenylene)bis(3-oxopropanoic acid)"),
        ("CC[C@H](C)C(=O)c1ccc(C(=O)[C@H](C)CC)cc1", "(2R,2'S)-1,1'-(1,4-phenylene)bis(2-methylbutan-1-one)"),
        ("CC(C)C(=O)OSSOC(=O)C(C)C", "1,1'-[disulfanediylbis(oxy)]bis(2-methylpropan-1-one)"),
    ],
)
def test_acyl_arms_cut_from_a_ring_or_chalcogen_core_are_ketones_not_aldehydes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCNCc1ccc(CNCC)cc1", "N,N'-[1,4-phenylenebis(methylene)]diethanamine"),
        ("CCN(C)Cc1ccc(CN(C)CC)cc1", "N,N'-[1,4-phenylenebis(methylene)]bis(N-methylethanamine)"),
        ("c1ccccc1NCc1ccc(CNc2ccccc2)cc1", "N,N'-[1,4-phenylenebis(methylene)]dianiline"),
        ("CCNCCCCNCC", "N1,N4-diethylbutane-1,4-diamine"),
    ],
)
def test_amine_parents_joined_through_nitrogen_are_multiplied_only_when_they_outnumber_the_substitutive_parent(
    smiles, expected
):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCC1CCC([2H])(CO)CC1", "[(1-2H)cyclohexane-1,4-diyl]dimethanol"),
        ("[2H]OCC1CCC(CO[2H])CC1", "(cyclohexane-1,4-diyl)di[methan(2H)ol]"),
        ("[13CH3]C(=O)c1ccc(C([13CH3])=O)cc1", "1,1'-(1,4-phenylene)di[(2-13C)ethan-1-one]"),
        ("[2H]OCC1CCC(CO)CC1", "[4-(hydroxymethyl)cyclohexyl]methan(2H)ol"),
        ("CC(=O)c1ccc(C([11CH3])=O)cc1", "1-(4-acetylphenyl)(2-11C)ethan-1-one"),
        ("OC[13CH]1CCC([13CH2]O)CC1", "[4-(hydroxymethyl)(4-13C)cyclohexyl](13C)methanol"),
        ("[2H]N(CCO)CCO", "2,2'-[(2H)azanediyl]di(ethan-1-ol)"),
        ("OCC[15N](CCO)CCO", "2,2',2''-[(15N)nitrilo]tri(ethan-1-ol)"),
        ("OC([2H])([2H])COCC([2H])([2H])O", "2,2'-oxydi[(1,1-2H2)ethan-1-ol]"),
    ],
)
def test_isotopic_modification_of_the_linker_or_of_identical_arms_keeps_the_multiplicative_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H](O)CN(C)C[C@H](C)O", "(2S,2'S)-1,1'-(methylazanediyl)di(propan-2-ol)"),
        ("S(/C1=C/CCCCCC1)/C1=C/CCCCCC/1", "(1E,1'E)-1,1'-sulfanediyldi(cyclooct-1-ene)"),
        ("C1=C/[C@@H](CC[C@@H]2/C=C/CCCCC2)CCCCC/1", "(1E,1'E,3R,3'R)-3,3'-(ethane-1,2-diyl)di(cyclooct-1-ene)"),
        ("[C@H]1(CC[C@H]3CCCC34CCCC4)CCCC12CCCC2", "(1R,1'R)-1,1'-(ethane-1,2-diyl)di(spiro[4.4]nonane)"),
        ("CC[C@@H](C)c1ccc(Sc2ccc([C@H](C)CC)cc2)cc1", "1,1'-sulfanediylbis{4-[(2R)-butan-2-yl]benzene}"),
    ],
)
def test_stereodescriptors_of_multiplied_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "C1=C\\[C@H](CC[C@@H]2/C=C/CCCCC2)CCCCC/1",
            "(1Z,3S)-3-{2-[(1R,2E)-cyclooct-2-en-1-yl]ethyl}cyclooct-1-ene",
        ),
        (
            "C1=C\\[C@@H](CC[C@H]2/C=C/CCCCC2)CCCCC/1",
            "(1Z,3R)-3-{2-[(1S,2E)-cyclooct-2-en-1-yl]ethyl}cyclooct-1-ene",
        ),
        ("[C@@H]1(CC[C@H]3CCCC34CCCC4)CCCC12CCCC2", "(1R)-1-{2-[(1S)-spiro[4.4]nonan-1-yl]ethyl}spiro[4.4]nonane"),
        ("[C@H]1(CC[C@@H]3CCCC34CCCC4)CCCC12CCCC2", "(1R)-1-{2-[(1S)-spiro[4.4]nonan-1-yl]ethyl}spiro[4.4]nonane"),
    ],
)
def test_units_of_different_configuration_take_the_senior_one_as_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(=O)NCCOCCNC(C)=O", "N,N'-[oxydi(ethane-2,1-diyl)]diacetamide"),
    ],
)
def test_amide_units_joined_through_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCN(CC)OBCCBON(CC)CC", "N,N'-[ethane-1,2-diylbis(boranediyloxy)]bis(N-ethylethanamine)"),
    ],
)
def test_amine_units_joined_through_nitrogen_to_a_heteroatom_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=C(O)CNC(=O)C(=O)NCC(=O)O", "2,2'-[oxalylbis(azanediyl)]diacetic acid"),
        ("OC(=O)CNC(=O)NCC(O)=O", "2,2'-[carbonylbis(azanediyl)]diacetic acid"),
        ("O=C(O)CNC(=O)CCC(=O)NCC(=O)O", "2,2'-butanediamidodiacetic acid"),
        ("OC(=O)c1ccccc1C(=O)CCC(=S)c1ccccc1C(=O)O", "2,2'-(1-oxo-4-sulfanylidenebutane-1,4-diyl)dibenzoic acid"),
        ("O=C(N)C(=O)NCCC(=O)O", "3-(oxamoylamino)propanoic acid"),
        ("OC(=O)c1ccc(cc1)C(=O)C1CCCCC1C(=O)c1ccc(cc1)C(O)=O", "4,4'-(cyclohexane-1,2-dicarbonyl)dibenzoic acid"),
        ("OC(=O)c1ccc(cc1)C(=S)CCC(=S)c1ccc(cc1)C(O)=O", "4,4'-[butanebis(thioyl)]dibenzoic acid"),
        ("OC(=O)c1ccc(cc1)NC(=S)CCC(=S)Nc1ccc(cc1)C(O)=O", "4,4'-(butanedithioamido)dibenzoic acid"),
        ("OC(=O)CNC(=O)c1ccccc1C(=O)NCC(=O)O", "2,2'-(benzene-1,2-dicarboxamido)diacetic acid"),
        ("OC(=O)CNC(=O)CC(C)C(=O)NCC(=O)O", "2,2'-(2-methylbutanediamido)diacetic acid"),
        ("O=C(CCC(=O)Nc1ccc(cc1)C(=O)O)Nc1ccc(cc1)C(=O)O", "4,4'-butanediamidodibenzoic acid"),
        ("O=S(=O)(CCS(=O)(=O)Nc1ccc(cc1)C(=O)O)Nc1ccc(cc1)C(=O)O", "4,4'-(ethane-1,2-disulfonamido)dibenzoic acid"),
        ("O=C(CCC(=O)N(C)c1ccc(cc1)C(=O)O)N(C)c1ccc(cc1)C(=O)O", "4,4'-[butanedioylbis(methylazanediyl)]dibenzoic acid"),
    ],
)
def test_diacyl_and_diamido_linking_groups_and_oxamoylamino(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(=O)N(CC(=O)O)CC(=O)O", "2,2'-(acetylazanediyl)diacetic acid"),
        ("CC(=O)N(CC(N)=O)CC(N)=O", "2-[acetyl(2-amino-2-oxoethyl)amino]acetamide"),
    ],
)
def test_acyl_group_on_a_linking_nitrogen_is_a_prefix_only_under_a_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("[SiH3][SiH2]C[SiH2][SiH3]", "1,1'-methylenebis(disilane)"),
        ("F[Si](F)(F)[Si](F)(F)C[Si](F)(F)[Si](F)(F)F", "1,1'-methylenebis(pentafluorodisilane)"),
        ("C([SiH2][SiH2]C)[SiH2][SiH2]C", "1,1'-methylenebis(2-methyldisilane)"),
        ("C([SiH](F)[SiH3])[SiH](F)[SiH3]", "1,1'-methylenebis(1-fluorodisilane)"),
        ("C([SiH2][SiH2]C)[SiH2][SiH2]CC", "1-[(2-ethyldisilanyl)methyl]-2-methyldisilane"),
        ("[SiH3]CC[SiH2]CC[SiH3]", "[silanediyldi(ethane-2,1-diyl)]bis(silane)"),
        ("C[Si](C)(C)C[Si](C)(C)C", "methylenebis(trimethylsilane)"),
        ("c1ccccc1P(c1ccccc1)CCP(c1ccccc1)c1ccccc1", "(ethane-1,2-diyl)bis(diphenylphosphane)"),
        ("CP(C)(=O)CCP(C)(C)=O", "(ethane-1,2-diyl)bis(dimethyl-\u03bb5-phosphanone)"),
        ("CP(C)(=S)CP(C)(C)=S", "methylenebis(dimethyl-\u03bb5-phosphanethione)"),
        (
            "O=P(c1ccccc1)(c1ccccc1)c1cc2ccccc2c(c1)-c1cc(P(=O)(c2ccccc2)c2ccccc2)cc2ccccc12",
            "[1,1'-binaphthalene]-3,3'-diylbis(diphenyl-\u03bb5-phosphanone)",
        ),
        (
            "c1ccccc1P(c1ccccc1)CP(c1ccccc1)CP(c1ccccc1)c1ccccc1",
            "[(phenylphosphanediyl)bis(methylene)]bis(diphenylphosphane)",
        ),
        ("C[Si](C)(C)C[Si](C)(C)C", "methylenebis(trimethylsilane)"),
    ],
)
def test_heteroatom_hydride_units_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_linker_with_a_substituent_that_carries_the_principal_group():
    assert (
        smiles_to_iupac("OC(=O)CN(CC(O)=O)CCN(CC(=O)O)CCN(CC(O)=O)CC(O)=O")
        == "2,2',2'',2'''-{[(carboxymethyl)azanediyl]bis(ethane-2,1-diylnitrilo)}tetraacetic acid"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCN(CCN(C)C)C", "N1-(2-aminoethyl)-N1,N2,N2-trimethylethane-1,2-diamine"),
        ("NCCNCCNCCNCN", "N1-{2-[(2-aminoethyl)amino]ethyl}-N2-(aminomethyl)ethane-1,2-diamine"),
        (
            "CC(CN(CC(=C)C)CC(=C)C)(C)C",
            "N-(2,2-dimethylpropyl)-2-methyl-N-(2-methylprop-2-en-1-yl)prop-2-en-1-amine",
        ),
        ("C1CNCCNCCN1", "1,4,7-triazonane"),
        ("C1NONNONNO1", "1,4,7,2,3,5,6,8-trioxapentazonane"),
    ],
)
def test_hetero_chain_and_macrocycle_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NCCNCCNCCNCCN", "N1,N1'-[azanediyldi(ethane-2,1-diyl)]di(ethane-1,2-diamine)", id="heteroatom_in_the_linker"),
        pytest.param("NCCCNCCc1ccc(cc1)CCNCCCN", "N1,N1'-[1,4-phenylenedi(ethane-2,1-diyl)]di(propane-1,3-diamine)", id="ring_in_the_linker"),
        pytest.param("NCCNCCNCCN", "N1,N1'-(ethane-1,2-diyl)di(ethane-1,2-diamine)", id="alkanediyl_linker"),
        pytest.param("NCCNCCOCCNCCN", "N1,N1'-[oxydi(ethane-2,1-diyl)]di(ethane-1,2-diamine)", id="oxygen_in_the_linker"),
        pytest.param("C(NCCN)NCCN", "N1,N1'-methylenedi(ethane-1,2-diamine)", id="methylene_linker_unchanged"),
        pytest.param("C(CNCCNCCNCCNCCN)N", "3,6,9,12-tetraazatetradecane-1,14-diamine", id="skeletal_replacement_for_four_heteroatoms"),
        pytest.param(
            "CNCCNCCNCCNC",
            "N1,N1'-(ethane-1,2-diyl)di(N2-methylethane-1,2-diamine)",
            id="amine_suffixes_kept_instead_of_skeletal_nitrogen_atoms_and_unit_numbered_from_the_junction",
        ),
        pytest.param(
            "CNCNCNCNC", "N',N'''-methylenebis(N-methylmethanediamine)", id="methylene_joined_units_numbered_from_the_junction"
        ),
        pytest.param(
            "CNCCOCCNCCOCCNC",
            "2,2'-[azanediylbis(ethane-2,1-diyloxy)]bis(N-methylethanamine)",
            id="terminal_secondary_amines_are_suffixes_not_skeletal_atoms",
        ),
    ],
)
def test_diamine_parents_joined_through_nitrogen_by_a_linking_group(smiles, expected):
    # P-62.2.4.1.3, P-62.2.5.2
    assert smiles_to_iupac(smiles) == expected


def test_amine_suffix_name_is_preferred_to_the_skeletal_replacement_name():
    assert smiles_to_iupac("CNCNCNCNC[SiH3]") == "N-({[(methylamino)methyl]amino}methyl)-N'-(silylmethyl)methanediamine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OB(O)c1ccc(B(O)O)cc1", "(1,4-phenylene)bis(boronic acid)", id="boronic_acids"),
        pytest.param("CB(O)CCB(C)O", "(ethane-1,2-diyl)bis(methylborinic acid)", id="borinic_acids"),
        pytest.param("CP(=O)(O)CCP(=O)(O)C", "(ethane-1,2-diyl)bis(methylphosphinic acid)", id="phosphinic_acids"),
        pytest.param("CP(O)CCP(C)O", "(ethane-1,2-diyl)bis(methylphosphinous acid)", id="phosphinous_acids"),
        pytest.param("OB(Cl)CCB(Cl)O", "(ethane-1,2-diyl)bis(boronochloridic acid)", id="halide_replaced_boronic_acids"),
        pytest.param("OB(N)CCB(N)O", "(ethane-1,2-diyl)bis(boronamidic acid)", id="amide_replaced_boronic_acids"),
        pytest.param("OBCCBO", "(ethane-1,2-diyl)bis(borinic acid)", id="borinic_acids_with_a_boron_hydrogen"),
        pytest.param("OP(=O)(Cl)CCP(=O)(O)Cl", "(ethane-1,2-diyl)bis(phosphonochloridic acid)", id="halide_replaced_phosphonic_acids"),
        pytest.param("OP(=S)(O)CCP(=S)(O)O", "(ethane-1,2-diyl)bis(phosphonothioic O,O-acid)", id="thio_replaced_phosphonic_acids"),
        pytest.param("OC(=O)N(C)CCN(C)C(=O)O", "(ethane-1,2-diyl)bis(methylcarbamic acid)", id="carbamic_acids_with_a_substituent"),
        pytest.param("SC(=S)Nc1ccc(NC(=S)S)cc1", "(1,4-phenylene)bis(carbamodithioic acid)", id="carbamodithioic_acids"),
        pytest.param("NC(=O)NCCNC(N)=O", "N,N''-(ethane-1,2-diyl)diurea", id="ureas"),
        pytest.param("CNC(=O)NCCNC(=O)NC", "N,N''-(ethane-1,2-diyl)bis(N'-methylurea)", id="far_nitrogen_substituent"),
        pytest.param("NC(=O)N(C)CCN(C)C(N)=O", "N,N''-(ethane-1,2-diyl)bis(N-methylurea)", id="joined_nitrogen_substituent"),
        pytest.param("NC(=S)Nc1ccc(NC(N)=S)cc1", "N,N''-(1,4-phenylene)bis(thiourea)", id="thioureas_take_bis"),
    ],
)
def test_identical_functional_parent_units_on_one_skeleton_are_multiplied(smiles, expected):
    # P-15.3.1.1(d), P-16.3.6; P-66.1.6.1.1.1: urea takes the locants N and N', not numbers
    assert smiles_to_iupac(smiles) == expected


def test_dithiocarboxy_prefix_is_enclosed_before_amino():
    assert smiles_to_iupac("SC(=S)NCC(=S)S") == "[(dithiocarboxy)amino]ethane(dithioic acid)"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[O]1CC[O](C)->[Pt]<-1([Cl])[Cl]", "dichlorido(1,2-dimethoxyethane-κ2O,O')platinum"),
        (
            "C1C[O]2->[K]3456<-[O]1CC[O]->3CC[O]->4CC[O]->5CC[O]->6CC2",
            "(1,4,7,10,13,16-hexaoxacyclooctadecane-κ6O,O',O'',O''',O''' ',O''' '')potassium",
        ),
    ],
)
def test_polydentate_ligand_complexes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C=C1CCC(=CC(=O)O)CC1", "2,2'-(cyclohexane-1,4-diylidene)diacetic acid"),
        (
            "OC(=O)C=C1C(Cl)=C2C=CC=CC2=CC1=CC(=O)O",
            "2,2'-(1-chloronaphthalene-2,3-diylidene)diacetic acid",
        ),
    ],
)
def test_linker_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CC(c1ccccc1)c1ccccc1", "diphenylacetaldehyde"),
        ("N#CC(c1ccccc1)c1ccccc1", "diphenylacetonitrile"),
        ("SC(c1ccccc1)c1ccccc1", "diphenylmethanethiol"),
        ("NC(=O)C(c1ccccc1)c1ccccc1", "2,2-diphenylacetamide"),
        ("OS(=O)(=O)C(c1ccccc1)c1ccccc1", "diphenylmethanesulfonic acid"),
    ],
)
def test_chain_parent_cites_each_aromatic_ring_as_a_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=Cc1ccc(Cc2ccccc2)cc1", "4-benzylbenzaldehyde"),
        ("N#Cc1ccc(Cc2ccccc2)cc1", "4-benzylbenzonitrile"),
        ("Sc1ccc(Cc2ccccc2)cc1", "4-benzylbenzene-1-thiol"),
        ("NC(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzamide"),
        ("OS(=O)(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzene-1-sulfonic acid"),
        ("NS(=O)(=O)c1ccc(Cc2ccccc2)cc1", "4-benzylbenzene-1-sulfonamide"),
        ("Nc1ccc(Cc2ccccc2)cc1", "4-benzylaniline"),
        ("Oc1ccc(Cc2ccccc2)cc1", "4-benzylphenol"),
        (
            "CC=1C=C(C=CC1)C=NNC1=CC(=NC(=N1)OCCC1=NC=CC=C1)N1CCOCC1",
            "4-(6-{[(3-methylphenyl)methylidene]hydrazinyl}-2-[2-(pyridin-2-yl)ethoxy]pyrimidin-4-yl)morpholine",
        ),
    ],
)
def test_ring_parent_cites_other_aromatic_rings_as_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "OCCc1cccc(CCO)c1SSc1c(C)cccc1C",
            "2,2'-{2-[(2,6-dimethylphenyl)disulfanyl]-1,3-phenylene}di(ethan-1-ol)",
        ),
        (
            "CC(O)CCc1cccc(CCC(C)O)c1SSc1c(COC)cccc1COC",
            "4,4'-(2-{[2,6-bis(methoxymethyl)phenyl]disulfanyl}-1,3-phenylene)di(butan-2-ol)",
        ),
    ],
)
def test_substituted_aryl_on_a_prefix_of_a_multiplied_ring_parent_is_its_own_compound_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_functional_group_on_a_heteroaromatic_ring_beside_another_ring_is_rejected_is_named():
    assert smiles_to_iupac("OC(=O)c1ccncc1Cc1ccccc1") == "3-benzylpyridine-4-carboxylic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1S(=O)(=O)c1ccccc1", "1,1'-sulfonyldibenzene"),
        ("c1ccccc1[Se](=O)(=O)c1ccccc1", "1,1'-selenonyldibenzene"),
        ("c1ccccc1C#Cc1ccccc1", "1,1'-(ethyne-1,2-diyl)dibenzene"),
        ("c1ccc2ccccc2c1Cc1cccc2ccccc12", "1,1'-methylenedinaphthalene"),
        ("BrC1=CC(=CC=C1)CC1=CC(=CC=C1)Cl", "1-bromo-3-[(3-chlorophenyl)methyl]benzene"),
    ],
)
def test_unsubstituted_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCc1ccc(Cc2ccc(CCl)c(Br)c2)cc1Br", "1,1'-methylenebis[3-bromo-4-(chloromethyl)benzene]"),
    ],
)
def test_substituted_units_and_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(Oc2ccc(C(=O)O)c(Br)c2)cc1Br", "4,4'-oxybis(2-bromobenzoic acid)"),
        ("Nc1ccc(Cc2ccc(N)cc2)cc1", "4,4'-methylenedianiline"),
        ("CN(C)c1ccc(Cc2ccc(N(C)C)cc2)cc1", "4,4'-methylenebis(N,N-dimethylaniline)"),
        ("CCN(C)c1ccc(Oc2ccc(N(C)CC)c(Cl)c2)cc1Cl", "4,4'-oxybis(2-chloro-N-ethyl-N-methylaniline)"),
        ("CNC1CCC(CC1)CC1CCC(CC1)NC", "4,4'-methylenebis(N-methylcyclohexan-1-amine)"),
        ("CNc1cc(NC)ccc1Cc1ccc(NC)cc1NC", "4,4'-methylenebis(N1,N3-dimethylbenzene-1,3-diamine)"),
        (
            "OC(=O)c1ccc(Cc2ccc(C(=O)O)c(C(=O)O)c2)cc1C(=O)O",
            "4,4'-methylenedi(benzene-1,2-dicarboxylic acid)",
        ),
        ("OC(=O)c1ccc(C(=O)c2ccc(C(=O)O)cc2)cc1", "4,4'-carbonyldibenzoic acid"),
        ("OC(=O)c1ccc(C(O)c2ccc(C(=O)O)cc2)cc1", "4,4'-(hydroxymethylene)dibenzoic acid"),
        ("O=C1CCC(CC1)CC1CCC(=O)CC1", "4,4'-methylenedi(cyclohexan-1-one)"),
        (
            "NC(=N)c1ccc(OCCCCCOc2ccc(C(N)=N)cc2)cc1",
            "4,4'-[pentane-1,5-diylbis(oxy)]di(benzene-1-carboximidamide)",
        ),
        (
            "NC(=N)c1ccc(CC(S(N)(=O)=O)c2ccc(C(N)=N)cc2)cc1",
            "4,4'-(1-sulfamoylethane-1,2-diyl)di(benzene-1-carboximidamide)",
        ),
    ],
)
def test_functionalized_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C1CCCN1CN1CCCC1=O", "1,1'-methylenedi(pyrrolidin-2-one)", id="saturated_lactam_units"),
        pytest.param("O=C(N1C=CC=CC1=O)N1C=CC=CC1=O", "1,1'-carbonyldi(pyridin-2(1H)-one)", id="mancude_lactam_units_with_added_hydrogen"),
        pytest.param("S=C(N1C=CC=CC1=O)N1C=CC=CC1=O", "1,1'-carbonothioyldi(pyridin-2(1H)-one)", id="thione_linker_between_ring_ketones_is_carbonothioyl"),
        pytest.param("O=c1ccc2ccccc2n1Cn1c(=O)ccc2ccccc21", "1,1'-methylenedi(quinolin-2(1H)-one)", id="fused_lactam_units"),
        pytest.param("O=C1CCN(C(=O)N2CCC(=O)C2)C1", "1,1'-carbonyldi(pyrrolidin-3-one)", id="saturated_heterocycle_unit_without_a_monocycle_spec"),
        pytest.param("OC(=O)C1CCCN1C(=O)N1CCCC1C(=O)O", "1,1'-carbonyldi(pyrrolidine-2-carboxylic acid)", id="acid_suffix_outranks_the_hidden_amide_linker"),
        pytest.param("S=C(c1ccc(C(=O)O)cc1)c1ccc(C(=O)O)cc1", "4,4'-carbonothioyldibenzoic acid", id="thione_linker_between_acid_units"),
    ],
)
def test_units_that_are_ring_pseudoketones_and_chalcogen_carbonyl_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("S=C1CCC(CC2CCC(=S)CC2)CC1", "4,4'-methylenedi(cyclohexane-1-thione)", id="carbocyclic_thione_units"),
        pytest.param("S=C1CCCN1CN1CCCC1=S", "1,1'-methylenedi(pyrrolidine-2-thione)", id="thiolactam_units"),
        pytest.param("S=C(N1C=CC=CC1=S)N1C=CC=CC1=S", "1,1'-carbonothioyldi(pyridine-2(1H)-thione)", id="carbonothioyl_between_thiolactam_units"),
        pytest.param("S=C1CC(=O)CC(CC2CC(=O)CC(=S)C2)C1", "3,3'-methylenebis(5-sulfanylidenecyclohexan-1-one)", id="thione_below_ketone_is_a_prefix"),
        pytest.param("S=C1CCC(CC2CCC(=O)CC2)CC1", "4-[(4-sulfanylidenecyclohexyl)methyl]cyclohexan-1-one", id="ketone_unit_outranks_thione_unit"),
    ],
)
def test_units_whose_principal_group_is_a_ring_thione_selone_or_tellone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(N(c2ccc(C(=O)O)cc2)c2ccc(C(=O)O)cc2)cc1", "4,4',4''-nitrilotribenzoic acid"),
    ],
)
def test_three_and_more_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CC=CCc1ccccc1", "1,1'-(but-2-ene-1,4-diyl)dibenzene"),
        ("c1ccccc1C(C1CCCCC1)c1ccccc1", "1,1'-(cyclohexylmethylene)dibenzene"),
        (
            "c1ccccc1CC(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1",
            "1,1'-(2,2-dibenzylpropane-1,3-diyl)dibenzene",
        ),
        (
            "c1ccccc1C(SCc1ccccc1)SCc1ccccc1",
            "1,1'-[(phenylmethylene)bis(sulfanediylmethylene)]dibenzene",
        ),
    ],
)
def test_central_group_forms(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Oc1ccc(OCC(C)COc2ccc(O)cc2)cc1", "4,4'-[(2-methylpropane-1,3-diyl)bis(oxy)]diphenol"),
        ("O=CCCCCOCCCCC=O", "5,5'-oxydi(pentanal)"),
        ("OC(=O)CCCCCCOCCCCCCC(O)=O", "7,7'-oxydi(heptanoic acid)"),
        ("O=CCCCOCCCC=O", "4,4'-oxydibutanal"),
        (
            "Oc1ccc(CCCCCCCCCCCCCCOCCCCCCCCCCCCCCc2ccc(O)cc2)cc1",
            "4,4'-[oxydi(tetradecane-14,1-diyl)]diphenol",
        ),
        (
            "C(COc1ccoc1)Oc1cocc1OCCOc1ccoc1",
            "3,3'-[furan-3,4-diylbis(oxyethane-2,1-diyloxy)]difuran",
        ),
    ],
)
def test_concatenated_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1/C=C/c1ccccc1", "1,1'-[(1E)-ethene-1,2-diyl]dibenzene"),
        ("C[C@H](c1ccccc1)[C@@H](C)c1ccccc1", "1,1'-[(2R,3S)-butane-2,3-diyl]dibenzene"),
    ],
)
def test_linker_stereodescriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccc(Cc2ccc(Cc3ccc(Cc4ccccc4)nc3)cc2)cc1", "3(2,5)-pyridina-1,7(1),5(1,4)-tribenzenaheptaphane"),
        pytest.param("c1ccc(C2CCC(C3CCC(c4ccc(C5CCC(C6CCC(C7CCCCC7)CC6)CC5)cc4)CC3)CC2)cc1",
            "1(1),4(1,4)-dibenzena-2,3,5,6(1,4),7(1)-pentacyclohexanaheptaphane", marks=pytest.mark.slow),
    ],
)
def test_linear_phane_with_any_ring_amplificants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_linear_phane_is_named_as_a_phane():
    assert (
        smiles_to_iupac("c1ccccc1Oc1cccc(Oc2cccc(Oc3ccccc3)c2)c1")
        == "2,4,6-trioxa-1,7(1),3,5(1,3)-tetrabenzenaheptaphane"
    )


def test_substituted_fused_unit():
    assert smiles_to_iupac("Clc1ccc2ccccc2c1Cc1c(Cl)ccc2ccccc12") == "1,1'-methylenebis(2-chloronaphthalene)"
    assert (
        smiles_to_iupac("C(=O)(NC1=CC=C2C=CC(=CC2=C1)S(=O)(=O)O)NC1=CC=C2C=CC(=CC2=C1)S(=O)(=O)O")
        == "7,7'-[carbonylbis(azanediyl)]di(naphthalene-2-sulfonic acid)"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccsc1Cc1ccco1", "2-[(thiophen-2-yl)methyl]furan"),
    ],
)
def test_unequal_rings_are_ranked_by_unsaturation_heteroatom_and_substituent_count(
    smiles, expected
):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OCCc1ccc(cc1)Oc1ccc(CCO)cc1", "2,2'-[oxybis(4,1-phenylene)]di(ethan-1-ol)"),
        ("OCCc1cc2ccccc2cc1CCO", "2,2'-(naphthalene-2,3-diyl)di(ethan-1-ol)"),
    ],
)
def test_chain_units_on_fused_and_composite_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC[C@H]1CCC[C@H]1CO", "[(1R,2S)-cyclopentane-1,2-diyl]dimethanol", id="r_takes_the_lower_locant"),
        pytest.param("NC[C@H]1CCC[C@@H]1CN", "[(1S,2S)-cyclopentane-1,2-diyl]dimethanamine", id="amine_units"),
        pytest.param("OC[C@H]1CCO[C@@H]1CO", "[(2S,3R)-oxolane-2,3-diyl]dimethanol", id="heteromonocycle_fixes_the_numbering"),
    ],
)
def test_stereodescriptors_of_a_ring_linker_stand_inside_its_enclosing_marks(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


_ARYL_ACID = "OC(=O)c1ccc(cc1)"
_ACID_ARYL = "c1ccc(cc1)C(=O)O"


@pytest.mark.parametrize(
    "linker, expected",
    [
        ("C1CCNCC1", "piperidine-3,4-diyl"),
        ("C1CC(NCC1)", "piperidine-2,4-diyl"),
        ("C1CCOC1", "oxolane-2,3-diyl"),
        ("C1COCC1", "oxolane-3,4-diyl"),
        ("c1ccc2cc(ccc2c1)", "naphthalene-2,6-diyl"),
        ("c1cc2cc(ccc2[nH]1)", "1H-indole-2,5-diyl"),
        ("c1ccc2nc(ccc2c1)", "quinoline-2,6-diyl"),
        ("C1CC2CCC1CC2", "bicyclo[2.2.2]octane-2,5-diyl"),
        ("C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2,7-diyl"),
        ("C1CCC2(CC1)CCCC2", "spiro[4.5]decane-1,8-diyl"),
        ("C1CC(C)CCC1", "4-methylcyclohexane-1,2-diyl"),
    ],
)
def test_ring_system_diyl_components_of_multiplicative_names(linker, expected):
    assert smiles_to_iupac(_ARYL_ACID + linker + _ACID_ARYL) == f"4,4'-({expected})dibenzoic acid"


@pytest.mark.parametrize(
    "linker, expected",
    [
        ("NN", "hydrazine-1,2-diyl"),
        ("N(C)N", "1-methylhydrazine-1,2-diyl"),
        ("N(C)N(C)", "1,2-dimethylhydrazine-1,2-diyl"),
        ("NNN", "triazane-1,3-diyl"),
        ("[PH][PH]", "diphosphane-1,2-diyl"),
        ("[PH][PH][PH]", "triphosphane-1,3-diyl"),
    ],
)
def test_homonuclear_heteroatom_hydride_linkers(linker, expected):
    assert smiles_to_iupac(_ARYL_ACID + linker + _ACID_ARYL) == f"4,4'-({expected})dibenzoic acid"


def test_diazenediyl_linker_cites_no_locants():
    assert smiles_to_iupac(_ARYL_ACID + "N=N" + _ACID_ARYL) == "4,4'-diazenediyldibenzoic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Oc1ccc(cc1)[As]=[As]c1ccc(O)cc1", "4,4'-diarsenediyldiphenol", id="ring_units"),
        pytest.param("OC(=O)CCP=PCCC(O)=O", "3,3'-diphosphenediyldipropanoic acid", id="chain_units"),
    ],
)
def test_doubly_bonded_group_15_pair_is_a_diylene_linker_without_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hydrazine_and_diazene_without_a_principal_group_are_the_parent():
    assert smiles_to_iupac("c1ccccc1NNc1ccccc1") == "1,2-diphenylhydrazine"
    assert smiles_to_iupac("c1ccccc1N=Nc1ccccc1") == "diphenyldiazene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccccc1/N=N\\c1ccccc1", "(Z)-diphenyldiazene"),
        ("C/N=N/c1ccccc1", "(E)-methyl(phenyl)diazene"),
        ("O/N=N\\O", "(Z)-diazene-1,2-diol"),
        ("C/N=N/O", "(E)-2-methyldiazen-1-ol"),
        ("N/N=N/N", "(2E)-tetraaz-2-ene"),
        ("O[SiH2][SiH2]O", "disilane-1,2-diol"),
    ],
)
def test_double_bond_descriptor_and_hydroxy_suffix_of_a_diazene_or_chain_hydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-41, P-44.1.2.2, P-68.3.1.3: a nitrogen heterocycle and a diazene are of one class, and the ring is the parent beside the diazenyl prefix
@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("c1ccc(N=Nc2ccccn2)cc1", "2-(phenyldiazenyl)pyridine", id="pyridine_position_2"),
        pytest.param("c1ccc(N=Nc2ccncc2)cc1", "4-(phenyldiazenyl)pyridine", id="pyridine_position_4"),
        pytest.param("c1ccc(N=Nc2cccnc2)cc1", "3-(phenyldiazenyl)pyridine", id="pyridine_position_3"),
        pytest.param("CN=Nc1ccccn1", "2-(methyldiazenyl)pyridine", id="alkyl_on_the_far_nitrogen"),
        pytest.param("N=Nc1ccccn1", "2-diazenylpyridine", id="unsubstituted_diazenyl"),
        pytest.param("Cc1ccc(N=Nc2ccccn2)cc1", "2-[(4-methylphenyl)diazenyl]pyridine", id="substituted_aryl_group"),
        pytest.param("Cc1cccc(N=Nc2ccccc2)n1", "2-methyl-6-(phenyldiazenyl)pyridine", id="ring_substituent_beside_the_diazenyl"),
        pytest.param("c1ccc(/N=N/c2ccccn2)cc1", "2-[(E)-phenyldiazenyl]pyridine", id="e_descriptor_inside_the_prefix"),
        pytest.param("c1ccc(/N=N\\c2ccccn2)cc1", "2-[(Z)-phenyldiazenyl]pyridine", id="z_descriptor_inside_the_prefix"),
        pytest.param("c1ccc(N=Nc2nccs2)cc1", "2-(phenyldiazenyl)-1,3-thiazole", id="thiazole"),
        pytest.param("c1ccc(N=Nc2ccc3ccccc3n2)cc1", "2-(phenyldiazenyl)quinoline", id="fused_nitrogen_ring"),
        pytest.param("c1ccc(N=Nc2ccc[nH]2)cc1", "2-(phenyldiazenyl)-1H-pyrrole", id="indicated_hydrogen"),
        pytest.param("c1ccc(N=Nc2cnc3ccccn23)cc1", "3-(phenyldiazenyl)imidazo[1,2-a]pyridine", id="fused_ring_with_a_bridgehead_nitrogen"),
        pytest.param("c1ccc(N=Nc2ccc(N=Nc3ccccc3)nc2)cc1", "2,5-bis(phenyldiazenyl)pyridine", id="two_diazenyl_groups"),
        pytest.param("c1ccc(N=Nc2ccccn2)cc1C(=O)O", "3-[(pyridin-2-yl)diazenyl]benzoic acid", id="acid_outranks_the_ring"),
        pytest.param("c1ccc(N=Nc2ccccc2)cc1", "diphenyldiazene", id="carbocycles_lose_to_the_diazene"),
        pytest.param("c1ccc(N=Nc2ccco2)cc1", "(furan-2-yl)(phenyl)diazene", id="oxygen_ring_loses_to_the_diazene"),
        pytest.param("CN(C)N=Nc1cccnc1", "3-[(dimethylamino)diazenyl]pyridine", id="amino_on_the_far_nitrogen_of_a_triazene"),
        pytest.param("CN(C)N=Nc1ccoc1", "1-(furan-3-yl)-3,3-dimethyltriaz-1-ene", id="oxygen_ring_loses_to_the_triazene"),
    ],
)
def test_diazene_beside_a_nitrogen_ring_is_a_prefix_of_the_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)CC1Cc2ccc(cc2)C(CC(=O)O)c2ccc1cc2",
            "2,2'-(1,3(1,4)-dibenzenacyclopentaphane-2,4-diyl)diacetic acid",
        ),
        (
            "OC(=O)CC1Cc2ccc(cc2)CC(CC(=O)O)c2ccc1cc2",
            "2,2'-(1,4(1,4)-dibenzenacyclohexaphane-2,6-diyl)diacetic acid",
        ),
    ],
)
def test_phane_diyl_between_chain_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "linker, expected",
    [
        ("[AsH][AsH]", "diarsane-1,2-diyl"),
        ("[SbH][SbH]", "distibane-1,2-diyl"),
        ("[SnH2][SnH2]", "distannane-1,2-diyl"),
    ],
)
def test_heavier_group_14_and_15_hydride_linkers(linker, expected):
    assert smiles_to_iupac(_ARYL_ACID + linker + _ACID_ARYL) == f"4,4'-({expected})dibenzoic acid"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)CC[AsH][AsH]CCC(=O)O", "3,3'-(diarsane-1,2-diyl)dipropanoic acid"),
        ("OC(=O)CC[PH][PH]CCC(=O)O", "3,3'-(diphosphane-1,2-diyl)dipropanoic acid"),
        ("OC(=O)CCNNCCC(=O)O", "3,3'-(hydrazine-1,2-diyl)dipropanoic acid"),
        ("OC(=O)CCN=NCCC(=O)O", "3,3'-diazenediyldipropanoic acid"),
    ],
)
def test_heteroatom_hydride_linkers_between_chain_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "OC(=O)c1ccc(cc1)C[SiH2]CC[SiH2]c1ccc(cc1)C(=O)O",
            "4-[(2-{[(4-carboxyphenyl)methyl]silyl}ethyl)silyl]benzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)[SiH2]C(F)C[SiH2]C(F)Cc1ccc(cc1)C(=O)O",
            "4-[(2-{[2-(4-carboxyphenyl)-1-fluoroethyl]silyl}-1-fluoroethyl)silyl]benzoic acid",
        ),
        (
            "OC(=O)c1ccc(cc1)OCC(Cl)c1ccc(C)cc1C(=O)O",
            "2-[2-(4-carboxyphenoxy)-1-chloroethyl]-5-methylbenzoic acid",
        ),
    ],
)
def test_unsymmetrical_linkers_and_units_are_named_substitutively(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1ccccc1C1CPCC(C1)c1ccccc1", "3,5-diphenylphosphinane"),
        ("c1ccccc1C1CC[BH]C1c1ccccc1", "2,3-diphenylborolane"),
    ],
)
def test_ring_seniority_ranks_every_heteroatom_of_the_element_order(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCCCC1[Sn](C)(C)C1CCCCC1O", "2,2'-(dimethylstannanediyl)di(cyclohexan-1-ol)"),
        ("OC1CCCCC1C(O)C1CCCCC1O", "2,2'-(hydroxymethylene)di(cyclohexan-1-ol)"),
        ("OC1CCCCC1C(O)C1CCCCC1", "2-[cyclohexyl(hydroxy)methyl]cyclohexan-1-ol"),
    ],
)
def test_substituted_metal_linker_and_principal_group_count_between_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_ring_unit_with_three_different_attachments_is_not_cited_as_a_triyl():
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("CC(c1cc(CCCCS)c(CCCCS)cc1)C(C)c1cc(CCCCS)c(CCCCS)cc1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CN(C)C(=O)CCSCCC(=O)N(C)C", "3,3'-sulfanediylbis(N,N-dimethylpropanamide)"),
        ("C(N1CCCCCCCCCCC1)N1CCCCCCCCCCC1", "1,1'-methylenebis(1-azacyclododecane)"),
        ("C(N1CCOCC1)N1CCOCC1", "4,4'-methylenedimorpholine"),
        ("C(CN=CCC(=O)O)N=CCC(=O)O", "3,3'-[ethane-1,2-diylbis(azanylylidene)]dipropanoic acid"),
        (
            "OC(=O)c1ccc(CC(C(=O)OC)c2ccc(C(O)=O)cc2)cc1",
            "4,4'-(3-methoxy-3-oxopropane-1,2-diyl)dibenzoic acid",
        ),
    ],
)
def test_identical_parents_with_substituted_units_hetero_rings_imine_and_ester_linkers_are_multiplied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CCP(CC1)CCP1CCCCC1", "1,1'-(ethane-1,2-diyl)bis(phosphinane)"),
        ("C1CCP(C1)CP1CCCC1", "1,1'-methylenebis(phospholane)"),
        ("C1CC[As](CC1)CC[As]1CCCCC1", "1,1'-(ethane-1,2-diyl)bis(arsinane)"),
        ("C1CCOC1CCC1CCCO1", "2,2'-(ethane-1,2-diyl)bis(oxolane)"),
        ("c1ccsc1CCc1cccs1", "2,2'-(ethane-1,2-diyl)dithiophene"),
        ("CP(C)CCP1CCCCC1", "1-[2-(dimethylphosphanyl)ethyl]phosphinane"),
    ],
)
def test_ring_units_with_a_replacement_name_are_multiplied_with_bis_unless_a_senior_hydride_atom_lies_outside_them(
    smiles, expected
):
    # P-15.3.2, P-16.3.6(c), P-44.1.2.2: 'diphosphinane' would be a ring of two phosphorus atoms
    assert smiles_to_iupac(smiles) == expected


def test_units_joined_to_the_linker_by_different_bonds_are_named_substitutively():
    # P-15.3.1.1, P-15.3.4.1.1: an -N=CH-CH2-NH- linker is unsymmetrical
    assert smiles_to_iupac("OCCN=CCNCCO") == "2-({2-[(2-hydroxyethyl)amino]ethylidene}amino)ethan-1-ol"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NCCCS(=O)CCCN", "3,3'-sulfinyldi(propan-1-amine)"),
        ("N(S(=O)CCCN)S(=O)CCCN", "3,3'-[azanediylbis(sulfinyl)]di(propan-1-amine)"),
    ],
)
def test_sulfinyl_linkers_of_multiplicative_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-15.3.1.2, P-41, P-67.1.4.1.1.2 to P-67.1.4.1.1.5: a phosphoryl-type group joins chains or rings that carry the principal group
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)CCP(=O)(CCC(=O)O)CCC(=O)O", "3,3',3''-phosphoryltripropanoic acid", id="three_chains_on_phosphoryl"),
        pytest.param("CP(=O)(CCC(=O)O)CCC(=O)O", "3,3'-(methylphosphonoyl)dipropanoic acid", id="organyl_group_makes_a_phosphonoyl"),
        pytest.param("OC(=O)CCP(C)(CCC(=O)O)=S", "3,3'-(methylphosphonothioyl)dipropanoic acid", id="chalcogen_analogue_of_phosphonoyl"),
        pytest.param("COC(=O)CCP(=O)(CCC(=O)OC)CCC(=O)OC", "trimethyl 3,3',3''-phosphoryltripropanoate", id="ester_chains"),
        pytest.param("OC(=O)CCP(=O)(Cl)CCC(=O)O", "3,3'-(chlorophosphoryl)dipropanoic acid", id="halogen_on_phosphoryl_beside_acid_chains"),
        pytest.param("OC(=O)c1ccc(cc1)P(=S)(O)c1ccc(cc1)C(=O)O", "4,4'-(hydroxyphosphorothioyl)dibenzoic acid", id="rings_on_a_substituted_phosphorothioyl"),
        pytest.param("C[As](=O)(CCC(=O)O)CCC(=O)O", "3,3'-(methylarsonoyl)dipropanoic acid", id="arsenic_analogue"),
        pytest.param("P(=O)(O)(CCC#N)CCC#N", "bis(2-cyanoethyl)phosphinic acid", id="phosphinic_acid_outranks_nitrile_chains"),
        pytest.param("P(=O)(OC)(c1ccc(cc1)O)c1ccc(cc1)O", "methyl bis(4-hydroxyphenyl)phosphinate", id="phosphinate_outranks_phenol_rings"),
    ],
)
def test_phosphoryl_type_groups_as_multiplying_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("OC(=O)c1ccc(cc1)S(=O)(=O)S(=O)(=O)c1ccc(cc1)C(O)=O", "4,4'-(disulfonyl)dibenzoic acid"),
        ("OC(=O)c1ccc(cc1)S(=O)SS(=O)c1ccc(cc1)C(O)=O", "4,4'-[sulfanediylbis(sulfinyl)]dibenzoic acid"),
        ("OC(=O)c1ccc(cc1)S(=O)(=O)SSS(=O)(=O)c1ccc(cc1)C(O)=O", "4,4'-[disulfanediylbis(sulfonyl)]dibenzoic acid"),
    ],
)
def test_runs_of_one_chalcogen_with_ligated_ends_are_multiplicative_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC[C@@H](C)c1ccc(Sc2ccc([C@H](C)CC)cc2)cc1",
            "1,1'-sulfanediylbis{4-[(2R)-butan-2-yl]benzene}",
        ),
        (
            "CC[C@@H](C)c1ccc(Sc2ccc([C@@H](C)CC)cc2)cc1",
            "1-[(2R)-butan-2-yl]-4-({4-[(2S)-butan-2-yl]phenyl}sulfanyl)benzene",
        ),
        (
            "C[C@@H](Br)C1([C@H](C)Br)CCCC1",
            "1-[(1R)-1-bromoethyl]-1-[(1S)-1-bromoethyl]cyclopentane",
        ),
    ],
)
def test_identical_configuration_units_are_multiplied_and_r_is_cited_before_s(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("n1nc(CO)c(Oc2c(CO)nncc2CO)c(CO)c1", "[oxydi(pyridazine-4,3,5-triyl)]tetramethanol"),
        ("n1nc(CO)c(CO)c(Oc2cnnc(CO)c2CO)c1", "[oxydi(pyridazine-5,3,4-triyl)]tetramethanol"),
        ("[SiH3]CCCCCCCCCCCCCCOCCCCCCCCCCCCCC[SiH3]", "[oxydi(tetradecane-14,1-diyl)]bis(silane)"),
        ("[SiH3]c1cc([SiH3])cc([SiH3])c1", "(benzene-1,3,5-triyl)tris(silane)"),
        (
            "OC(=O)Cc1c(CC(O)=O)c(CC(O)=O)c(CC(O)=O)c(CC(O)=O)c1",
            "2,2',2'',2''',2''' '-(benzene-1,2,3,4,5-pentayl)pentaacetic acid",
        ),
        ("CC(=O)N(C)[Si](C)(C=C)N(C)C(C)=O", "N,N'-[ethenyl(methyl)silanediyl]bis(N-methylacetamide)"),
    ],
)
def test_multiplicative_linker_locants_enclosure_and_multiplying_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Cl[SiH2]CC[SiH3]", "chloro(2-silylethyl)silane"),
        ("O=C(N=Nc1ccccc1)N=Nc1ccccc1", "bis(phenyldiazenyl)methanone"),
        ("S=C(N=N)N=N", "bis(diazenyl)methanethione"),
        ("[SiH3][SiH2]OCCS[SiH2][SiH3]", "3-oxa-6-thia-1,2,7,8-tetrasilaoctane"),
        ("[SiH3][SiH2]COO[SiH2][SiH3]", "[(disilanylmethyl)peroxy]disilane"),
        (
            "OC1CCC(CC1)=CCN=CCC1CCC(O)CC1",
            "4-(2-{[2-(4-hydroxycyclohexyl)ethylidene]amino}ethylidene)cyclohexan-1-ol",
        ),
        (
            "OC(=O)c1ccccc1CCOCCOCCOCCOCCc1ccccc1C(O)=O",
            "2,2'-(3,6,9,12-tetraoxatetradecane-1,14-diyl)dibenzoic acid",
        ),
        (
            "c1ccccc1OCCOCCOCCOCCOCCOc1ccccc1",
            "1,1'-[(3,6,9,12-tetraoxatetradecane-1,14-diyl)bis(oxy)]dibenzene",
        ),
        (
            "O(CCOCCOCCOCCOCCOc1ccoc1)c1cocc1OCCOCCOCCOCCOCCOc1ccoc1",
            "3,3'-{furan-3,4-diylbis[oxy(3,6,9,12-tetraoxatetradecane-14,1-diyl)oxy]}difuran",
        ),
        (
            "c1ccccc1OCCOCCOCCOCCOc1ccccc1",
            "1,1'-[oxybis(ethane-2,1-diyloxyethane-2,1-diyloxy)]dibenzene",
        ),
        (
            "OC(=O)c1ccc(cc1)C(c1cccc(c1)C(O)=O)(c1cccc(c1)C(O)=O)c1ccc(cc1)C(O)=O",
            "3,3'-[bis(4-carboxyphenyl)methylene]dibenzoic acid",
        ),
    ],
)
def test_parent_choice_beside_identical_hydrides_acyl_diazenes_and_heterounits(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Oc1ccc(cc1)[Se]O[Se]c1ccc(O)cc1", "4,4'-diselenoxanediyldiphenol", id="alternating_chalcogen_linker"),
        pytest.param("c1ccccc1OSOc1ccccc1", "1,1'-[sulfanediylbis(oxy)]dibenzene", id="senior_ends_keep_the_central_atom"),
    ],
)
def test_alternating_chalcogen_chains_are_divalent_parent_hydride_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param(
            "OC(=O)CCc1ccc(OCCOc2ccc(CCC(O)=O)cc2)cc1",
            "3,3'-[ethane-1,2-diylbis(oxy-4,1-phenylene)]dipropanoic acid",
            id="linking_group_then_arylene",
        ),
    ],
)
def test_locants_of_an_arylene_after_a_linking_group_follow_a_hyphen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param(
            "OC(=O)CCc1ccc(OCCOc2cccc(CCC(O)=O)c2)cc1",
            "3-(3-{2-[4-(2-carboxyethyl)phenoxy]ethoxy}phenyl)propanoic acid",
            id="linker_components_differ_on_the_two_sides",
        ),
        pytest.param(
            "OC(=O)CCOCC(C)OCCC(O)=O",
            "3,3'-[propane-1,2-diylbis(oxy)]dipropanoic acid",
            id="unsymmetrical_central_group_with_identical_branches",
        ),
    ],
)
def test_multiplicative_name_needs_identical_branches_on_the_linking_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param(
            "OC(=O)CCc1ccc(/C=C/O/C=C/c2ccc(CCC(O)=O)cc2)cc1",
            "3,3'-{oxybis[(1E)-ethene-2,1-diyl-4,1-phenylene]}dipropanoic acid",
            id="descriptor_of_the_linker_is_cited_once_before_its_component",
        ),
        pytest.param(
            "C[C@H](C(=O)O)c1ccc(OCCOc2ccc([C@@H](C)C(=O)O)cc2)cc1",
            "(2R,2'S)-2,2'-[ethane-1,2-diylbis(oxy-4,1-phenylene)]dipropanoic acid",
            id="units_cut_at_a_stereocentre",
        ),
        pytest.param(
            "OC(=O)CCc1ccc(/C=C/O/C=C\\c2ccc(CCC(O)=O)cc2)cc1",
            "3-{4-[(E)-2-({(Z)-2-[4-(2-carboxyethyl)phenyl]ethen-1-yl}oxy)ethen-1-yl]phenyl}propanoic acid",
            id="different_descriptors_in_the_linker_forbid_a_multiplicative_name",
        ),
        pytest.param(
            "C[C@H](O)c1ccc([C@@H](C)O)cc1",
            "(1R,1'S)-1,1'-(1,4-phenylene)di(ethan-1-ol)",
            id="arm_centres_beside_a_ring_linker",
        ),
        pytest.param(
            "C[C@H](O)CC1CCC(CC1)C[C@H](C)O",
            "(2S,2'S)-1,1'-(cyclohexane-1,4-diyl)di(propan-2-ol)",
            id="arm_centres_away_from_the_ring_linker",
        ),
        pytest.param(
            "OC[C@H](C)c1ccc(Cc2ccc([C@@H](C)CO)cc2)cc1",
            "(2R,2'R)-2,2'-[methylenebis(4,1-phenylene)]di(propan-1-ol)",
            id="arm_centres_beside_a_composite_linker",
        ),
        pytest.param(
            "OCCc1ccc(/C=C/C=C\\c2ccc(CCO)cc2)cc1",
            "2,2'-[(1Z,3E)-buta-1,3-diene-1,4-diyldi(4,1-phenylene)]di(ethan-1-ol)",
            id="double_bonds_of_a_central_chain",
        ),
        pytest.param(
            "OCCc1ccc([C@H](C)[C@@H](C)c2ccc(CCO)cc2)cc1",
            "2,2'-[(2R,3R)-butane-2,3-diyldi(4,1-phenylene)]di(ethan-1-ol)",
            id="centres_of_a_central_chain",
        ),
        pytest.param(
            "C[C@H](O)C(=O)OOOC(=O)[C@@H](C)O",
            "(2R,2'S)-1,1'-trioxidanediylbis(2-hydroxypropan-1-one)",
            id="arm_centres_beside_a_chalcogen_chain",
        ),
    ],
)
def test_stereodescriptors_in_a_multiplicative_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param(
            "OC(=O)c1ccc(cc1)S(=N)(=N)c1ccc(cc1)C(O)=O",
            "4,4'-sulfonodiimidoyldibenzoic acid",
            id="sulfonyl_centre_with_replaced_positions",
        ),
        pytest.param(
            "OC(=O)c1ccc(cc1)[Se](=S)(=N)c1ccc(cc1)C(O)=O",
            "4,4'-selenonimidothioyldibenzoic acid",
            id="centre_of_another_chalcogen_with_infixes_in_alphabetical_order",
        ),
        pytest.param(
            "CN=S(=O)(CC(=O)O)CC(=O)O",
            "2,2'-(N-methylsulfonimidoyl)diacetic acid",
            id="n_substituent_of_an_imidoyl_linker_between_chains",
        ),
        pytest.param(
            "OC(=O)c1ccc(cc1)S(=N)c1ccc(cc1)C(O)=O",
            "4,4'-sulfinimidoyldibenzoic acid",
            id="sulfinyl_type_centre_with_one_replaced_position",
        ),
        pytest.param(
            "OC(=O)c1ccc(cc1)S(=O)(=N)CCS(=O)(=N)c1ccc(cc1)C(O)=O",
            "4,4'-[ethane-1,2-diylbis(sulfonimidoyl)]dibenzoic acid",
            id="replaced_sulfonyl_group_in_an_arm",
        ),
        pytest.param(
            "OC(=O)c1ccc(cc1)N(c1ccc(cc1)C(O)=O)[N+](=O)[O-]",
            "4,4'-(nitroazanediyl)dibenzoic acid",
            id="nitro_group_on_the_linking_nitrogen",
        ),
        pytest.param(
            "OC(=O)c1ccc(cc1)[Sb](O)c1ccc(cc1)C(O)=O",
            "4,4'-(hydroxystibanediyl)dibenzoic acid",
            id="hydroxy_group_on_a_linking_group_15_atom",
        ),
    ],
)
def test_divalent_central_groups_with_replaced_chalcogen_positions_and_hetero_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("c1ncn(NCCNn2cncn2)n1", "N1,N2-bis(1H-1,2,4-triazol-1-yl)ethane-1,2-diamine", id="two_azol_1_amines_on_one_carbon_chain_are_its_diamine"),
        pytest.param("c1ncn(NCCOCCNn2cncn2)n1", "1,1'-[oxybis(ethane-2,1-diylazanediyl)]bis(1H-1,2,4-triazole)", id="two_azol_1_amines_joined_through_a_heteroatom_chain"),
        pytest.param("N(n1cnnc1)n1cnnc1", "N-(4H-1,2,4-triazol-4-yl)-4H-1,2,4-triazol-4-amine", id="single_nh_between_two_ring_nitrogens_is_an_amine"),
        pytest.param("c1ncn(N(C)n2cncn2)n1", "N-methyl-N-(1H-1,2,4-triazol-1-yl)-1H-1,2,4-triazol-1-amine", id="substituted_nitrogen_between_two_ring_nitrogens"),
    ],
)
def test_azolyl_units_joined_through_a_nitrogen_of_the_linker(smiles, expected):
    # P-62.2.1.3: a nitrogen on a ring nitrogen is the amine of the ring; P-44.1.1: a carbon chain holding both has the most principal groups
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(C)(c1ccc(S(=O)(=O)O)cc1)c1ccc(S(=O)(=O)O)cc1", "4,4'-(propane-2,2-diyl)di(benzene-1-sulfonic acid)", id="substituted_carbon_linker"),
        pytest.param("CN(c1ccc(S(=O)(=O)Cl)cc1)c1ccc(S(=O)(=O)Cl)cc1", "4,4'-(methylazanediyl)di(benzene-1-sulfonyl) dichloride", id="substituted_nitrogen_linker_of_acid_halide_units"),
    ],
)
def test_sulfonic_acid_units_joined_by_a_linker_that_carries_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C1CCCN1CC(C)CN1CCCC1=O", "1,1'-(2-methylpropane-1,3-diyl)di(pyrrolidin-2-one)", id="branched_alkanediyl_between_lactams"),
        pytest.param("O=C1CCCN1C[Si](C)(C)CN1CCCC1=O", "1,1'-[(dimethylsilanediyl)bis(methylene)]di(pyrrolidin-2-one)", id="substituted_silanediyl_between_lactams"),
        pytest.param("O=C1CCCN1c1ccc(-c2ccc(N3CCCC3=O)cc2)cc1", "1,1'-([1,1'-biphenyl]-4,4'-diyl)di(pyrrolidin-2-one)", id="ring_assembly_linker"),
        pytest.param("O=C1CCCN1c1ccc(-c2ccc(N3CCCC3=O)cc2)c(c1)-c1ccccc1", "1,1'-([1¹,2¹:2²,3¹-terphenyl]-1⁴,2⁴-diyl)di(pyrrolidin-2-one)", id="identical_ring_pendant_on_a_linking_ring_joins_the_assembly"),
        pytest.param("O=C1CCCN1c1ccc(-c2ccc(N3CCCC3=O)cc2)c(c1)-c1ccccn1", "1,1'-[2-(pyridin-2-yl)[1,1'-biphenyl]-4,4'-diyl]di(pyrrolidin-2-one)", id="other_ring_pendant_on_a_linking_assembly_is_a_substituent"),
        pytest.param("O=C(O)c1ccc(Oc2ccc(-c3ccc(Oc4ccc(C(=O)O)cc4)cc3)cc2)cc1", "4,4'-[[1,1'-biphenyl]-4,4'-diylbis(oxy)]dibenzoic acid", id="ring_assembly_between_oxy_arms"),
        pytest.param("c1ccc(cc1)Oc1ccc(-c2ccc(Oc3ccccc3)cc2)cc1", "4,4'-diphenoxy-1,1'-biphenyl", id="assembly_outranks_benzene_units_without_a_principal_group"),
    ],
)
def test_lactam_units_with_a_substituted_linker_and_ring_assembly_linkers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NNC(=O)c1ccc(Oc2ccc(C(=O)NN)cc2)cc1", "4,4'-oxydibenzohydrazide", id="hydrazide_units"),
        pytest.param("CC(C)(c1ccc(C(=O)NN)cc1)c1ccc(C(=O)NN)cc1", "4,4'-(propane-2,2-diyl)dibenzohydrazide", id="hydrazide_units_substituted_linker"),
        pytest.param("NC(=NO)c1ccc(Oc2ccc(C(N)=NO)cc2)cc1", "4,4'-oxybis(N'-hydroxybenzene-1-carboximidamide)", id="amidoxime_units"),
        pytest.param("CN(c1ccc(C(N)=NO)cc1)c1ccc(C(N)=NO)cc1", "4,4'-(methylazanediyl)bis(N'-hydroxybenzene-1-carboximidamide)", id="amidoxime_units_substituted_linker"),
        pytest.param("CNNC(=O)CCOCCC(=O)NNC", "3,3'-oxybis(N'-methylpropanehydrazide)", id="n_substituted_hydrazide_units"),
        pytest.param("CC(=NC)CCOCCC(C)=NC", "4,4'-oxybis(N-methylbutan-2-imine)", id="n_substituted_imine_units"),
    ],
)
def test_hydrazide_and_amidoxime_units_joined_by_a_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("Nc1ccc(N=Nc2ccc(N)cc2)cc1", "4,4'-diazenediyldianiline", id="azo_linker"),
        pytest.param("Cc1cc(ccc1N)N=Nc1ccc(N)c(C)c1", "4,4'-diazenediylbis(2-methylaniline)", id="azo_linker_between_substituted_units"),
        pytest.param("Nc1ccc(NNc2ccc(N)cc2)cc1", "4,4'-(hydrazine-1,2-diyl)dianiline", id="hydrazine_linker"),
        pytest.param("Nc1ccc(N=Nc2ccc(N)c(C)c2)cc1", "4-[(4-aminophenyl)diazenyl]-2-methylaniline", id="different_units_are_not_multiplied"),
    ],
)
def test_aniline_units_joined_through_a_nitrogen_linker_that_is_not_an_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)CC[18O][18O]CCC(=O)O", "3,3'-[(18O2)peroxy]dipropanoic acid", id="nuclides_of_a_linking_chalcogen_chain"),
        pytest.param(
            r"OC(=O)CC1CC/C=C\CC(CC(=O)O)CC1", "2,2'-[(6Z)-cyclonon-6-ene-1,4-diyl]diacetic acid", id="double_bond_of_a_ring_linker"
        ),
        pytest.param(
            "OC(=O)CC([C@H](F)Cl)COCC([C@H](F)Cl)CC(=O)O",
            "4,4'-oxybis{3-[(R)-chloro(fluoro)methyl]butanoic acid}",
            id="stereocentre_inside_the_substituent_of_each_unit",
        ),
    ],
)
def test_multiplied_units_with_nuclides_and_stereo_outside_their_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC(=O)CC(CCl)COCC(CCl)CC(=O)O", "4,4'-oxybis[3-(chloromethyl)butanoic acid]", id="unit_with_parentheses"),
    ],
)
def test_multiplied_unit_with_enclosing_marks_nests_them(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
