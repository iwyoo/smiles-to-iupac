import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._amino_acid import _SIDE_CHAIN_TABLE
from smiles_to_iupac._common import UnsupportedStructure


def test_side_chain_table_has_no_collisions():
    assert len(_SIDE_CHAIN_TABLE) == 24
    assert len(set(_SIDE_CHAIN_TABLE.values())) == 24


def test_glycine():
    assert smiles_to_iupac("NCC(=O)O") == "glycine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("NC(=N)NCCC[C@H](N)C(=O)O", "L-arginine"),
    ],
)
def test_arginine_specified_stereo_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C([C@H]([C@@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O)O", "β-D-glucofuranose"),  # CID 11309871
    ],
)
def test_cyclic_aldofuranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_pentofuranose_is_a_named_sugar():
    assert smiles_to_iupac("C([C@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O") == "α-L-arabinofuranose"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1[C@H]([C@H]([C@H]([C@](O1)(CO)O)O)O)O", "β-D-psicopyranose"),
    ],
)
def test_cyclic_ketohexopyranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", "α-D-glucopyranose", id="aldopyranose_still_resolves"),
        pytest.param("C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O", "β-D-fructofuranose", id="ketofuranose_is_a_named_sugar"),
        pytest.param("C1=C(NC=N1)C[C@H](C(=O)O)N", "D-histidine", id="d_histidine"),
        pytest.param("C1=C(NC=N1)CC(C(=O)O)N", "histidine", id="histidine_unspecified_stereocenter_no_ld_prefix"),
        pytest.param("Cn1cnc(CC(N)C(=O)O)c1", "2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid", id="ring_substituted_histidine_is_named_as_an_amino_acid_with_a_ring_prefix"),
        pytest.param("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "L-tryptophan", id="tryptophan_still_resolves"),
    ],
)
def test_aldopyranose_still_resolves_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1O", "neo-inositol"),  # (1,2,3/4,5,6-)
        ("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "myo-inositol"),  # (1,2,3,5/4,6-)
        ("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "scyllo-inositol"),  # (1,3,5/2,4,6-)
    ],
)
def test_inositol_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CO[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "1D-1-O-methyl-myo-inositol", id="o_substituent_lowest_locant_and_dl"),
        pytest.param("N[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "1D-1-amino-1-deoxy-myo-inositol", id="amino_deoxy_pair"),
        pytest.param("O[C@H]1[C@H](O)[C@@H](OC)[C@H](O)[C@@H](O)[C@H]1O", "5-O-methyl-myo-inositol", id="achiral_derivative_omits_dl"),
        pytest.param("F[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1Cl", "1D-2-chloro-1,2-dideoxy-1-fluoro-myo-inositol", id="halogeno_deoxy_pair_merged"),
        pytest.param("S[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1OC", "1D-1-deoxy-2-O-methyl-1-sulfanyl-allo-inositol", id="sulfanyl_deoxy_with_ether"),
        pytest.param("O[C@H]1[C@H](O)[C@H](OC)[C@@H](O)[C@H](O)[C@H]1O", "1D-2-O-methyl-chiro-inositol", id="chiro_enantiomer_from_numbering"),
        pytest.param("O[C@H]1[C@H](OCC)[C@@H](O)[C@H](O)[C@@H](OCC)[C@H]1O", "1L-1,4-di-O-ethyl-myo-inositol", id="multiplied_o_substituent"),
        pytest.param("O[C@H]1[C@H](OC(C)=O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "1L-myo-inositol 4-acetate", id="carboxylic_ester_named_as_alkanoate"),
        pytest.param("O[C@H]1[C@H](OS(=O)(=O)O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "1L-myo-inositol 4-(hydrogen sulfate)", id="sulfate_ester"),
        pytest.param(
            "O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O",
            "myo-inositol 4,5,6-tris(dihydrogen phosphate)",
            id="multiplied_phosphate_ester",
        ),
    ],
)
def test_inositol_derivatives(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_nucleoside_base_attached_via_oxygen_is_not_matched():
    name = smiles_to_iupac("C1=CN(C(=O)NC1=O)O[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O")
    assert name != "uridine"
    assert name == "1-{[(2S,3R,4S,5R)-3,4-dihydroxy-5-(hydroxymethyl)oxolan-2-yl]oxy}pyrimidine-2,4(1H,3H)-dione"


@pytest.mark.slow
def test_nucleoside_wrong_stereoisomer_is_not_matched():
    name = smiles_to_iupac("C1=CN(C(=O)NC1=O)[C@@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O")
    assert name != "uridine"
    assert name == "1-[(2S,3R,4S,5R)-3,4-dihydroxy-5-(hydroxymethyl)oxolan-2-yl]pyrimidine-2,4(1H,3H)-dione"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O",
            "5′-inosinic acid",
        ),
    ],
)
def test_nucleotide_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C1=NC2=C(N1[C@H]3C[C@@H]([C@H](O3)CO)O)N=C(N)N(C)C2=O",
            "2′-deoxy-1-methylguanosine",
        ),
        (
            "CCC1=CN(C(=O)NC1=S)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
            "5-ethyl-4-thiouridine",
        ),
        (
            "CC(=O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)OC(=O)C)OC(=O)C",
            "adenosine 2′,3′,5′-triacetate",
        ),
        (
            "CSC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C2N=C(NC3=O)NCCO)O)O",
            "N2-(2-hydroxyethyl)-5′-S-methyl-5′-thioguanosine",
        ),
        (
            "COC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=C(C(=NC2=O)N)I)F)O",
            "2′-deoxy-2′-fluoro-5-iodo-5′-O-methylcytidine",
        ),
        (
            "CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)COC)O",
            "5′-O-methylthymidine",
        ),
        (
            "C1C[C@@H](O[C@@H]1CO)N2C=NC3=C2N=C(NC3=O)N",
            "2′,3′-dideoxyguanosine",
        ),
        (
            "C1=CN(C(=O)N=C1N)[C@H]2/C(=C/F)[C@@H]([C@H](O2)CO)O",
            "(2′E)-2′-deoxy-2′-(fluoromethylidene)cytidine",
        ),
        (
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@H]4[C@@H]([C@H](O3)CO)OC(=O)O4)N",
            "2′,3′-dideoxyadenosine-2′,3′-diyl carbonate",
        ),
        (
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)[C@@H](C)O)O)O)N",
            "(5′R)-5′-C-methyladenosine",
        ),
        (
            "C[C@@]1(O)[C@H](O)[C@@H](CO)O[C@H]1n1cnc2c(N)ncnc21",
            "2′-C-methyladenosine",
        ),
        ("Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)C2(F)F)c(=O)n1", "2′-deoxy-2′,2′-difluorocytidine"),
        ("Nc1ncnc2c1ncn2[C@@H]1S[C@H](CO)[C@@H](O)[C@H]1O", "4′-thioadenosine"),
        ("COc1nc(N)nc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "O6-methylguanosine"),
        ("Cn1cnc2c(c1=N)ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "1-methyladenosine"),
        ("Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O", "3′-azido-3′-deoxythymidine"),
        (
            "CNC1=NC(=O)N(C=C1CCC(=O)O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
            "3-[4-(methylamino)-2-oxo-1-β-D-ribofuranosyl-1,2-dihydropyrimidin-5-yl]propanoic acid",
        ),
        (
            "CC(=O)NCc1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",
            "N-[(2,4-dioxo-1-β-D-ribofuranosyl-1,2,3,4-tetrahydropyrimidin-5-yl)methyl]acetamide",
        ),
        (
            "NC1=NC(=O)N(C=C1C#N)[C@H]2C[C@@H]([C@H](O2)COC(C)=O)O",
            "{(2R,3S,5R)-5-[4-amino-5-cyano-2-oxopyrimidin-1(2H)-yl]-3-hydroxyoxolan-2-yl}methyl acetate",
        ),
    ],
)
def test_substituted_nucleoside_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_substituted_nucleoside_wrong_sugar_stereo_is_not_matched():
    assert (
        smiles_to_iupac("CC1=CN(C(=O)NC1=O)[C@@H]2[C@H]([C@@H]([C@H](O2)CO)O)O")
        == "1-[(2S,3S,4S,5R)-3,4-dihydroxy-5-(hydroxymethyl)oxolan-2-yl]-5-methylpyrimidine-2,4(1H,3H)-dione"
    )


def test_plain_nucleoside_unaffected():
    assert (
        smiles_to_iupac("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N")
        == "adenosine"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C([C@H](C(=O)CO)O)O", "D-erythrulose"),  # CID 5460177
        ("C([C@@H]([C@H](C(=O)CO)O)O)O", "L-xylulose"),  # CID 22253
    ],
)
def test_open_chain_2_ketose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereo_2_ketose_still_falls_through_unchanged():
    assert smiles_to_iupac("C(C(C(C(C(=O)CO)O)O)O)O") == "1,3,4,5,6-pentahydroxyhexan-2-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C([C@H](C=O)O)O", "D-glyceraldehyde"),  # CID 79014
        ("C([C@@H](C=O)O)O", "L-glyceraldehyde"),  # CID 439723
    ],
)
def test_open_chain_aldose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C([C@H]([C@@H]([C@@H]([C@H]([C@@H](C=O)O)O)O)O)O)O",
            "D-glycero-L-gluco-heptose",
        ),  # CID 21139463 (mixed-series case)
    ],
)
def test_open_chain_heptose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C(C(C(C(C(C=O)O)O)O)O)O", "2,3,4,5,6-pentahydroxyhexanal", id="unspecified_stereo_aldose_still_falls_through_unchanged"),
        pytest.param("C1C[C@@H](NC1)C(=O)O", "D-proline", id="d_proline"),
        pytest.param("C1CC(NC1)C(=O)O", "proline", id="proline_unspecified_stereocenter_no_ld_prefix"),
        pytest.param("OC1C[C@H](NC1)C(=O)O", "4-hydroxy-L-proline", id="hydroxyproline_cites_the_specified_center"),
        pytest.param("O=C1CCCCO1", "oxan-2-one", id="hetero_ring_ketone_still_resolves"),
    ],
)
def test_unspecified_stereo_aldose_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-107.2 glycerides and glycerol phosphates
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('CCCCCCCCCCCCCCCCCC(=O)OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCC', 'propane-1,2,3-triyl trioctadecanoate'),
        ('CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)CO', '(2S)-2,3-dihydroxypropyl hexadecanoate'),
        ('CCCCCCCCCCCCCCCC(=O)OC(CO)CO', '1,3-dihydroxypropan-2-yl hexadecanoate'),
        ('CCCCCCCCCCCCCCCC(=O)OCC(O)COC(=O)CCCCCCCCCCCCCCC', '2-hydroxypropane-1,3-diyl dihexadecanoate'),
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@@H](OC(=O)CCCCCCCCCCCCCCC)CO', '(2S)-3-hydroxypropane-1,2-diyl dihexadecanoate', marks=pytest.mark.slow),
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@H](OC(C)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC', '(2S)-propane-1,2,3-triyl 2-acetate 1-hexadecanoate 3-[(9Z)-octadec-9-enoate]', marks=pytest.mark.slow),
        ('OC[C@H](O)COP(O)(O)=O', '(2S)-2,3-dihydroxypropyl dihydrogen phosphate'),
        ('OC[C@@H](O)COP(O)(O)=O', '(2R)-2,3-dihydroxypropyl dihydrogen phosphate'),
        ('OCCOP(O)(=O)OCCO', 'bis(2-hydroxyethyl) hydrogen phosphate'),
    ],
)
def test_glycerides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-107.3 phosphatidic acids and phosphoglycerides
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(O)=O)OC(=O)CCCCCCCCCCCCCCC', '(2S)-3-(phosphonooxy)propane-1,2-diyl dihexadecanoate', marks=pytest.mark.slow),
        ('CCCCCCCCCCCCCCCC(=O)OCC(OP(O)(O)=O)COC(=O)CCCCCCCCCCCCCCC', '2-(phosphonooxy)propane-1,3-diyl dihexadecanoate'),
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OCCN)OC(=O)CCCCCCCCCCCCCCC', '(2S)-3-{[(2-aminoethoxy)hydroxyphosphoryl]oxy}propane-1,2-diyl dihexadecanoate', marks=pytest.mark.slow),
        ('CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC', 'O-{[(2R)-2,3-bis(octadecanoyloxy)propoxy]hydroxyphosphoryl}-L-serine'),
        ('CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCC', 'O-[((2R)-2-(hexadecanoyloxy)-3-{[(9Z)-octadec-9-enoyl]oxy}propoxy)hydroxyphosphoryl]-L-serine'),
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)OC[C@H](O)CO)OC(=O)CCCCCCCCCCCCCCC', '(2S)-3-({[(2R)-2,3-dihydroxypropoxy]hydroxyphosphoryl}oxy)propane-1,2-diyl dihexadecanoate', marks=pytest.mark.slow),
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCC', '(2S)-3-[(hydroxy{[(1S,2R,3R,4S,5S,6R)-2,3,4,5,6-pentahydroxycyclohexyl]oxy}phosphoryl)oxy]propane-1,2-diyl dihexadecanoate', marks=pytest.mark.slow),
        pytest.param('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O)OC(=O)CCCCCCCCCCCCCCC', '(2S)-3-[(hydroxy{[(1s,2R,3S,4s,5R,6S)-2,3,4,5,6-pentahydroxycyclohexyl]oxy}phosphoryl)oxy]propane-1,2-diyl dihexadecanoate', marks=pytest.mark.slow),
        ('CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O', '(2R)-2-hydroxy-3-(phosphonooxy)propyl hexadecanoate'),
        ('CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN', '(2R)-3-{[(2-aminoethoxy)hydroxyphosphoryl]oxy}-2-hydroxypropyl hexadecanoate'),
        ('CCCCCCCCCCCCCCCCOC[C@H](OC(C)=O)COP(O)(O)=O', '(2S)-1-(hexadecyloxy)-3-(phosphonooxy)propan-2-yl acetate'),
        ('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP(O)([O-])=O)OC(=O)CCCCCCCCCCCCCCC.[Na+]', 'sodium (2S)-2,3-bis(hexadecanoyloxy)propyl hydrogen phosphate'),
        ('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP([O-])(=O)OC[C@H](O)CO)OC(=O)CCCCCCCCCCCCCCC', '(2S)-2,3-bis(hexadecanoyloxy)propyl (2R)-2,3-dihydroxypropyl phosphate'),
        ('CCCCCCCCCCCCCCCC(=O)OC[C@@H](COP([O-])(=O)OCC[NH3+])OC(=O)CCCCCCCCCCCCCCC', '2-azaniumylethyl (2S)-2,3-bis(hexadecanoyloxy)propyl phosphate'),
    ],
)
def test_phosphoglycerides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-107.3.4 phosphatidylcholines, cationic parent with skeletal replacement
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('C[N+](C)(C)CCOP(O)(=O)OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC', '(7R)-7-(hexadecanoyloxy)-4-hydroxy-N,N,N-trimethyl-4,10-dioxo-3,5,9-trioxa-4λ5-phosphapentacosan-1-aminium'),
        ('C[N+](C)(C)CCOP(O)(=O)OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC.[OH-]', '(7R)-7-(hexadecanoyloxy)-4-hydroxy-N,N,N-trimethyl-4,10-dioxo-3,5,9-trioxa-4λ5-phosphapentacosan-1-aminium hydroxide'),
        ('C[N+](C)(C)CCOP([O-])(=O)OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC', '(2R)-2,3-bis(hexadecanoyloxy)propyl 2-(trimethylazaniumyl)ethyl phosphate'),
        ('C[N+](C)(C)CCOP([O-])(=O)OC[C@@H](O)COC(=O)CCCCCCCCCCCCCCC', '(2S)-3-(hexadecanoyloxy)-2-hydroxypropyl 2-(trimethylazaniumyl)ethyl phosphate'),
        ('CCCCCCCCCCCCCCCCOC[C@H](OC(C)=O)COP([O-])(=O)OCC[N+](C)(C)C', '(2S)-2-(acetyloxy)-3-(hexadecyloxy)propyl 2-(trimethylazaniumyl)ethyl phosphate'),
        ('C[N+](C)(C)CC.[OH-]', 'N,N,N-trimethylethanaminium hydroxide'),
        ('C[N+](C)(C)CCOP(O)(O)=O', 'N,N,N-trimethyl-2-(phosphonooxy)ethan-1-aminium'),
    ],
)
def test_phosphatidylcholines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-107.4.3 sphingoids, ceramides and sphingophospholipids
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('CCCCCCCCCCCCCCC[C@H]([C@H](CO)N)O', 'sphinganine'),
        ('CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO)N)O', '(4E)-sphing-4-enine'),
        ('CCC/C=C/CCCCCCCC/C=C/[C@H]([C@H](CO)N)O', '(4E,14E)-sphinga-4,14-dienine'),
        ('CCCCCCCCCCCCCCC[C@H]([C@H](CO)NC)O', 'N-methylsphinganine'),
        ('CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO)N(C)C)O', '(4E)-N,N-dimethylsphing-4-enine'),
        ('CCCCCCCCCCCCCCC[C@H]([C@H](COC)N)OC', '1,3-di-O-methylsphinganine'),
        ('CCCCCCCCCCCCCCC[C@H]([C@H](COC)NCC)O', 'N-ethyl-1-O-methylsphinganine'),
        ('CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)N)O', '(4E)-1-O-(β-D-galactopyranosyl)sphing-4-enine'),
        ('CCCCCCCCCCCCCCCCC[C@H]([C@H](CO)N)O', '(2S,3R)-2-aminoicosane-1,3-diol'),
        ('CCCCCCCCCCCCCCC[C@@H]([C@H](CO)N)O', '(2S,3S)-2-aminooctadecane-1,3-diol'),
        pytest.param('CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO)NC(=O)CCCCCCCCCCCCCCC)O', 'N-[(2S,3R,4E)-1,3-dihydroxyoctadec-4-en-2-yl]hexadecanamide', marks=pytest.mark.slow),
        pytest.param('CCCCCCCCCCCCC/C=C/[C@H]([C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCCCCCCCCC)O', 'N-[(2S,3R,4E)-1-(β-D-galactopyranosyloxy)-3-hydroxyoctadec-4-en-2-yl]hexadecanamide', marks=pytest.mark.slow),
        ('CCCCCCCCCCCCC/C=C/[C@H]([C@H](COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCCC)O', '(2S,3R,4E)-2-(hexadecanoylamino)-3-hydroxyoctadec-4-en-1-yl 2-(trimethylazaniumyl)ethyl phosphate'),
        ('CCCCCCCCCCCCC/C=C/[C@H]([C@H](COP(O)(O)=O)N)O', '(2S,3R,4E)-2-amino-3-hydroxyoctadec-4-en-1-yl dihydrogen phosphate'),
        ('CCCCCCCCCCCCC/C=C/[C@H]([C@H](COP(O)(O)=O)NC(=O)CCCCCCCCCCCCCCC)O', '(2S,3R,4E)-2-(hexadecanoylamino)-3-hydroxyoctadec-4-en-1-yl dihydrogen phosphate'),
    ],
)
def test_sphingolipids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-107.4.2 glycoglycerolipids
@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('CCCCCCCCCCCCCCCCCC(=O)OC[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)OC(=O)CCCCCCCCCCCCCCCCC', '(2S)-3-(β-D-galactopyranosyloxy)propane-1,2-diyl dioctadecanoate'),
    ],
)
def test_glycoglycerolipids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-102.5.6.2, P-102.6.1, P-102.7 glycosides, glycosyl compounds and oligosaccharides
@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC[C@H]1O[C@@H](OC)[C@H](O)[C@@H](O)[C@@H]1O", "methyl β-D-glucopyranoside", id="methyl_glycoside"),
        pytest.param("CCO[C@]1(CO)OC[C@@H](O)[C@@H](O)[C@@H]1O", "ethyl β-D-fructopyranoside", id="ketopyranoside"),
        pytest.param("Oc1ccc(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc1", "4-hydroxyphenyl β-D-glucopyranoside", id="aglycone_with_hydroxy_group"),
        pytest.param("CC(=O)c1ccc(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc1", "1-[4-(β-D-glucopyranosyloxy)phenyl]ethan-1-one", id="senior_aglycone_keeps_prefix_form"),
        pytest.param("OC[C@H]1O[C@@H](F)[C@H](O)[C@@H](O)[C@@H]1O", "β-D-glucopyranosyl fluoride", id="glycosyl_halide"),
        pytest.param("OC[C@H]1O[C@@H](N)[C@H](O)[C@@H](O)[C@@H]1O", "β-D-glucopyranosylamine", id="glycosylamine"),
        pytest.param("OC[C@H]1O[C@@H](Nc2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O", "N-phenyl-β-D-glucopyranosylamine", id="n_substituted_glycosylamine"),
        pytest.param("OC[C@H]1O[C@H](c2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O", "(β-D-glucopyranosyl)benzene", id="c_glycosyl_ring"),
        pytest.param("C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@H](O2)O[C@]3([C@H]([C@@H]([C@H](O3)CO)O)O)CO)O)O)O)O)O)O)O", "β-D-fructofuranosyl α-D-galactopyranosyl-(1→6)-α-D-glucopyranoside", id="trisaccharide_without_hemiacetal"),
        pytest.param("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CO)O)O)CO)O)O)O)O", "β-D-fructofuranosyl α-D-glucopyranoside", id="aldose_is_the_glycoside_parent"),
        pytest.param("OC[C@H]1O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "α-D-glucopyranosyl α-D-glucopyranoside", id="identical_units_anomeric_bond"),
        pytest.param("OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](OC)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O", "methyl α-D-glucopyranosyl-(1→4)-β-D-glucopyranoside", id="aglycone_on_oligosaccharide"),
        pytest.param("OC[C@H]1O[C@H](OC[C@H]2O[C@H](O[C@H]3[C@H](O)[C@@H](O)C(O)O[C@@H]3CO)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "α-D-glucopyranosyl-(1→6)-α-D-glucopyranosyl-(1→4)-D-glucopyranose", id="reducing_end_without_anomeric_descriptor"),
        pytest.param("C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O", "β-D-fructofuranose", id="ketofuranose"),
    ],
)
def test_glycosides_and_oligosaccharides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-102.5.3 to P-102.5.6.4 deoxy, amino, thio, halo, O- and C-substituted monosaccharides
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("C([C@H]([C@H](CC=O)O)O)O", "2-deoxy-D-erythro-pentose", id="deoxy_removes_a_centre_systematic_prefix"),
        pytest.param("C([C@H]([C@H]([C@@H](CC=O)O)O)O)O", "2-deoxy-D-arabino-hexose", id="deoxy_open_chain_hexose"),
        pytest.param("C([C@@H]1[C@H](CC(O1)O)O)O", "2-deoxy-D-erythro-pentofuranose", id="furanose_with_unspecified_anomeric_carbon"),
        pytest.param("C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)O)O)O)O", "α-L-rhamnopyranose", id="six_deoxy_retained_names"),
        pytest.param("C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)N)O)O)O", "2-amino-2-deoxy-β-D-glucopyranose", id="amino_deoxy_pair"),
        pytest.param("CC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1O)CO)O)O", "2-(acetylamino)-2-deoxy-β-D-glucopyranose", id="substituted_amino_group"),
        pytest.param("BrC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "6-bromo-6-deoxy-β-D-glucopyranose", id="halogen_with_deoxy"),
        pytest.param("OC[C@H]1S[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "5-thio-β-D-glucopyranose", id="ring_sulfur"),
        pytest.param("COC[C@H]1O[C@@H](O)[C@H](OC)[C@@H](OC)[C@@H]1OC", "2,3,4,6-tetra-O-methyl-β-D-glucopyranose", id="multiplied_o_alkyl_prefix"),
        pytest.param("COC[C@H]1O[C@@H](O)[C@H](OC(C)=O)[C@@H](O)[C@@H]1O", "6-O-methyl-β-D-glucopyranose 2-acetate", id="ester_after_the_name"),
        pytest.param("OC1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H](O)[C@H]1O", "D-glucopyranose 6-(dihydrogen phosphate)", id="phosphate_ester"),
        pytest.param("OC[C@H]1O[C@@H](O)[C@@](O)(c2ccccc2)[C@@H](O)[C@@H]1O", "2-C-phenyl-β-D-mannopyranose", id="c_substituent_on_a_nonterminal_carbon"),
        pytest.param("O=C[C@@H](OC(=O)c1ccccc1)[C@@H](OC(=O)c1ccccc1)[C@H](OC(=O)c1ccccc1)[C@H](OC(=O)c1ccccc1)COC(=O)c1ccccc1", "D-mannose 2,3,4,5,6-pentabenzoate", id="open_chain_polyester"),
        pytest.param("C([C@@H]1[C@H]([C@H]([C@@H](O1)O)O)O)O", "β-D-ribofuranose", id="plain_pentofuranose"),
        pytest.param("OP(O)(=O)OC[C@H]1O[C@](O)(COP(O)(O)=O)[C@@H](O)[C@@H]1O", "β-D-fructofuranose 1,6-bis(dihydrogen phosphate)", id="ketose_with_multiplied_esters"),
    ],
)
def test_substituted_monosaccharides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-66.1.1.1.2.5, P-66.3.1.2.4, P-66.5.1.2.4 amides, hydrazides, nitriles and esters of carbohydrate acids
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NC(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "D-gluconamide", id="aldonamide"),
        pytest.param("NNC(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "D-gluconohydrazide", id="aldonohydrazide"),
        pytest.param("N#C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "D-glucononitrile", id="aldononitrile"),
        pytest.param("CC(C)OC(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "propan-2-yl D-gluconate", id="aldonate_ester"),
        pytest.param("NC(=O)[C@H]1O[C@@H](OC)[C@H](O)[C@@H](O)[C@@H]1O", "methyl β-D-glucopyranosiduronamide", id="glycosiduronamide"),
        pytest.param("NC(=O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "α-D-glucopyranuronamide", id="pyranuronamide"),
    ],
)
def test_carbohydrate_acid_derivatives(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-102.5.6.5, P-102.5.6.6 alditols and monosaccharide carboxylic acids
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO", "D-glucitol", id="alditol_parent_chosen_by_alphabetical_stem"),
        pytest.param("OC[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO", "galactitol", id="meso_alditol_has_no_dl"),
        pytest.param("O=C(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "D-gluconic acid", id="aldonic_acid"),
        pytest.param("O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C(=O)O", "D-glucuronic acid", id="uronic_acid"),
        pytest.param("O=C(O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O", "D-glucaric acid", id="aldaric_acid"),
        pytest.param("OC(=O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "α-D-glucopyranuronic acid", id="cyclic_uronic_acid"),
        pytest.param("OC(=O)[C@H]1O[C@@H](OC)[C@H](O)[C@@H](O)[C@@H]1O", "methyl β-D-glucopyranosiduronic acid", id="uronic_acid_glycoside"),
    ],
)
def test_alditols_and_sugar_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-103.2.3 N-, S- and O-substitution of retained amino acid names
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('C[C@H](NC(C)=O)C(=O)O', 'N-acetyl-L-alanine'),
        ('N[C@@H](CSCc1ccccc1)C(=O)O', 'S-benzyl-L-cysteine'),
        ('N[C@@H](Cc1ccc(OC)cc1)C(=O)O', 'O-methyl-L-tyrosine'),
        ('CC(=O)NCCCC[C@H](NC)C(=O)O', 'N6-acetyl-N2-methyl-L-lysine'),
        ('CN(C)[C@@H](C)C(=O)O', 'N,N-dimethyl-L-alanine'),
        ('CC(=O)NC(=O)CC[C@H](NC(=O)OCc1ccccc1)C(=O)O', 'N5-acetyl-N2-[(benzyloxy)carbonyl]-L-glutamine'),
        ('CNC(=O)C[C@H](N)C(=O)O', 'N4-methyl-L-asparagine'),
        ('CN=C(NC)NCCC[C@H](NC(C)=O)C(=O)O', 'Nα-acetyl-Nω,Nω′-dimethyl-L-arginine'),
        ('NC(=N)N(C)CCC[C@H](N)C(=O)O', 'Nδ-methyl-L-arginine'),
    ],
)
def test_substituted_amino_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-103.1.3.2.2 allo prefix; P-103.2.4.2 ions; P-103.2.6 esters
@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC[C@H](C)[C@H](N)C(=O)O", "L-isoleucine", id="isoleucine"),
        pytest.param("CC[C@H](C)[C@@H](N)C(=O)O", "D-alloisoleucine", id="allo_inverts_beta_in_d_series"),
        pytest.param("C[C@@H](O)[C@H](N)C(=O)O", "L-threonine", id="threonine"),
        pytest.param("CC(=O)N[C@@H]([C@H](C)O)C(=O)O", "N-acetyl-L-allothreonine", id="substituted_allo"),
        pytest.param("CC[C@H](C)C(N)C(=O)O", "(3S)-2-amino-3-methylpentanoic acid", id="beta_without_alpha_stays_systematic"),
        pytest.param("C[C@H](NC(C)=O)C(=O)OC", "methyl N-acetyl-L-alaninate", id="ester_of_substituted_acid"),
        pytest.param("CC(C)[C@H](N)C(=O)OC[C@@H](N)C(C)C", "(2S)-2-amino-3-methylbutyl L-valinate", id="ester_group_with_stereo"),
        pytest.param("OC(=O)C[C@H](N)C(=O)OC", "1-methyl L-aspartate", id="diacid_monoester_locant"),
        pytest.param("COC(=O)C[C@H](N)C(=O)OCC", "1-ethyl 4-methyl L-aspartate", id="diacid_mixed_esters"),
        pytest.param("COC(=O)C[C@H](N)C(=O)OC", "dimethyl L-aspartate", id="diacid_identical_esters"),
        pytest.param("C[C@H](N)C(=O)[O-]", "L-alaninate", id="anion"),
        pytest.param("[NH3+]CC(=O)O", "glycinium", id="cation"),
        pytest.param("N[C@@H](CCC(=O)[O-])C(=O)[O-]", "L-glutamate", id="diacid_dianion"),
    ],
)
def test_amino_acid_esters_ions_and_allo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)c1ccc(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc1", "4-(β-D-glucopyranosyloxy)benzoic acid"),
        (
            "OC(=O)c1ccc(O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)cc1",
            "4-{[(2R,3S,4R,5R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy}benzoic acid",
        ),
        (
            "OC(=O)c1ccc(O[C@H]2[C@H](O)C(O)O[C@@H](CO)[C@@H]2O)cc1",
            "4-{[(3S,4R,5S,6S)-2,3,5-trihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy}benzoic acid",
        ),
    ],
)
def test_sugar_substituent_is_glycosyl_only_when_linked_at_the_anomeric_carbon(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1[C@@H]2[C@H]([C@@H]([C@H]([C@H](O1)O2)O)O)O", "1,6-anhydro-β-D-glucopyranose", id="levoglucosan"),
        pytest.param("C([C@@H]1[C@H]([C@@H]([C@H](O1)C=O)O)O)O", "2,5-anhydro-D-mannose", id="aldose_with_a_furan_bridge"),
        pytest.param("CO[C@H]1[C@@H]([C@H](C=O)OC)OC[C@H]1OC", "3,6-anhydro-2,4,5-tri-O-methyl-D-glucose", id="anhydro_sorts_before_the_O_substituents"),
        pytest.param("CO[C@H]1[C@@H]([C@@H](C=O)OC)OC[C@H]1OC", "3,6-anhydro-2,4,5-tri-O-methyl-D-mannose", id="anhydro_mannose_ether"),
        pytest.param("O=C[C@H](O)[C@H]1OC[C@@H](O)[C@@H]1O", "3,6-anhydro-D-galactose", id="anhydro_galactose"),
    ],
)
def test_intramolecular_anhydro_sugars(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(C)C[C@H](N)C(=O)N[C@H](CCC(=O)O)C(=O)N[C@@H]([C@H](C)O)C(=O)N[C@H](C(C)C)C(=O)N[C@@H](CC(C)C)C(=O)O", "L-leucyl-D-glutamyl-L-allothreonyl-D-valyl-L-leucine", id="d_residues_and_allo"),
        pytest.param("NCC(=O)NCC(=O)O", "glycylglycine", id="glycine_residues_have_no_descriptor"),
        pytest.param("CC(N)C(=O)N[C@@H](C)C(=O)O", "ξ-alanyl-L-alanine", id="unspecified_residue_xi"),
        pytest.param("N[C@@H](CS)C(=O)N[C@@H](CC(=O)O)C(=O)N1CCC[C@H]1C(=O)O", "L-cysteinyl-L-aspartyl-L-proline", id="irregular_acyl_endings"),
    ],
)
def test_peptide_acyl_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("NCCC[C@H](N)C(=O)O", "L-ornithine", id="side_chain_table_entry"),
        pytest.param("N[C@@H](CS(=O)(=O)O)C(=O)O", "L-cysteic acid", id="sulfur_on_c3_reverses_cip_to_ld"),
        pytest.param("CC(=O)NCCC[C@H](N)C(=O)O", "N5-acetyl-L-ornithine", id="side_chain_nitrogen_locant"),
        pytest.param("NCCC(=O)O", "β-alanine", id="whole_molecule_retained_name"),
        pytest.param("N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O", "L-cystine", id="two_centres_one_descriptor"),
        pytest.param("N[C@@H](CSC[C@@H](N)C(=O)O)C(=O)O", "(2R,2'S)-3,3'-sulfanediylbis(2-aminopropanoic acid)", id="disagreeing_centres_fall_through"),
    ],
)
def test_less_common_amino_acid_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O[C@@H]1CCN[C@@H]1C(=O)O", "(3R)-3-hydroxy-L-proline", id="ring_centre_by_cip_alpha_by_ld"),
        pytest.param("OC1CNC(C1)C(=O)O", "4-hydroxyproline", id="no_stereo_no_hyphen"),
    ],
)
def test_substituted_proline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-]C(=O)CC[C@H](N)C(=O)O", "L-glutamate(1–)", id="diacid_monoanion_charge_cited"),
        pytest.param("NCCCC[C@H]([NH3+])C(=O)O", "L-lysinium(1+)", id="two_amino_groups_charge_cited"),
        pytest.param("NC(=[NH2+])NCCC[C@H](N)C(=O)O", "L-argininium(1+)", id="side_chain_cation_neutralized"),
        pytest.param("NCCCC[C@H]([NH3+])C(=O)O.[Cl-]", "L-lysinium(1+) chloride", id="cation_with_anion"),
        pytest.param("[Na+].[O-]C(=O)CC[C@H](N)C(=O)O", "sodium L-glutamate", id="counter_ion_fixes_the_charge"),
    ],
)
def test_ionized_amino_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("CC(=O)OC[C@H]1O[C@H](Br)[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@@H]1OC(C)=O", "2,3,4,6-tetra-O-acetyl-α-D-glucopyranosyl bromide", id="glycosyl_halide_esters_as_acyl_prefixes"),
        pytest.param("OC[C@H]1O[C@@H](OC)[C@H](OC)[C@@H](O)[C@@H]1O", "methyl 2-O-methyl-β-D-glucopyranoside", id="substituted_glycoside"),
        pytest.param("CC(=O)OC[C@H]1O[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@@H]1OC(C)=O", "β-D-glucopyranose 1,2,3,4,6-pentaacetate", id="anomeric_ester_like_other_esters"),
        pytest.param("OC[C@H]1O[C@@H](N)[C@H](OC(C)=O)[C@@H](O)[C@@H]1O", "2-O-acetyl-β-D-glucopyranosylamine", id="substituted_glycosylamine"),
    ],
)
def test_substituted_glycosides_and_glycosyl_derivatives(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("O=C1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "D-glucono-1,5-lactone", id="six_membered_ring_ester"),
        pytest.param("O=C1O[C@H]([C@H](O)CO)[C@H](O)[C@H]1O", "D-glucono-1,4-lactone", id="five_membered_ring_ester"),
    ],
)
def test_aldonic_acid_lactones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected

