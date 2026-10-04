import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._amino_acid import _SIDE_CHAIN_TABLE
from smiles_to_iupac._common import UnsupportedStructure


def test_side_chain_table_has_no_collisions():
    assert len(_SIDE_CHAIN_TABLE) == 15
    assert len(set(_SIDE_CHAIN_TABLE.values())) == 15


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
        ("C([C@H]([C@@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O)O", "beta-D-glucofuranose"),  # CID 11309871
    ],
)
def test_cyclic_aldofuranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_furanose_cites_the_specified_elements():
    assert (
        smiles_to_iupac("C([C@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O")
        == "(2R,3R,4R,5S)-5-(hydroxymethyl)oxolane-2,3,4-triol"
    )


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1[C@H]([C@H]([C@H]([C@](O1)(CO)O)O)O)O", "beta-D-psicopyranose"),
    ],
)
def test_cyclic_ketohexopyranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aldopyranose_still_resolves():
    assert (
        smiles_to_iupac("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O") == "alpha-D-glucopyranose"
    )


def test_ketofuranose_cites_the_specified_elements():
    assert (
        smiles_to_iupac("C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O")
        == "(2R,3S,4S,5R)-2,5-bis(hydroxymethyl)oxolane-2,3,4-triol"
    )


def test_d_histidine():
    assert smiles_to_iupac("C1=C(NC=N1)C[C@H](C(=O)O)N") == "D-histidine"


def test_histidine_unspecified_stereocenter_no_ld_prefix():
    assert smiles_to_iupac("C1=C(NC=N1)CC(C(=O)O)N") == "histidine"


def test_ring_substituted_histidine_is_named_as_an_amino_acid_with_a_ring_prefix():
    assert (
        smiles_to_iupac("Cn1cnc(CC(N)C(=O)O)c1")
        == "2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid"
    )


def test_tryptophan_still_resolves():
    assert smiles_to_iupac("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O") == "L-tryptophan"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1O", "neo-inositol"),  # (1,2,3/4,5,6-)
    ],
)
def test_inositol_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nucleoside_base_attached_via_oxygen_is_not_matched():
    with pytest.raises(Exception):
        smiles_to_iupac("C1=CN(C(=O)NC1=O)O[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O")


def test_nucleoside_wrong_stereoisomer_is_not_matched():
    with pytest.raises(Exception):
        smiles_to_iupac("C1=CN(C(=O)NC1=O)[C@@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O",
            "5'-inosinic acid",
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
            "NC1=NC(=O)N(C=C1C#N)[C@H]2C[C@@H]([C@H](O2)COC(C)=O)O",
            "4-amino-1-(5-O-acetyl-2-deoxy-β-D-erythro-pentofuranosyl)-2-oxo-1,2-dihydropyrimidine-5-carbonitrile",
        ),
    ],
)
def test_substituted_nucleoside_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_nucleoside_wrong_sugar_stereo_is_not_matched():
    with pytest.raises(Exception):
        smiles_to_iupac("CC1=CN(C(=O)NC1=O)[C@@H]2[C@H]([C@@H]([C@H](O2)CO)O)O")


def test_plain_nucleoside_unaffected():
    assert (
        smiles_to_iupac("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N")
        == "adenosine"
    )


def test_nucleoside_phosphorylated_at_wrong_position_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O)O)N")


def test_nucleoside_triphosphate_chain_raises():
    # A chain longer than triphosphate is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)"
            "COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
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


def test_unspecified_stereo_aldose_still_falls_through_unchanged():
    assert smiles_to_iupac("C(C(C(C(C(C=O)O)O)O)O)O") == "2,3,4,5,6-pentahydroxyhexanal"


def test_d_proline():
    assert smiles_to_iupac("C1C[C@@H](NC1)C(=O)O") == "D-proline"


def test_proline_unspecified_stereocenter_no_ld_prefix():
    assert smiles_to_iupac("C1CC(NC1)C(=O)O") == "proline"


def test_hydroxyproline_cites_the_specified_center():
    assert smiles_to_iupac("OC1C[C@H](NC1)C(=O)O") == "(2S)-4-hydroxypyrrolidine-2-carboxylic acid"


def test_hetero_ring_ketone_still_resolves():
    assert smiles_to_iupac("O=C1CCCCO1") == "oxan-2-one"
