import pytest

from smiles_to_iupac import smiles_to_iupac

ANDROSTANE = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2{x})CC[C@@H]4[C@@]3(CC{a}C4)C"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-101.7.1.1.2: substitution is preferred to a 'nor' parent.
        (ANDROSTANE.format(x="C", a="C"), "17β-methyl-5α-androstane"),
        # P-101.7.1.1.3: the two 17-methyls carry no α/β.
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2(C)C)CC[C@@H]4[C@@]3(CCCC4)C", "17,17-dimethyl-5α-androstane"),
        # P-101.7.1.2: the acid and ester suffixes outrank the ketone.
        (ANDROSTANE.format(x="C(=O)O", a="C(=O)"), "3-oxo-5α-androstane-17β-carboxylic acid"),
        (ANDROSTANE.format(x="C(=O)OC", a="C(=O)"), "methyl 3-oxo-5α-androstane-17β-carboxylate"),
        (ANDROSTANE.format(x="C(=O)N", a="C"), "5α-androstane-17β-carboxamide"),
        # P-66.6.1.1.3, P-66.5.1.1.2: a ring -CHO or -CN is a carbaldehyde or carbonitrile suffix, otherwise a formyl or cyano prefix.
        (ANDROSTANE.format(x="C=O", a="C(=O)"), "3-oxo-5α-androstane-17β-carbaldehyde"),
        (ANDROSTANE.format(x="C#N", a="C(C=O)"), "3ξ-formyl-5α-androstane-17β-carbonitrile"),
        ("O=Cc1ccc2c(c1)CC[C@@H]1[C@@H]2CC[C@]2(C)C(=O)CC[C@@H]12", "17-oxoestra-1,3,5(10)-triene-3-carbaldehyde"),
        # P-101.7.3: an ester of a steroid alcohol names the steroid as a substituent group.
        (ANDROSTANE.format(x="C(C)=O", a="[C@@H](OC(C)=O)"), "20-oxo-5α-pregnan-3α-yl acetate"),
        (ANDROSTANE.format(x="OC(C)=O", a="C"), "5α-androstan-17β-yl acetate"),
    ],
)
def test_derivatives_on_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-101.7.1.3: a carbon already in the terminal segment is the acid carbon.
        ("OC(=O)C12CCCC1C1CCC3=CC(=O)CCC3(C)C1CC2", "3-oxoandrost-4-en-18-oic acid"),
        ("OC(=O)CC1CCC2C1(C)CCC1C2CCC2=CC(=O)CCC12C", "3-oxopregn-4-en-21-oic acid"),
    ],
)
def test_terminal_segment_acid_is_an_oic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-101.7.3, P-67.1.3.2: the sulfate ester outranks the ketone, which is cited in the steroid group.
        ("C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](OS(=O)(=O)O)CC[C@]34C)[C@@H]1CCC2=O", "17-oxoandrost-5-en-3\u03b2-yl hydrogen sulfate"),
        # P-101.7.3, P-66.1.1.4: the steroid is the N-substituent of an amide.
        ("C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](NC(C)=O)CC[C@]34C)[C@@H]1CC[C@@H]2O", "N-(17\u03b2-hydroxyandrost-5-en-3\u03b2-yl)acetamide"),
        # P-101.7.3, P-15.3: one acid on two alcohol groups of the steroid.
        ("C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](OC(C)=O)CC[C@]34C)[C@@H]1CC[C@@H]2OC(C)=O", "androst-5-ene-3\u03b2,17\u03b2-diyl diacetate"),
    ],
)
def test_steroid_as_substituent_group_of_esters_and_amides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-101.2.6.1.3: a side-chain centre is cited with its CIP descriptor.
        (
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2[C@H](C)O)CC[C@@H]4[C@@]3(CC[C@@H](N(C)C)C4)C",
            "(20S)-3α-(dimethylamino)-5α-pregnan-20-ol",
        ),
        (
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2[C@@H](C)O)CC[C@@H]4[C@@]3(CC[C@@H](N(C)C)C4)C",
            "(20R)-3α-(dimethylamino)-5α-pregnan-20-ol",
        ),
    ],
)
def test_side_chain_cip_descriptor(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "geometry, expected",
    [("/C=C/", "(23E)-5α-cholest-23-ene"), ("/C=C\\", "(23Z)-5α-cholest-23-ene")],
)
def test_side_chain_double_bond_geometry(geometry, expected):
    # P-101.6.2
    smiles = f"C[C@H](C{geometry}C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C"
    assert smiles_to_iupac(smiles) == expected




ANDROSTANE_EXO = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC{c}2{x})CC[C@@H]4[C@@]3(CC{a}C4)C"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-64.7.1: the ketone outranks the ylidene group, which is cited as a prefix.
        (
            "C[C@]12CC[C@H]3[C@@H](CC(=C)C4=CC(=O)C=C[C@]34C)[C@@H]1CCC2=O",
            "6-methylideneandrosta-1,4-diene-3,17-dione",
        ),
        (ANDROSTANE_EXO.format(c="C", x="=C", a="C"), "17-methylidene-5α-androstane"),
        # P-62.3.1.1: an oxime is an N-hydroxy imine named with the suffix 'imine'.
        (ANDROSTANE_EXO.format(c="C", x="=NO", a="C"), "N-hydroxy-5α-androstan-17-imine"),
        (ANDROSTANE_EXO.format(c="C", x="=N", a="C"), "5α-androstan-17-imine"),
        # P-62.3.1.2, P-64.7.1: below a ketone the imine is cited as a prefix.
        (ANDROSTANE_EXO.format(c="C", x="=NO", a="C(=O)"), "17-(hydroxyimino)-5α-androstan-3-one"),
        (ANDROSTANE_EXO.format(c="C", x="=NN", a="C(=O)"), "17-hydrazinylidene-5α-androstan-3-one"),
        # P-62.3.1.1: several imines take N locants numbered by the position of their carbon.
        (ANDROSTANE_EXO.format(c="C", x="=NO", a="C(=NOC)"), "N17-hydroxy-N3-methoxy-5α-androstane-3,17-diimine"),
        # P-64.6.1, P-64.6.2: the chalcogen analogues of the ketone are suffixes below it.
        (ANDROSTANE_EXO.format(c="C", x="=S", a="C(=S)"), "5α-androstane-3,17-dithione"),
        (ANDROSTANE_EXO.format(c="C", x="=S", a="C(=O)"), "17-sulfanylidene-5α-androstan-3-one"),
        # P-93.5.1.4.2: the geometry of the double bond to an ylidene group or an imine nitrogen is cited with the ring locant.
        (ANDROSTANE_EXO.format(c="/C", x="=N/O", a="C"), "(17Z)-N-hydroxy-5α-androstan-17-imine"),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1C/C(=C/c1ccccc1)C2=O)CC[C@@H]4[C@@]3(CCCC4)C", "(16Z)-16-benzylidene-5α-androstan-17-one"),
    ],
)
def test_exocyclic_double_bonds_on_a_steroid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "x",
    ["=NC(=O)C", "=NOC(C)=O", "=NN", "=NC#N"],
)
def test_imine_nitrogen_bearing_a_senior_group_is_not_an_imine_suffix(x):
    name = smiles_to_iupac(ANDROSTANE_EXO.format(c="C", x=x, a="C"))
    assert not name.endswith("-17-imine")


def test_stereo_centre_inside_a_steroid_ester_alkyl_is_cited():
    smiles = "C[C@H](CC)OC(=O)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    assert smiles_to_iupac(smiles) == "(2R)-butan-2-yl 3β-hydroxyandrost-5-ene-17β-carboxylate"


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles",
    [
        "NC(=O)c1ccc([C@@H]2O[C@@H]3C[C@H]4[C@@H]5CCC6=CC(=O)C=C[C@]6(C)[C@H]5[C@@H](O)C[C@]4(C)[C@]3(C(=O)CO)O2)cc1",
        "OC(=O)c1ccc([C@@H]2O[C@@H]3C[C@H]4[C@@H]5CCC6=CC(=O)C=C[C@]6(C)[C@H]5[C@@H](O)C[C@]4(C)[C@]3(C(=O)CO)O2)cc1",
    ],
)
def test_cyclic_acetal_is_not_named_with_class_word_after_a_senior_carbonyl_component(smiles):
    assert "acetal" not in smiles_to_iupac(smiles)


def test_o_acyl_groups_of_different_acids_are_cited_with_the_locants_of_the_diyl_groups():
    smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@H]4C[C@@H](OC(C)=O)CC[C@]34C)[C@@H]1CC[C@@H]2OC(=O)CC"
    assert smiles_to_iupac(smiles) == "5α-androstane-3β,17β-diyl 3-acetate 17-propanoate"


def test_esters_of_different_alkyl_groups_are_cited_alphanumerically_with_locants():
    smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@H]4C[C@@H](C(=O)OC)CC[C@]34C)[C@@H]1CC[C@@H]2C(=O)OCC"
    assert smiles_to_iupac(smiles) == "17-ethyl 3-methyl 5α-androstane-3β,17β-dicarboxylate"


def test_terminal_ester_below_the_principal_group_is_cited_as_alkoxy_and_oxo():
    smiles = "C[C@]12CC[C@H]3[C@@H](CC[C@H]4C[C@@H](C(=O)O)CC[C@]34C)[C@@H]1CC[C@@H]2CC(=O)OC"
    assert smiles_to_iupac(smiles) == "21-methoxy-21-oxo-5α-pregnane-3β-carboxylic acid"
