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
