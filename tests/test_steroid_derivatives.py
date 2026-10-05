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
        # P-101.7.3: an ester of a steroid alcohol names the steroid as a substituent group.
        (ANDROSTANE.format(x="C(C)=O", a="[C@@H](OC(C)=O)"), "(3R,5S,8R,9S,10S,13S,14S,17S)-17-acetyl-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-3-yl acetate"),
        (ANDROSTANE.format(x="OC(C)=O", a="C"), "(5R,8R,9S,10S,13S,14S,17S)-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl acetate"),
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
