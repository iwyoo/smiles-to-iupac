import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@@]12CCC[C@H]1[C@@H]1CC[C@@H]3CCC[C@]3(C)[C@H]1CC2", "4-nor-5α-androstane"),
        ("C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3COCC[C@]3(C)[C@H]1CC2", "3-oxa-5α-androstane"),
        ("C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CCCN[C@]3(C)[C@H]1CC2", "1-aza-5α-androstane"),
        ("C1CCC2C(C1)CCC1C3CCCC3CCC21", "gonane"),
        ("C[C@@]12CCC[C@H]1[C@@H]1CCC34CCCCC3(C4)[C@H]1CC2", "5,19-cyclo-5ξ,10ξ-androstane"),
        ("CC1CCCCC1CCC1CCC[C@]2(C)CCC[C@@H]12", "9,10-seco-5ξ,8ξ,10ξ-androstane"),
    ],
)
def test_modified_steroid_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCC2(C)C)C(C)(C)CC1",
            "2,2′-dinor-β,β-carotene",
        ),
        (
            "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
            "17-methyl-7,8-didehydro-4,5α-epoxymorphinan-3,6α-diol",
        ),
    ],
)
def test_other_retained_natural_product_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
            "3α,7α,12α-trihydroxy-5β-cholan-24-oic acid",
        ),
        (
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)O",
            "17α-ethynylestra-1,3,5(10)-triene-3,17β-diol",
        ),
        (
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2N)CCC4=CC(=O)CC[C@]34C",
            "17β-aminoandrost-4-en-3-one",
        ),
    ],
)
def test_natural_steroids_with_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "CC(=O)[C@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@@]21C",
            "17-hydroxy-17α-pregn-4-ene-3,20-dione",
        ),
        ("C[C@]12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34", "17ξ-hydroxy-8ξ,9ξ,10ξ,14ξ-estr-4-en-3-one"),
    ],
)
def test_configuration_that_differs_from_the_parent_or_is_unspecified(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_partly_specified_parent_cites_unspecified_centers_as_xi():
    assert smiles_to_iupac("C1CCC2[C@H](C1)CCC3C2CCC4C3CCC4") == "5α,8ξ,9ξ,10ξ,13ξ,14ξ-gonane"


def test_extra_methyl_on_a_ring_is_a_substituent_of_the_retained_parent():
    assert smiles_to_iupac("CC12CCCC1(C)CCC1C2CCC2C1CCCC2") == "14-methylestrane"
    assert smiles_to_iupac("CC12CCCC1C3CC(C)C4CCCCC4(C3CC2)C") == "6-methylandrostane"
    assert smiles_to_iupac("CC12CCCC1C1CC(C)C3CCCCC3C1CC2") == "6-methylestrane"
    assert smiles_to_iupac("CCC1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C") == "15-methylpregnane"
    assert smiles_to_iupac("CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C") == "15-methylcholane"
    assert smiles_to_iupac("CC(C)CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C") == "15-methylcholestane"


def test_side_chain_that_fits_no_retained_parent_is_a_substituent_of_the_ring_parent():
    assert smiles_to_iupac("CC(C)C(C)CC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "17-(4,5-dimethylhexan-2-yl)androstane"


def test_androstane_other_diastereomer_gets_its_own_alpha_beta_citation():
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "5α,8α-androstane"
    )


def test_androst_4_ene():
    assert smiles_to_iupac("CC12CCCC1C1CCC3=CCCCC3(C)C1CC2") == "androst-4-ene"
