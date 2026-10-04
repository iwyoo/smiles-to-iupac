import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C", "17β-hydroxyandrost-4-en-3-one"),
        (
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
            "cholest-5-en-3β-ol",
        ),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O", "estra-1,3,5(10)-triene-3,17β-diol"),
        (
            "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O",
            "11β,17,21-trihydroxypregn-4-ene-3,20-dione",
        ),
        ("CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C", "pregn-4-ene-3,20-dione"),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CC[C@@H]4[C@@]3(CCC(=O)C4)C", "17β-hydroxy-5α-androstan-3-one"),
        (
            "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
            "3α,7α,12α-trihydroxy-5β-cholan-24-oic acid",
        ),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CCC4=C3C=CC(=C4)O", "3-hydroxyestra-1,3,5(10)-trien-17-one"),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)O", "17α-ethynylestra-1,3,5(10)-triene-3,17β-diol"),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)OC", "17α-ethynyl-3-methoxyestra-1,3,5(10)-trien-17β-ol"),
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2N)CCC4=CC(=O)CC[C@]34C", "17β-aminoandrost-4-en-3-one"),
    ],
)
def test_natural_steroids_with_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@H]2O)CC[C@@H]4[C@@]3(CCC(=O)C4)C", "17α-hydroxy-5α-androstan-3-one"),
        ("C[C@]12CC[C@@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C", "17β-hydroxy-9β-androst-4-en-3-one"),
        ("CC(=O)[C@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@@]21C", "17-hydroxy-17α-pregn-4-ene-3,20-dione"),
        ("C[C@]12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34", "17ξ-hydroxy-8ξ,9ξ,10ξ,14ξ-estr-4-en-3-one"),
    ],
)
def test_configuration_that_differs_from_the_parent_or_is_unspecified(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34", "17-hydroxyestr-4-en-3-one"),
        ("CC12CCC3C(C1CCC2O)CCc1cc(O)ccc13", "estra-1,3,5(10)-triene-3,17-diol"),
    ],
)
def test_steroids_without_stereochemistry_carry_no_descriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
