import pytest

from chemonym import smiles_to_iupac
from chemonym._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-14.3.4.2(d): the locant '1' is omitted for unsubstituted two- and
        # three-carbon alkenes/alkynes.
        ("C=C", "ethene"),
        ("C#C", "acetylene"),  # P-31.1.2.1: 'acetylene' is the retained PIN, not 'ethyne'.
        ("CC=C", "propene"),
        ("CC#C", "propyne"),
        # Four carbons and up: the locant is always essential, since more
        # than one position is possible.
        ("CC=CC", "but-2-ene"),
        ("C=CCC", "but-1-ene"),
        ("CC#CC", "but-2-yne"),
        ("C#CCC", "but-1-yne"),
        ("C#CCCC", "pent-1-yne"),
        ("CCCC=C", "pent-1-ene"),
        ("CCC#CC", "pent-2-yne"),
        # P-14.3.3 / P-14.3.4.2(d): once a substituent is present, even a
        # three-carbon chain must cite the otherwise-omittable locant.
        ("C=C(C)C", "2-methylprop-1-ene"),
        ("C=C(C)CC", "2-methylbut-1-ene"),
        ("CC(C)=CC", "2-methylbut-2-ene"),
        ("CC(C)C=C", "3-methylbut-1-ene"),
        # P-14.4(e): lowest locant for the double bond outranks lowest locant
        # for substituents when choosing numbering direction. Numbering left
        # to right gives the double bond locant 2 and the methyl locant 6;
        # numbering right to left would give the methyl locant 2 but pushes
        # the double bond to locant 5. The double bond wins the tie-break.
        ("CC=CCCC(C)C", "6-methylhept-2-ene"),
    ],
)
def test_smiles_to_iupac_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C=CC=C",  # buta-1,3-diene: two double bonds, out of scope.
        "C=CC#C",  # a double bond and a triple bond together, out of scope.
        "C#CC#C",  # two triple bonds, out of scope.
        "C1=CCCCC1",  # cyclohexene: unsaturation in a ring, out of scope.
    ],
)
def test_out_of_scope_unsaturation_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_unsaturated_compound_substituent_raises():
    # Same decane-plus-sec-butyl-like-branch shape as
    # test_acyclic.test_compound_substituent_raises, but with a double bond
    # placed in the main chain: the branch still forks (P-29.4), which this
    # module doesn't support any more than `_acyclic.py` does.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCCC(C(C)CC)CCCCC")
