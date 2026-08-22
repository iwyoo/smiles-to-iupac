import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C", "methane"),
        ("CC", "ethane"),
        ("CCC", "propane"),
        ("CCCC", "butane"),
        ("CCCCC", "pentane"),
        ("CCCCCC", "hexane"),
        ("CC(C)C", "2-methylpropane"),
        ("CC(C)(C)C", "2,2-dimethylpropane"),
        ("CCC(C)CC", "3-methylpentane"),
        ("CC(C)C(CC)CCC", "3-ethyl-2-methylhexane"),
        ("CC(C)CCCCC(C)C", "2,7-dimethyloctane"),
        # P-45.2.3 example: citation-order tiebreak among equal, equally-placed
        # locant sets.
        ("CCC(C)CC(CCCC)C(CCC)CC(CC)CC", "5-butyl-8-ethyl-3-methyl-6-propyldecane"),
    ],
)
def test_smiles_to_iupac(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        # Alcohols now name via name_alcohol instead (see test_alcohol.py:
        # "CCO" -> "ethanol"), so it's no longer out of scope here.
        "CCN",
    ],
)
def test_out_of_scope_structures_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_compound_substituent():
    # A decane chain (the unique longest chain) carrying a sec-butyl-like
    # branch at C5: the branch itself forks (root -> a methyl, and an ethyl
    # continuation), so it is named as a compound substituent (P-29.4), with
    # the longer (ethyl) continuation chosen as the branch's own chain and
    # the leftover methyl cited at its locant 1: '1-methylpropyl'. The branch
    # can't be absorbed into a longer main chain because the decane backbone
    # is strictly longer than any path running through the branch instead.
    assert smiles_to_iupac("CCCCC(C(C)CC)CCCCC") == "5-(1-methylpropyl)decane"
