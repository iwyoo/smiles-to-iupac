import pytest

from chemonym import smiles_to_iupac
from chemonym._common import UnsupportedStructure


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
        "c1ccccc1",
        # C=C and C#C now name as ethene/acetylene (see test_unsaturated.py);
        # two multiple bonds together remain out of scope here.
        "C=CC=C",
        "CCO",
    ],
)
def test_out_of_scope_structures_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


def test_compound_substituent_raises():
    # A decane chain (the unique longest chain) carrying a sec-butyl-like branch
    # at C5: the branch itself forks, which needs P-29.4 compound-substituent
    # naming, not yet implemented. The branch can't be absorbed into a longer
    # main chain because the decane backbone is strictly longer than any path
    # running through the branch instead.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC(C(C)CC)CCCCC")
