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
        # Alcohols/amines now name via name_alcohol/name_amine instead (see
        # test_alcohol.py: "CCO" -> "ethanol", test_amine.py: "CCN" ->
        # "ethanamine"), and a secondary/tertiary amine is supported too
        # (see test_amine.py: "CCNCC" -> "N-ethylethanamine"); a
        # secondary/tertiary amine nitrogen on a ring still is out of scope.
        "CN(C)C1CCCCC1",
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A halogenated chain's own stereocenter (P-91.3/P-92), same
        # mechanism as `_carboxylic_acid.py`'s/`_ether.py`'s stereocenter
        # tests (CIP computed entirely by RDKit's `rdCIPLabeler`, not
        # reimplemented here).
        ("C[C@H](Cl)CC", "(2S)-2-chlorobutane"),
        ("C[C@@H](Cl)CC", "(2R)-2-chlorobutane"),
    ],
)
def test_acyclic_alkane_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_alkane_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention (see `_common.py`'s `specified_stereocenters` docstring).
    assert smiles_to_iupac("CC(Cl)CC") == "2-chlorobutane"
