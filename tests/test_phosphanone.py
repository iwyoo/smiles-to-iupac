import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_single_substituent_no_lambda():
    # The Blue Book's own direct worked example is 'C6H5-P=O' ->
    # 'phenylphosphanone (PIN)' (aromatic, unverified here for a plain
    # alkyl chain -- see module docstring), so this alkyl case is a
    # reviewed, not independently confirmed, extension of the same
    # single-substituent/no-lambda mechanism.
    assert smiles_to_iupac("CP=O") == "methylphosphanone"
    assert smiles_to_iupac("CCP=O") == "ethylphosphanone"


def test_three_identical_substituents_needs_lambda5():
    # Mirrors the Blue Book's own direct worked example structure,
    # '(C6H5)3P=O' -> 'triphenyl-λ5-phosphanone (PIN)' [explicitly "not
    # oxotriphenyl-λ5-phosphane"; "triphenylphosphane oxide" is cited as
    # an acceptable but non-preferred functional-class name] -- reviewed
    # here for the plain alkyl analogue (unverified for methyl
    # specifically, but the same mechanism).
    assert smiles_to_iupac("CP(C)(C)=O") == "trimethyl-λ5-phosphanone"


def test_two_substituents_needs_lambda5_by_analogy():
    # Two explicit substituents plus P=O (RDKit fills the remaining
    # valence with one implicit H, reaching total bond order 5) is
    # treated the same way as the three-substituent case by direct
    # analogy -- no worked example of its own confirms this specific
    # count, so this is a reviewed, not independently verified, result.
    assert smiles_to_iupac("CP(C)=O") == "dimethyl-λ5-phosphanone"


def test_two_different_substituents_parenthesization():
    # P-16.5.1.3.1: the first (alphabetically) substituent is never
    # parenthesized, every subsequent one is -- same rule
    # `_phosphane.py`'s own mixed-substituent case already applies.
    assert smiles_to_iupac("CCP(C)=O") == "ethyl(methyl)-λ5-phosphanone"


def test_zero_substituents_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=P")


def test_aromatic_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=Pc1ccccc1")


def test_branched_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)P=O")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=P1CCCCC1")


def test_plain_phosphane_still_works():
    # Sanity check: a phosphine oxide's own phosphorus-bonded oxygen must
    # not be misrouted to `_phosphane.py`, and a plain phosphane (no
    # oxygen) must not be misrouted here.
    assert smiles_to_iupac("CP") == "methylphosphane"
    assert smiles_to_iupac("CP(C)C") == "trimethylphosphane"
