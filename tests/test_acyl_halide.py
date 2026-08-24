import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # -oyl halide is a mechanical stem+halide-word construction (P-65.1.5,
        # see _acyl_halide.py docstring) with a single, non-ambiguous
        # substitution point; no locant tie-break or alphabetization choice
        # is involved, so these are cross-checked directly against common
        # usage (acetyl chloride == 'ethanoyl chloride' etc.) rather than a
        # per-case PubChem lookup.
        ("C(=O)Cl", "methanoyl chloride"),
        ("CC(=O)Cl", "ethanoyl chloride"),
        ("CCC(=O)Cl", "propanoyl chloride"),
        ("CCCC(=O)F", "butanoyl fluoride"),
        ("CCCC(=O)Br", "butanoyl bromide"),
        ("CCCC(=O)I", "butanoyl iodide"),
    ],
)
def test_saturated_acyl_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_acyl_halide():
    assert smiles_to_iupac("C=CC(=O)Cl") == "prop-2-enoyl chloride"


def test_other_halogen_substituent_coexists():
    assert smiles_to_iupac("CC(Cl)C(=O)Cl") == "2-chloropropanoyl chloride"


def test_branched_substituent():
    assert smiles_to_iupac("CC(C)CC(=O)Cl") == "3-methylbutanoyl chloride"


def test_two_acyl_halides_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC(=O)CC(=O)Cl")


def test_cyclic_acyl_halide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC(=O)C1CCCCC1")


def test_coexisting_ketone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)CC(=O)Cl")
