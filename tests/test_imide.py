import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=CNC=O", "N-methanoylmethanamide"),
        ("CC(=O)NC(=O)C", "N-ethanoylethanamide"),
        ("CCC(=O)NC(=O)CC", "N-propanoylpropanamide"),
        ("CCCC(=O)NC(=O)CCC", "N-butanoylbutanamide"),
    ],
)
def test_symmetric_imide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_symmetric_imide():
    assert smiles_to_iupac("CC(C)C(=O)NC(=O)C(C)C") == "N-(2-methylpropanoyl)-2-methylpropanamide"


def test_halogen_substituent():
    assert smiles_to_iupac("ClCC(=O)NC(=O)CCl") == "N-(2-chloroethanoyl)-2-chloroethanamide"


def test_unsymmetric_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC(=O)CC")


def test_n_substituted_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(C)C(=O)C")


def test_unsaturated_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)NC(=O)C=C")


def test_cyclic_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCC(=O)N1")
