import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[Ti](Cl)(Cl)Cl", "trichlorido(methyl)titanium"),
        ("CC[Ti](C)(Cl)Cl", "dichlorido(ethyl)(methyl)titanium"),
        ("C[Hg]C", "dimethylmercury"),
        ("C[Hg]Cl", "chlorido(methyl)mercury"),
        ("C[Zn]C", "dimethylzinc"),
        ("C[Pt](C)(C)C", "tetramethylplatinum"),
        ("c1ccccc1[Pd]Br", "bromido(phenyl)palladium"),
    ],
)
def test_coordination_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "C[Ti](Cl)(Cl)Cl.[Na+]",
        "C[Pt](C)(C)(C)[Pt](C)(C)C",
        "C[Ti](=O)(Cl)Cl",
        "C[Hg]c1ccc(C(=O)O)cc1",
    ],
)
def test_coordination_out_of_scope_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
