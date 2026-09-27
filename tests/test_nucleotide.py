import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 6083's isomeric SMILES for AMP confirms this exact
        # structure.
        ("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N", "5'-adenylic acid"),
        ("C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N", "5'-guanylic acid"),
        ("C1=NC2=C(C(=O)N1)N=CN2[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O", "5'-inosinic acid"),
        # Xanthylic acid is esterified at 3', not 5', unlike the other six.
        ("C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O)O)NC(=O)NC2=O", "3'-xanthylic acid"),
        ("C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O", "5'-cytidylic acid"),
        ("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)COP(=O)(O)O)O", "5'-thymidylic acid"),
        ("C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O", "5'-uridylic acid"),
    ],
)
def test_nucleotide_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_nucleoside_unaffected():
    assert smiles_to_iupac("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N") == "adenosine"


def test_nucleoside_phosphorylated_at_wrong_position_raises():
    # AMP's phosphate esterifying the 3'-OH instead of the prescribed 5'
    # position is a different, unrecognized structure.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O)O)N")


def test_nucleoside_diphosphate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N"
        )
