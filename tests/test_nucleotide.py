import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 6083's isomeric SMILES for AMP confirms this exact
        # structure.
        ("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N", "5'-adenylic acid"),
        # PubChem CID 6804 (GMP).
        ("C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)NC(=NC2=O)N", "5'-guanylic acid"),
        # PubChem CID 8582 (IMP).
        ("C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O", "5'-inosinic acid"),
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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # ADP, PubChem CID 6022.
        (
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
            "adenosine 5'-(trihydrogen diphosphate)",
        ),
        # ATP, PubChem CID 5957.
        (
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
            "adenosine 5'-(tetrahydrogen triphosphate)",
        ),
        # P-106.2's own worked examples (tmp/bluebook/P10.txt lines
        # 4536-4553): uridine 5'-triphosphate, and xanthosine's own
        # 3'-diphosphate (esterified at 3', not 5', like its monophosphate).
        (
            "C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O",
            "uridine 5'-(tetrahydrogen triphosphate)",
        ),
        (
            "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)OP(=O)(O)O)O)NC(=O)NC2=O",
            "xanthosine 3'-(trihydrogen diphosphate)",
        ),
        # GDP, PubChem CID 8977.
        (
            "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)NC(=NC2=O)N",
            "guanosine 5'-(trihydrogen diphosphate)",
        ),
        # GTP, PubChem CID 6830.
        (
            "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)NC(=NC2=O)N",
            "guanosine 5'-(tetrahydrogen triphosphate)",
        ),
        # IDP.
        (
            "C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O",
            "inosine 5'-(trihydrogen diphosphate)",
        ),
        # ITP.
        (
            "C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O",
            "inosine 5'-(tetrahydrogen triphosphate)",
        ),
    ],
)
def test_nucleoside_di_triphosphate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nucleoside_triphosphate_chain_raises():
    # A chain longer than triphosphate is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(
            "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)"
            "COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N"
        )
