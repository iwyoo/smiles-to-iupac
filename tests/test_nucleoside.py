import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N", "adenosine"),  # PubChem CID 60961
        ("C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)NC(=NC2=O)N", "guanosine"),  # CID 6802
        ("C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O", "inosine"),  # CID 6021
        ("C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)NC(=O)NC2=O", "xanthosine"),  # CID 64959
        ("C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O", "cytidine"),  # CID 6175
        ("CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CO)O", "thymidine"),  # CID 5789
        ("C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O", "uridine"),  # CID 6029
    ],
)
def test_nucleoside_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nucleoside_base_attached_via_oxygen_is_not_matched():
    with pytest.raises(Exception):
        smiles_to_iupac("C1=CN(C(=O)NC1=O)O[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O")


def test_nucleoside_wrong_stereoisomer_is_not_matched():
    with pytest.raises(Exception):
        smiles_to_iupac("C1=CN(C(=O)NC1=O)[C@@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O")
