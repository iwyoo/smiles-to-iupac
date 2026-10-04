import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=c1cccc[nH]1", "pyridin-2(1H)-one"),  # PubChem CID 8871
        ("O=c1cc[nH]cc1", "pyridin-4(1H)-one"),  # PubChem CID 12290
    ],
)
def test_pyridinone_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_pyridine_unaffected():
    assert smiles_to_iupac("c1ccncc1") == "pyridine"


def test_non_aromatic_ketone_unaffected():
    assert smiles_to_iupac("O=C1C=CC=CC1") == "cyclohexa-2,4-dien-1-one"
