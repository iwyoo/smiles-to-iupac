import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CC=C2C(=C1)C=CS2=O", "1-benzothiophene 1-oxide"),  # PubChem CID 5383918
        ("c1ccc2cs(=O)cc2c1", "2-benzothiophene 2-oxide"),
    ],
)
def test_fused_hetero_ring_oxide_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_fused_ring_unaffected():
    assert smiles_to_iupac("c1ccc2sccc2c1") == "1-benzothiophene"


def test_nitrogen_fused_ring_still_raises():
    # Quinoline/isoquinoline/indole have a nitrogen, not a chalcogen -- an
    # oxide there is `_amine_oxide.py`'s territory (a charged N+/O- shape),
    # not this module's.
    with pytest.raises(Exception):
        smiles_to_iupac("O=[n+]1ccc2ccccc2c1")
