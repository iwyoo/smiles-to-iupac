import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1cn2nccnc2n1", "imidazo[1,2-b][1,2,4]triazine"),
        ("OC(=O)CN1CC=NC2=NC=CN12", "2-[imidazo[1,2-b][1,2,4]triazin-1(2H)-yl]ethanoic acid"),
        ("OC(=O)CC1=CN2N=CC=NC2=N1", "2-(imidazo[1,2-b][1,2,4]triazin-6-yl)ethanoic acid"),
        ("OC(=O)Cc1cn2ccccc2n1", "2-(imidazo[1,2-a]pyridin-2-yl)ethanoic acid"),
        ("OC(=O)Cc1ccn2ccnc2c1", "2-(imidazo[1,2-a]pyridin-7-yl)ethanoic acid"),
        ("OC(=O)CN1C=CN2C=CC=CC12", "2-[imidazo[1,2-a]pyridin-1(8aH)-yl]ethanoic acid"),
    ],
)
def test_bridgehead_heteroatom_fused_yl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
