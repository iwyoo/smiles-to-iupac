import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_oxidanium():
    # PubChem structure match: "oxidanium" (OH3+).
    assert smiles_to_iupac("[OH3+]") == "oxidanium"


def test_methyloxidanium():
    # PubChem structure match: "methyloxidanium".
    assert smiles_to_iupac("C[OH2+]") == "methyloxidanium"


def test_ethyloxidanium():
    # PubChem structure match: "ethyloxidanium".
    assert smiles_to_iupac("CC[OH2+]") == "ethyloxidanium"


def test_dimethyloxidanium():
    # PubChem structure match: "dimethyloxidanium".
    assert smiles_to_iupac("C[OH+]C") == "dimethyloxidanium"


def test_trimethyloxidanium():
    # PubChem structure match: "trimethyloxidanium".
    assert smiles_to_iupac("C[O+](C)C") == "trimethyloxidanium"


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[OH2+]")


def test_ring_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[OH2+]C1CCCCC1")


def test_alcohol_not_confused_with_oxonium():
    assert smiles_to_iupac("CO") == "methanol"


def test_ether_not_confused_with_oxonium():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_sulfonium_not_confused_with_oxonium():
    assert smiles_to_iupac("C[SH2+]") == "methylsulfanium"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem PUG REST confirmed: `c1ccccc1[OH2+]` ->
        # "phenyloxidanium" (CID 5152889),
        # `c1ccccc1[O+](c1ccccc1)c1ccccc1` -> "triphenyloxidanium" (CID
        # 3474029). A charged oxygen's cation valence is 3, so three
        # phenyls never trigger a lambda-convention label.
        ("c1ccccc1[OH2+]", "phenyloxidanium"),
        ("c1ccccc1[O+](c1ccccc1)c1ccccc1", "triphenyloxidanium"),
    ],
)
def test_phenyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_and_phenyl_substituents():
    assert smiles_to_iupac("C[O+](C)c1ccccc1") == "dimethyl(phenyl)oxidanium"


def test_substituted_phenyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[OH2+]")
