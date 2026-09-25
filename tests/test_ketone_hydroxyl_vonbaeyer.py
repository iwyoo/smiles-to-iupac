import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 13974583/85551307 -- the two bicyclo[2.2.1]heptanone
        # bridgehead-adjacent hydroxyl positions.
        ("C1C2CC(C1CC2=O)O", "5-hydroxybicyclo[2.2.1]heptan-2-one"),
        ("C1CC2C(C1CC2=O)O", "7-hydroxybicyclo[2.2.1]heptan-2-one"),
        # PubChem CID 115011371/130037946 -- bicyclo[3.2.1]octanone.
        ("C1CC2CC(=O)CC1C2O", "8-hydroxybicyclo[3.2.1]octan-3-one"),
        ("C1CC2CC1CC(C2=O)O", "3-hydroxybicyclo[3.2.1]octan-2-one"),
        # PubChem CID 440017, hydroxycamphor -- exercises alkyl substituents
        # coexisting with the hydroxyl.
        ("CC1(C2CC(=O)C1(CC2O)C)C", "5-hydroxy-1,7,7-trimethylbicyclo[2.2.1]heptan-2-one"),
        # PubChem CID 10964884/64184 -- tricyclic adamantanone, exercises
        # the `polycyclic_core` branch rather than `bicyclic_core`.
        ("C1C2CC3CC1CC(C2)(C3=O)O", "1-hydroxytricyclo[3.3.1.1^3,7]decan-2-one"),
        ("C1C2CC3CC(C2)(CC1C3=O)O", "5-hydroxytricyclo[3.3.1.1^3,7]decan-2-one"),
    ],
)
def test_hydroxy_von_baeyer_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_von_baeyer_ketone_without_hydroxyl_still_resolves():
    assert smiles_to_iupac("O=C1CC2CCC1CC2") == "bicyclo[2.2.2]octan-2-one"


def test_unsaturated_ring_alongside_hydroxyl_still_raises():
    # Testosterone's own ring unsaturation (androst-4-ene) is a separate,
    # still out-of-scope gate -- a coexisting hydroxyl doesn't bypass it.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34")
