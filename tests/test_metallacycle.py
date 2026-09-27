import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-69.4's own worked skeletal-replacement pattern: the metal
        # always sits at locant 1, cited explicitly (unlike the analogous
        # single-heteroatom Hantzsch-Widman case).
        ("[Ni]1CCCC1", "1-nickelacyclopentane"),
        ("[Ti]1CCC1", "1-titanacyclobutane"),
        ("[Fe]1CCCC1", "1-ferracyclopentane"),
        ("[Pt]1CCC1", "1-platinacyclobutane"),
    ],
)
def test_metallacycle_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metallacycle_with_ligand_raises():
    # A ligand on the metal (P-69.2's coordination-nomenclature naming) is
    # out of scope for this bare-parent-hydride pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ni]1(Cl)CCCC1")


def test_metallacycle_with_two_metals_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ni]1CC[Ni]C1")
