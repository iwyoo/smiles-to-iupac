import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,2-dihydronaphthalene: CAS 447-53-0, NIST WebBook and
        # Sigma-Aldrich both list "1,2-Dihydronaphthalene" as the name.
        ("C1CC=Cc2ccccc12", "1,2-dihydronaphthalene"),
        # 1,4-dihydronaphthalene: CAS 612-17-9, a distinct real compound
        # from the 1,2- isomer (NIST WebBook separate entry).
        ("C1C=CCc2ccccc12", "1,4-dihydronaphthalene"),
        # Same 1,4-dihydronaphthalene molecule as above, written starting
        # from the intact ring instead of the reduced ring -- the
        # lowest-locants tie-break (P-31.1.4.2) must still normalize this
        # to "1,4-", not "2,3-"/"5,8-"/"6,7-" (whichever the raw SMILES
        # atom order would otherwise suggest).
        ("c1ccc2c(c1)CC=CC2", "1,4-dihydronaphthalene"),
    ],
)
def test_smiles_to_iupac_dihydronaphthalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_naphthalene_itself_is_unaffected():
    # exactly zero hydro pairs still goes through _aromatic.py's own path,
    # not this module (find_dihydronaphthalene_core requires the reduced
    # ring to have exactly two sp3 atoms).
    assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


def test_tetrahydronaphthalene_raises():
    # more than one hydro pair (tetralin) is out of scope for this module
    # (see docstring) -- must not be mistaken for the single-pair case.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCc2ccccc12")


def test_substituted_dihydronaphthalene_raises():
    # any substituent is out of scope -- every ring atom must have exactly
    # its "bare" degree.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC=Cc2ccccc12")


def test_halogen_substituted_dihydronaphthalene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CC=Cc2ccccc12")
