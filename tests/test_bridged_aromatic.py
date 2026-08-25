import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles",
    [
        # benzonorbornadiene: CAS 4453-90-1, PubChem CID 97391, NIST WebBook
        # gives the IUPAC name "1,4-Dihydro-1,4-methanonaphthalene".
        "C1=CC2CC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring
        # instead of the bridged ring -- the lowest-locants tie-break must
        # still normalize this to "1,4-", not some other numbering.
        "c1ccc2c(c1)C1C=CC2C1",
    ],
)
def test_bridged_naphthalene(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-methanonaphthalene"


def test_naphthalene_itself_is_unaffected():
    assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


def test_dihydronaphthalene_itself_is_unaffected():
    assert smiles_to_iupac("C1CC=Cc2ccccc12") == "1,2-dihydronaphthalene"


def test_fully_saturated_bridge_raises():
    # the bridged reduced ring must still carry the "ene" double bond
    # (P-25.4's bridge doesn't itself imply further saturation) -- a fully
    # saturated bridgehead pair is a different, out-of-scope shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC2CC1c1ccccc12")


def test_substituted_bridge_atom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC2C(C)C1c1ccccc12")


def test_substituted_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2C1")
