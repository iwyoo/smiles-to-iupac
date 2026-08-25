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


@pytest.mark.parametrize(
    "smiles",
    [
        # 1,4-Epoxy-1,4-dihydronaphthalene ("7-oxabenzonorbornadiene"):
        # CAS 573-57-9, PubChem CID 97139 (structure/formula cross-checked
        # via ConnectivitySMILES "C1=CC=C2C3C=CC(C2=C1)O3" -- PubChem's own
        # computed IUPACName for this shape is a von Baeyer bridged-ring
        # name, not usable to verify the fusion+bridge name itself).
        "C1=CC2OC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring.
        "c1ccc2c(c1)C1C=CC2O1",
    ],
)
def test_bridged_naphthalene_epoxy(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-epoxynaphthalene"


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


def test_ether_itself_is_unaffected():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_substituted_epoxy_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2O1")


@pytest.mark.parametrize(
    "smiles",
    [
        # 9,10-Dihydro-9,10-methanoanthracene: literature name confirmed in
        # J. Org. Chem. ("9,10-Dihydro-9,10-methanoanthracene and Its
        # Perhydro Derivatives"); structure/formula (C15H12) cross-checked
        # against PubChem CID 12651785's ConnectivitySMILES.
        "C12c3ccccc3C(c3ccccc31)C2",
        # Same molecule, atom order starting from one of the intact
        # aromatic rings instead of a bridgehead -- the lowest-locants
        # tie-break must still normalize this to "9,10-", not some other
        # numbering (anthracene's meso positions are always 9,10 by
        # definition, but this still exercises the candidate search).
        "c12ccccc1C1c3ccccc3C2C1",
    ],
)
def test_bridged_anthracene(smiles):
    assert smiles_to_iupac(smiles) == "9,10-dihydro-9,10-methanoanthracene"


def test_anthracene_itself_is_unaffected():
    assert smiles_to_iupac("c1ccc2cc3ccccc3cc2c1") == "anthracene"


def test_bridged_anthracene_epoxy_raises():
    # Unlike the naphthalene case, an 'epoxy' bridge on anthracene's 9,10
    # positions is deliberately out of scope here -- no independently
    # verifiable name was found while scoping this (see
    # `find_bridged_anthracene_core`'s docstring), only the structure
    # (PubChem CID cross-checked via ConnectivitySMILES), which isn't
    # enough under this project's test-writing policy to assert a name.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12c3ccccc3C(c3ccccc31)O2")


def test_substituted_bridged_anthracene_bridge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12c3ccccc3C(c3ccccc31)C2(C)")


def test_substituted_bridged_anthracene_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1c3ccccc3C2C1")
