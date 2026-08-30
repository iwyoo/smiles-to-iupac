import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_formaldehyde_azine_name():
    # PubChem PUG REST CID 136329, auto-generated name matches exactly
    # (queried via POST -- GET 400s on a SMILES containing '/').
    assert smiles_to_iupac("C=NN=C") == "N-(methylideneamino)methanimine"


def test_acetaldehyde_azine_name():
    # PubChem PUG REST CID 68964, auto-generated name matches exactly.
    assert smiles_to_iupac("CC=NN=CC") == "N-(ethylideneamino)ethanimine"


def test_propionaldehyde_azine_name():
    # The imine (parent) side keeps `_imine.py`'s own "cite the locant even
    # at C1 for a 3+-carbon chain" rule. PubChem PUG REST CID 123355,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC=NN=CCC") == "N-(propylideneamino)propan-1-imine"


def test_acetone_azine_name():
    # A ketone-shaped attachment (locant 2) on both sides. PubChem PUG
    # REST CID 79085, auto-generated name matches exactly.
    assert smiles_to_iupac("CC(C)=NN=C(C)C") == "N-(propan-2-ylideneamino)propan-2-imine"


def test_asymmetric_azine_raises():
    # The two C=N substituents differ (acetaldehyde vs. acetone side) --
    # which side becomes the imine parent is unverified, out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NN=C(C)C")


def test_n_substituted_hydrazone_still_raises():
    # Sanity check: a plain (non-azine) N-substituted hydrazone is still
    # out of scope and must not be misrouted into this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NNC")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1=NN=C1CCCCC1")


def test_aromatic_carbon_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C=NN=Cc1ccccc1")


def test_azine_specified_chain_stereocenter_raises():
    # Both of this azine's C=N bonds are always unspecified potential
    # Bond_Double stereo elements to RDKit, regardless of substituents
    # (same conclusion as `_imine.py`/`_hydrazone.py`) -- so a specified
    # chain stereocenter here always coexists with at least one of them,
    # and `specified_stereocenters` correctly rejects the combination
    # (P-92/P-93) instead of the silent drop this project's
    # stereodescriptor safety net exists to fix.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(C)=NN=C(C)[C@H](C)CC")
