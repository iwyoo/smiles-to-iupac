import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_nitrone_n_methyl():
    # Real PubChem structure, CID 59937459.
    assert smiles_to_iupac("C[CH]=[N+](C)[O-]") == "N-methylethanimine N-oxide"


def test_nitrone_unsubstituted_nitrogen():
    assert smiles_to_iupac("C=[N+]([H])[O-]") == "methanimine N-oxide"


def test_nitrile_oxide_acetonitrile():
    # Real PubChem structure, CID 522391 -- the retained PIN 'acetonitrile'
    # (P-66.5.1.2.1), not `_nitrile.py`'s own systematic 'ethanenitrile'.
    assert smiles_to_iupac("CC#[N+][O-]") == "acetonitrile oxide"


def test_nitrile_oxide_propane():
    assert smiles_to_iupac("CCC#[N+][O-]") == "propanenitrile oxide"


def test_plain_imine_not_confused_with_nitrone():
    assert smiles_to_iupac("CC=NC") == "N-methylethanimine"


def test_plain_nitrile_not_confused_with_nitrile_oxide():
    assert smiles_to_iupac("CCC#N") == "propanenitrile"


def test_nitrile_imide_still_out_of_scope():
    # Nitrile imide (P-74.2.2.2.1.1) needs the zwitterionic hydrazinium-ide
    # method instead of this module's simple 'oxide' functional-class
    # naming -- deferred, not this step.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC#[N+][N-]C")
