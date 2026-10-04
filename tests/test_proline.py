import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_l_proline():
    # PubChem CID 145742.
    assert smiles_to_iupac("C1C[C@H](NC1)C(=O)O") == "L-proline"


def test_d_proline():
    # PubChem CID 25247.
    assert smiles_to_iupac("C1C[C@@H](NC1)C(=O)O") == "D-proline"


def test_proline_unspecified_stereocenter_no_ld_prefix():
    assert smiles_to_iupac("C1CC(NC1)C(=O)O") == "proline"


def test_hydroxyproline_cites_the_specified_center():
    assert smiles_to_iupac("OC1C[C@H](NC1)C(=O)O") == "(2S)-4-hydroxypyrrolidine-2-carboxylic acid"


def test_pipecolic_acid_cites_the_specified_elements():
    # A 6-membered ring analogue (piperidine-2-carboxylic acid) is a
    # different ring size, out of scope here.
    assert smiles_to_iupac("C1CC[C@H](NC1)C(=O)O") == '(2S)-piperidine-2-carboxylic acid'


def test_hetero_ring_ketone_still_resolves():
    # A plain saturated single-heteroatom ring bearing an actual ring
    # ketone (not proline's shape) must still route to `_ketone.py`
    # unchanged.
    assert smiles_to_iupac("O=C1CCCCO1") == "oxan-2-one"


def test_pyrrolidinone_lactam_still_resolves():
    # A pyrrolidine ring with the carbonyl *in* the ring (a lactam, not a
    # coexisting exocyclic -COOH) is a different shape, still routed to
    # `_ketone.py`'s hetero-ring-ketone path unchanged.
    assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"
