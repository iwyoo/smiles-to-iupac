import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_sodium_acetate():
    # PubChem's own computed IUPACName for this exact SMILES is "sodium
    # acetate"; this project's `_carboxylate.py` already deliberately uses
    # the systematic 'ethanoate' stem over the retained 'acetate' one for
    # consistency with `_carboxylic_acid.py`'s identical 'ethanoic acid'
    # choice (see that module's docstring) -- this module inherits the
    # same divergence, so "sodium ethanoate" is the expected PIN here, not
    # a mismatch.
    assert smiles_to_iupac("[Na+].CC(=O)[O-]") == "sodium ethanoate"


def test_potassium_propanoate():
    assert smiles_to_iupac("[K+].CCC(=O)[O-]") == "potassium propanoate"


def test_lithium_methanoate():
    assert smiles_to_iupac("[Li+].[O-]C=O") == "lithium methanoate"


def test_divalent_cation_raises():
    # a 2+ cation needs charge-balancing stoichiometry (e.g. "calcium
    # diethanoate") -- explicitly out of scope for this first, narrowest
    # pass (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ca+2].CC(=O)[O-].CC(=O)[O-]")


def test_polyatomic_cation_raises():
    # ammonium (a molecular, not monoatomic, cation) needs its own name
    # assembled via `_ammonium.py` -- combining that into a salt name is
    # unattempted here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH4+].CC(=O)[O-]")


def test_non_carboxylate_anion_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Na+].CCO")
