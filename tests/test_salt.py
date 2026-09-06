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


def test_calcium_diethanoate():
    # Blue Book P-65.6.2.1's own worked example is "calcium diacetate
    # (PIN)"; this project's systematic 'ethanoate' stem (see
    # test_sodium_acetate above) makes "calcium diethanoate" the expected
    # name here.
    assert smiles_to_iupac("[Ca+2].CC(=O)[O-].CC(=O)[O-]") == "calcium diethanoate"


def test_aluminium_tripropanoate():
    assert (
        smiles_to_iupac("[Al+3].CCC(=O)[O-].CCC(=O)[O-].CCC(=O)[O-]")
        == "aluminium tripropanoate"
    )


def test_unbalanced_divalent_cation_raises():
    # only one anion fragment for a 2+ cation is not charge-balanced.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ca+2].CC(=O)[O-]")


def test_mixed_anions_on_divalent_cation_raises():
    # two different carboxylate anions on the same cation (a mixed salt)
    # is out of scope -- both anions must be identical here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ca+2].CC(=O)[O-].CCC(=O)[O-]")


def test_transition_metal_cation_raises():
    # transition metals need Stock/oxidation-number disambiguation
    # (e.g. "iron(II)" vs "iron(III)"), not attempted here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Fe+2].CC(=O)[O-].CC(=O)[O-]")


def test_polyatomic_cation_raises():
    # ammonium (a molecular, not monoatomic, cation) needs its own name
    # assembled via `_ammonium.py` -- combining that into a salt name is
    # unattempted here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH4+].CC(=O)[O-]")


def test_non_carboxylate_anion_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Na+].CCO")
