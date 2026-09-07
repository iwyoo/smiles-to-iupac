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


def test_ammonium_ethanoate():
    # PubChem's own computed IUPACName for this exact SMILES is "azanium
    # acetate" -- 'ethanoate' is this project's own systematic stem (see
    # test_sodium_acetate above), and 'azanium' (not 'ammonium') is
    # `_ammonium.py`'s own PIN for unsubstituted NH4+.
    assert smiles_to_iupac("[NH4+].CC(=O)[O-]") == "azanium ethanoate"


def test_methanaminium_ethanoate():
    # PubChem's own computed IUPACName for this exact SMILES is
    # "methylazanium acetate"; `_ammonium.py`'s own PIN convention for a
    # substituted ammonium is the '-aminium' derivation, not the
    # 'azanium' substitutive style (see that module's docstring), giving
    # "methanaminium" here instead of PubChem's "methylazanium".
    assert smiles_to_iupac("C[NH3+].CC(=O)[O-]") == "methanaminium ethanoate"


def test_quaternary_ammonium_cation_ethanoate():
    assert (
        smiles_to_iupac("C[N+](C)(C)C.CC(=O)[O-]")
        == "N,N,N-trimethylmethanaminium ethanoate"
    )


def test_two_ammonium_cation_fragments_raises():
    # more than one cation fragment (a multi-cation salt) is unattempted
    # here, same as the existing monoatomic-cation scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH4+].[NH4+].CC(=O)[O-].CC(=O)[O-]")


def test_non_carboxylate_anion_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Na+].CCO")


def test_sodium_methoxide():
    # PubChem's own computed IUPACName for this exact SMILES is "sodium
    # methanolate"; `_alkoxide.py`'s own PIN convention retains the
    # traditional 'methoxide' name over the systematic 'methanolate' one
    # (P-72.2.2.2.2), so this module's output is "sodium methoxide".
    assert smiles_to_iupac("[Na+].C[O-]") == "sodium methoxide"


def test_potassium_ethoxide():
    assert smiles_to_iupac("[K+].CC[O-]") == "potassium ethoxide"


def test_ammonium_methoxide():
    # PubChem's own computed IUPACName for this exact SMILES is "azanium
    # methanolate" -- same 'methoxide' vs 'methanolate' divergence as
    # test_sodium_methoxide above.
    assert smiles_to_iupac("[NH4+].C[O-]") == "azanium methoxide"


def test_divalent_cation_alkoxide_anion_raises():
    # a 2+/3+ metal cation paired with an alkoxide anion is out of scope:
    # PubChem's own generator inconsistently omits the multiplying prefix
    # here (`[Ca+2].C[O-].C[O-]` -> "calcium methanolate", not "calcium
    # dimethanolate", unlike the carboxylate case's confirmed "calcium
    # diacetate"), so there is no reliable worked example to verify this
    # form against yet.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Ca+2].C[O-].C[O-]")
