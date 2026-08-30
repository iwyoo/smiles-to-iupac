import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_ditelluride():
    # PubChem PUG REST CID 88493, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Te][Te]C") == "(methylditellanyl)methane"


def test_diethyl_ditelluride():
    # PubChem PUG REST CID 141264, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Te][Te]CC") == "(ethylditellanyl)ethane"


def test_methyl_ethyl_ditelluride():
    # A 2-carbon parent omits the locant even though the sole substituent
    # is compound (parenthesized), same rule as `_diselenide.py`. PubChem
    # PUG REST CID 86011476, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Te][Te]CC") == "(methylditellanyl)ethane"


def test_methyl_propyl_ditelluride():
    # A 3-carbon parent needs the locant -- no PubChem-listed compound
    # found for this specific structure (CID 0), so this is a reviewed
    # result, not an independently verified one: the mechanism itself
    # already has independent confirmation via `_diselenide.py`'s own
    # identical '1-(methyldiselanyl)propane' case.
    assert smiles_to_iupac("CCC[Te][Te]C") == "1-(methylditellanyl)propane"


def test_branched_ditellanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Te][Te]C(C)C")


def test_tritellurium_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Te][Te][Te]C")


def test_stereocenter_on_parent_chain():
    # No PubChem-registered stereoisomer for this shape (specified SMILES
    # resolves to CID 0); cross-checked against RDKit's independent
    # rdCIPLabeler, which agrees with the R/S the parent chain's own
    # stereocenter should carry -- same evidentiary bar as `_peroxide.py`'s
    # PR #225.
    assert smiles_to_iupac("CC[C@H](C)[Te][Te]CC") == "(2S)-2-(ethylditellanyl)butane"


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CCC(C)[Te][Te]CC") == "2-(ethylditellanyl)butane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC[Te][Te][C@H](C)CC")
