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


def test_terminal_pertellurol_ethane():
    # PubChem PUG REST CID 173350364 -- a bare 'ditellanyl' with no locant
    # next to it needs no P-16.3.3 parentheses, unlike the alkyl-prefixed
    # 'methylditellanyl' cases above.
    assert smiles_to_iupac("CC[Te][TeH]") == "ditellanylethane"


def test_terminal_pertellurol_methane():
    # No PubChem-listed compound for this exact structure (CID 0) -- a
    # reviewed extension of the confirmed ethane case above and of
    # `_disulfide.py`/`_diselenide.py`'s identical mononuclear rule.
    assert smiles_to_iupac("C[Te][TeH]") == "ditellanylmethane"


def test_terminal_pertellurol_propane():
    # No PubChem-listed compound for this exact structure (CID 0) -- a
    # reviewed extension of `_disulfide.py`/`_diselenide.py`'s identical
    # confirmed 3-carbon case.
    assert smiles_to_iupac("CCC[Te][TeH]") == "1-(ditellanyl)propane"


def test_both_terminal_ditellane_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH][TeH]")


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


def test_phenyl_ditelluride_direct_bond():
    # A plain, unsubstituted benzene ring directly bonded to one
    # tellurium (P-44.1.2.2 rule (1), same "ring always wins" pattern as
    # `_disulfide.py`'s benzene-ring path -- 'ditellanyl' has no suffix
    # form). PubChem CID 101099279.
    assert smiles_to_iupac("c1ccccc1[Te][Te]C") == "(methylditellanyl)benzene"


def test_phenyl_ditelluride_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Te][Te]C")


def test_phenyl_ditelluride_teh_terminal_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Te][TeH]")


def test_phenyl_ditelluride_branched_other_side_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Te][Te]C(C)C")


def test_phenyl_ditelluride_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[Te][Te]C")
