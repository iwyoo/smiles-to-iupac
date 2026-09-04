import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_diselenide():
    # PubChem PUG REST CID 23496, auto-generated name matches exactly --
    # note the enclosing parens (unlike the plain 'methylselanylmethane'
    # of the corresponding selenide), see module docstring.
    assert smiles_to_iupac("C[Se][Se]C") == "(methyldiselanyl)methane"


def test_diethyl_diselenide():
    # PubChem PUG REST CID 69405, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Se][Se]CC") == "(ethyldiselanyl)ethane"


def test_methyl_propyl_diselenide():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 85591054,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[Se][Se]C") == "1-(methyldiselanyl)propane"


def test_methyl_ethyl_diselenide():
    # A 2-carbon parent omits the locant even though the sole substituent
    # is compound (parenthesized) -- see module docstring. PubChem PUG
    # REST CID 129678518, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Se][Se]CC") == "(methyldiselanyl)ethane"


def test_branched_diselanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[Se][Se]C(C)C")


def test_triselenium_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[Se][Se][Se]C")


def test_terminal_perselenol_methane():
    # PubChem PUG REST CID 101729611 -- a bare 'diselanyl' with no locant
    # next to it needs no P-16.3.3 parentheses, unlike the alkyl-prefixed
    # 'methyldiselanyl' cases above.
    assert smiles_to_iupac("C[Se][SeH]") == "diselanylmethane"


def test_terminal_perselenol_ethane():
    # PubChem PUG REST CID 173348765.
    assert smiles_to_iupac("CC[Se][SeH]") == "diselanylethane"


def test_terminal_perselenol_propane():
    # PubChem PUG REST CID 174964680 -- a 3-carbon parent needs a locant,
    # and the parentheses return once a digit sits directly in front of
    # the bare 'diselanyl' name.
    assert smiles_to_iupac("CCC[Se][SeH]") == "1-(diselanyl)propane"


def test_both_terminal_diselane_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH][SeH]")


def test_stereocenter_on_parent_chain():
    # No PubChem-registered stereoisomer for this shape (specified SMILES
    # resolves to CID 0); cross-checked against RDKit's independent
    # rdCIPLabeler, which agrees with the R/S the parent chain's own
    # stereocenter should carry -- same evidentiary bar as `_peroxide.py`'s
    # PR #225.
    assert smiles_to_iupac("CC[C@H](C)[Se][Se]CC") == "(2S)-2-(ethyldiselanyl)butane"


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CCC(C)[Se][Se]CC") == "2-(ethyldiselanyl)butane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC[Se][Se][C@H](C)CC")


def test_phenyl_diselenide_direct_bond():
    # A plain, unsubstituted benzene ring directly bonded to one
    # selenium (P-44.1.2.2 rule (1), same "ring always wins" pattern as
    # `_disulfide.py`'s benzene-ring path -- 'diselanyl' has no suffix
    # form). PubChem CID 59041517.
    assert smiles_to_iupac("c1ccccc1[Se][Se]C") == "(methyldiselanyl)benzene"
    # PubChem CID 71327844.
    assert smiles_to_iupac("c1ccccc1[Se][Se]CC") == "(ethyldiselanyl)benzene"


def test_phenyl_diselenide_chain_spacer_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C[Se][Se]C")


def test_phenyl_diselenide_seh_terminal_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Se][SeH]")


def test_phenyl_diselenide_branched_other_side_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Se][Se]C(C)C")


def test_phenyl_diselenide_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[Se][Se]C")
