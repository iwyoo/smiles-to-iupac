import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_disulfide():
    # PubChem PUG REST CID 12232, auto-generated name matches exactly --
    # note the enclosing parens (unlike the plain 'methylsulfanylmethane'
    # of the corresponding sulfide), same as `_diselenide.py`.
    assert smiles_to_iupac("CSSC") == "(methyldisulfanyl)methane"


def test_diethyl_disulfide():
    # PubChem PUG REST CID 8077, auto-generated name matches exactly.
    assert smiles_to_iupac("CCSSCC") == "(ethyldisulfanyl)ethane"


def test_methyl_ethyl_disulfide():
    # A 2-carbon parent omits the locant even though the sole substituent
    # is compound (parenthesized). PubChem PUG REST CID 123388,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CSSCC") == "(methyldisulfanyl)ethane"


def test_methyl_propyl_disulfide():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 16592,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCCSSC") == "1-(methyldisulfanyl)propane"


def test_branched_disulfanyl_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)SSC(C)C")


def test_trisulfur_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CSSSC")


def test_terminal_persulfide_methane():
    # PubChem PUG REST CID 522059 -- a bare 'disulfanyl' with no locant
    # next to it needs no P-16.3.3 parentheses, unlike the alkyl-prefixed
    # 'methyldisulfanyl' cases above.
    assert smiles_to_iupac("CSS") == "disulfanylmethane"


def test_terminal_persulfide_ethane():
    # PubChem PUG REST CID 94671.
    assert smiles_to_iupac("CCSS") == "disulfanylethane"


def test_terminal_persulfide_propane():
    # PubChem PUG REST CID 6428842 -- a 3-carbon parent needs a locant, and
    # the parentheses return once a digit sits directly in front of the
    # bare 'disulfanyl' name.
    assert smiles_to_iupac("CCCSS") == "1-(disulfanyl)propane"


def test_both_terminal_disulfane_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SS")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No PubChem-registered stereoisomer for this shape (specified
        # SMILES resolves to CID 0); cross-checked against RDKit's
        # independent rdCIPLabeler, which agrees with the R/S the parent
        # chain's own stereocenter should carry -- same evidentiary bar as
        # `_peroxide.py`'s PR #225.
        ("CC[C@H](C)SSCC", "(2S)-2-(ethyldisulfanyl)butane"),
        ("CC[C@@H](C)SSCC", "(2R)-2-(ethyldisulfanyl)butane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CCC(C)SSCC") == "2-(ethyldisulfanyl)butane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCSS[C@H](C)CC")


def test_phenyl_disulfide_direct_bond():
    # A plain, unsubstituted benzene ring directly bonded to one sulfur
    # (P-44.1.2.2 rule (1), same "ring always wins" pattern as
    # `_azide.py`'s benzene-ring path -- 'disulfanyl' has no suffix form).
    # PubChem CID 84234.
    assert smiles_to_iupac("c1ccccc1SSC") == "(methyldisulfanyl)benzene"
    # PubChem CID 257711.
    assert smiles_to_iupac("c1ccccc1SSCC") == "(ethyldisulfanyl)benzene"


def test_phenyl_disulfide_chain_spacer_not_supported():
    # A -CH2- spacer between the ring and the near sulfur is deliberately
    # out of scope for this narrow first slice (see module docstring) --
    # PubChem CID 12779 ('(methyldisulfanyl)methylbenzene') nests the
    # 'disulfanyl' prefix in a way this project's general-purpose
    # `name_branch` machinery doesn't reproduce without further research.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CSSC")


def test_phenyl_disulfide_perthiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1SS")


def test_phenyl_disulfide_branched_other_side_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1SSC(C)C")


def test_phenyl_disulfide_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1SSC")


def test_phenyl_disulfide_unsaturation_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CSSC")
