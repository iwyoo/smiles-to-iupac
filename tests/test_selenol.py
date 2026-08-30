import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methaneselenol():
    # PubChem PUG REST CID 440764, auto-generated name matches exactly.
    assert smiles_to_iupac("C[SeH]") == "methaneselenol"


def test_ethaneselenol():
    # PubChem PUG REST CID 5252527, auto-generated name matches exactly --
    # the same P-14.3.4.2(b) two-carbon locant omission as 'ethanethiol'.
    assert smiles_to_iupac("CC[SeH]") == "ethaneselenol"


def test_propane_1_selenol():
    # PubChem PUG REST CID 71373846, auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[SeH]") == "propane-1-selenol"


def test_2_methylpropane_1_selenol():
    # A branched chain. PubChem PUG REST CID 157368643, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CC(C)C[SeH]") == "2-methylpropane-1-selenol"


def test_prop_2_ene_1_selenol():
    # An unsaturated chain. PubChem PUG REST CID 15821407, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("C=CC[SeH]") == "prop-2-ene-1-selenol"


def test_2_chloroethane_1_selenol():
    # Halogen coexistence (P-35.2.1) -- no PubChem-listed compound found
    # for this specific structure (every halogenated selenol SMILES tried
    # came back as CID 0), so this is a reviewed result, not an
    # independently verified one: the mechanism itself already has
    # independent confirmation via `_thiol.py`'s own identical
    # '2-chloroethane-1-thiol' case (same locant-citation rule once a
    # substituent is present on the two-carbon chain).
    assert smiles_to_iupac("ClCC[SeH]") == "2-chloroethane-1-selenol"


def test_ethane_1_2_diselenol():
    # PubChem PUG REST auto-generated name matches exactly.
    assert smiles_to_iupac("[SeH]CC[SeH]") == "ethane-1,2-diselenol"


def test_propane_1_3_diselenol():
    # PubChem PUG REST auto-generated name matches exactly.
    assert smiles_to_iupac("[SeH]CCC[SeH]") == "propane-1,3-diselenol"


def test_propane_1_2_3_triselenol():
    # `_thiol.py`'s own precedent generalizes group-count support beyond
    # two without a specific 3+-group PubChem worked example of its own
    # (the shared `_alcohol.py`-style locant/suffix machinery already
    # handles an arbitrary-length locant list) -- same reasoning applied
    # here.
    assert smiles_to_iupac("[SeH]CC([SeH])C[SeH]") == "propane-1,2,3-triselenol"


def test_butane_1_2_4_triselenol():
    assert smiles_to_iupac("[SeH]CC([SeH])CC[SeH]") == "butane-1,2,4-triselenol"


def test_cyclohexaneselenol():
    # PubChem structure match: "cyclohexaneselenol".
    assert smiles_to_iupac("C1CCCCC1[SeH]") == "cyclohexaneselenol"


def test_cyclopentaneselenol():
    # PubChem structure match: "cyclopentaneselenol".
    assert smiles_to_iupac("C1CCCC1[SeH]") == "cyclopentaneselenol"


def test_substituted_ring_selenol_locant_cited():
    assert smiles_to_iupac("CC1CCCCC1[SeH]") == "2-methylcyclohexane-1-selenol"


def test_polycyclic_selenol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12(CCC(CC1)CC2)[SeH]")


def test_selenol_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C[SeH]")


def test_selenol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[SeH]")


def test_acyclic_selenol_stereocenter():
    # A single specified tetrahedral stereocenter (P-92): like a thiol's
    # -SH sulfur, a selenol's -SeH selenium is monovalent and can never
    # itself be a stereocenter, so this mirrors `_thiol.py` cleanly (CIP
    # computed entirely by RDKit's `rdCIPLabeler`). PubChem CID 175184280.
    assert smiles_to_iupac("C[C@@H](CC)C[SeH]") == "(2S)-2-methylbutane-1-selenol"


def test_cyclic_selenol_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_thiol.py`'s `_name_cyclic_thiol`. PubChem has no registered
    # stereoisomer for this exact ring (even the non-stereo parent isn't
    # registered), so this is a structural/regression check on the
    # already-proven mechanism ported verbatim from `_thiol.py`, not an
    # independent PubChem cross-check.
    assert smiles_to_iupac("[SeH][C@H]1CCCC[C@H]1C") == "(1S,2R)-2-methylcyclohexane-1-selenol"


def test_selenol_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)[SeH]") == "butane-2-selenol"


def test_selenol_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH][C@H]1CCCCC1Cl")
