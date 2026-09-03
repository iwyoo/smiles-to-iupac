import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_propane_2_tellone():
    # PubChem structure match: "propane-2-tellone" (same structure as the
    # Blue Book's own "propane-2-thione (PIN)" for the sulfur case).
    assert smiles_to_iupac("CC(=[Te])C") == "propane-2-tellone"


def test_pentane_2_4_ditellone():
    # `_thione.py`'s/`_selone.py`'s own precedent generalizes group-count
    # support without a specific PubChem worked example for every
    # chalcogen (PubChem has no computed IUPACName for this exact
    # structure) -- the shared locant/suffix machinery already handles it.
    assert smiles_to_iupac("CC(=[Te])CC(=[Te])C") == "pentane-2,4-ditellone"


def test_cyclohexanetellone():
    assert smiles_to_iupac("C1CCC(=[Te])CC1") == "cyclohexanetellone"


def test_telluroaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Te]")


def test_polycyclic_tellone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCC2(CCCCC2)CC1")


def test_tellurol_not_confused_with_tellone():
    assert smiles_to_iupac("C[TeH]") == "methanetellurol"


def test_telluride_not_confused_with_tellone():
    assert smiles_to_iupac("C[Te]C") == "methyltellanylmethane"


def test_selone_not_confused_with_tellone():
    assert smiles_to_iupac("CC(=[Se])C") == "propane-2-selone"


def test_thione_not_confused_with_tellone():
    assert smiles_to_iupac("CC(=S)C") == "propane-2-thione"


def test_acyclic_tellone_stereocenter():
    # A single specified tetrahedral stereocenter (P-92): a tellone's
    # C=Te carbon is double-bonded to tellurium exactly like a ketone's
    # C=O carbon, so this mirrors `_thione.py`/`_ketone.py` cleanly (CIP
    # computed entirely by RDKit's `rdCIPLabeler`). PubChem has no
    # registered stereoisomer for this molecule, so this is a
    # structural/regression check on the already-proven mechanism, not an
    # independent PubChem cross-check.
    assert smiles_to_iupac("CC[C@@H](C)C(C)=[Te]") == "(3R)-3-methylpentane-2-tellone"


def test_tellone_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(C)=[Te]") == "3-methylpentane-2-tellone"


def test_phenyl_chain_tellone():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#284's
    # carboxylic-acid/.../selone chains): the ring is cited as a "phenyl"
    # substituent prefix, mirroring `_thione.py`'s
    # '1-phenylpropane-2-thione'.
    assert smiles_to_iupac("c1ccccc1CC(=[Te])C") == "1-phenylpropane-2-tellone"


def test_phenyl_chain_tellone_longer_chain():
    assert smiles_to_iupac("c1ccccc1CC(=[Te])CC") == "1-phenylbutane-2-tellone"


def test_phenyl_directly_attached_tellone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=[Te])C")


def test_phenyl_substituted_benzene_ring_tellone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=[Te])C")
