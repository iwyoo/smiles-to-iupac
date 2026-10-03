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


def test_unsaturated_ring_tellone():
    # Monocyclic ring, single ring double bond (P-31.1.3): same pattern
    # already confirmed for `_thione.py`/`_selone.py`. The exact ring
    # structures aren't PubChem-registered (sparse tellone ring coverage,
    # same reason `test_cyclohexanetellone` above has no PubChem name
    # citation), but the identical ring shape is PubChem-confirmed for the
    # selenium analogue (`cyclohex-2-ene-1-selone`, CID 102393135).
    assert smiles_to_iupac("[Te]=C1CCCC=C1") == "cyclohex-2-ene-1-tellone"
    assert smiles_to_iupac("[Te]=C1CC=CCC1") == "cyclohex-3-ene-1-tellone"


def test_unsaturated_ring_tellone_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC=C1C")


def test_unsaturated_ring_tellone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C1CCCC#C1")


def test_ring_substituent_chain_tellone():
    # A tellone entirely on a chain hanging off a plain saturated ring
    # (the ring itself bears no tellone) -- mirrors `_thione.py`'s
    # `test_ring_substituent_chain_thione`. No PubChem-registered
    # structure for this exact molecule (sparse tellone coverage), so
    # this is a structural/regression check on the mechanism ported
    # verbatim from `_ketone.py`/`_thione.py`.
    assert smiles_to_iupac("CC(=[Te])C1CCCCC1") == "1-cyclohexylethanetellone"


def test_ring_substituent_chain_tellone_ring_with_substituent():
    assert smiles_to_iupac("CC(=[Te])C1CCC(C)CC1") == "1-(4-methylcyclohexyl)ethanetellone"


def test_ring_substituent_chain_tellone_unsaturated_ring():
    assert smiles_to_iupac("CC(=[Te])C1CCCC=C1") == "1-(cyclohex-2-en-1-yl)ethanetellone"


def test_ring_with_tellone_chain_tellone_tie():
    # Ring and chain each carry exactly one tellone (P-44.1.1 tie),
    # mirrors `_thione.py`'s `test_ring_with_thione_chain_thione_tie`.
    # No PubChem-registered structure for the exact molecule or even the
    # "tellanylidene" substituent prefix alone (sparse tellurium
    # coverage), so this is a structural/regression check on the
    # mechanism ported verbatim from `_thione.py`/`_selone.py`.
    assert (
        smiles_to_iupac("[Te]=C1CCCCC1C(=[Te])C")
        == "2-(1-tellanylideneethyl)cyclohexane-1-tellone"
    )


def test_telluroaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Te]")


def test_von_baeyer_spiro_tellone_now_supported():
    # A single tellone on a von Baeyer bicyclic/polycyclic or monospiro
    # ring system is supported (see `test_von_baeyer_spiro_tellone.py`);
    # this SMILES is a spiro shape.
    assert smiles_to_iupac("[Te]=C1CCC2(CCCCC2)CC1") == "spiro[5.5]undecane-3-tellone"


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


def test_cyclic_tellone_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_ketone.py`'s
    # `_name_cyclic_ketone` -- the branch sits at C4 (para to the C=Te
    # carbon, which -- like a ketone's carbonyl carbon -- is sp2 and
    # cannot itself bear a branch). PubChem has no registered structure
    # here (tellurium compounds are sparse), so this is a
    # structural/regression check on the already-proven mechanism ported
    # from `_ketone.py`.
    assert (
        smiles_to_iupac("[Te]=C1CCC(CC1)[C@@H](C)CC")
        == "4-[(2S)-butan-2-yl]cyclohexane-1-tellone"
    )


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
