import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem structure match: "propane-2-selone" (same structure as
        # the Blue Book's own "propane-2-thione (PIN)" for the sulfur case).
        ("CC(=[Se])C", "propane-2-selone"),
        # PubChem structure match: "pentane-2,4-diselone".
        ("CC(=[Se])CC(=[Se])C", "pentane-2,4-diselone"),
        # PubChem structure match: "cyclohexaneselone".
        ("C1CCC(=[Se])CC1", "cyclohexaneselone"),
        # PubChem structure match: "1-chloropropane-2-selone".
        ("ClCC(=[Se])C", "1-chloropropane-2-selone"),
        # Monocyclic ring, single ring double bond (P-31.1.3): the selone
        # always gets locant 1 (suffix priority), the ring double bond's
        # locant is minimized by choosing direction. Cross-checked against
        # PubChem (CID 102393135/102393138).
        ("[Se]=C1CCCC=C1", "cyclohex-2-ene-1-selone"),
        ("[Se]=C1CC=CCC1", "cyclohex-3-ene-1-selone"),
    ],
)
def test_selone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_selone_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCCC=C1C")


def test_unsaturated_ring_selone_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCCC#C1")


def test_ring_substituent_chain_selone():
    # A selone entirely on a chain hanging off a plain saturated ring (the
    # ring itself bears no selone) -- mirrors `_thione.py`'s
    # `test_ring_substituent_chain_thione`. No PubChem-registered
    # structure for this exact molecule (sparse selone coverage), so this
    # is a structural/regression check on the mechanism ported verbatim
    # from `_ketone.py`/`_thione.py`.
    assert smiles_to_iupac("CC(=[Se])C1CCCCC1") == "1-cyclohexylethane-1-selone"


def test_ring_substituent_chain_selone_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=[Se])C1CCC(C)CC1")


def test_ring_substituent_chain_selone_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=[Se])C1CCCC=C1")


def test_ring_with_selone_chain_selone_tie():
    # Ring and chain each carry exactly one selone (P-44.1.1 tie), mirrors
    # `_thione.py`'s `test_ring_with_thione_chain_thione_tie`. PubChem
    # PUG REST-confirmed "selanylidene" substituent-prefix naming
    # ("4-selanylidenepentan-2-one" CID 175299406) applied here as a
    # structural/regression check (no exact-structure PubChem record).
    assert (
        smiles_to_iupac("[Se]=C1CCCCC1C(=[Se])C")
        == "2-(1-selanylideneethyl)cyclohexane-1-selone"
    )


def test_selenoaldehyde_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=[Se]")


def test_polycyclic_selone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Se]=C1CCC2(CCCCC2)CC1")


def test_selone_with_hydroxyl_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=[Se])C")


def test_selenol_not_confused_with_selone():
    assert smiles_to_iupac("C[SeH]") == "methaneselenol"


def test_selenide_not_confused_with_selone():
    assert smiles_to_iupac("C[Se]C") == "methylselanylmethane"


def test_diselenide_not_confused_with_selone():
    assert smiles_to_iupac("C[Se][Se]C") == "(methyldiselanyl)methane"


def test_thione_not_confused_with_selone():
    assert smiles_to_iupac("CC(=S)C") == "propane-2-thione"


def test_ketone_not_confused_with_selone():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_acyclic_selone_stereocenter():
    # A single specified tetrahedral stereocenter (P-92): a selone's C=Se
    # carbon is double-bonded to selenium exactly like a ketone's C=O
    # carbon, so this mirrors `_thione.py`/`_ketone.py` cleanly (CIP
    # computed entirely by RDKit's `rdCIPLabeler`). PubChem has no
    # registered stereoisomer for this molecule, so this is a
    # structural/regression check on the already-proven mechanism, not an
    # independent PubChem cross-check.
    assert smiles_to_iupac("CC[C@@H](C)C(C)=[Se]") == "(3R)-3-methylpentane-2-selone"


def test_cyclic_selone_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_ketone.py`'s
    # `_name_cyclic_ketone` -- the branch sits at C4 (para to the C=Se
    # carbon, which -- like a ketone's carbonyl carbon -- is sp2 and
    # cannot itself bear a branch). PubChem has no registered structure
    # here (selenium compounds are sparse), so this is a
    # structural/regression check on the already-proven mechanism ported
    # from `_ketone.py`.
    assert (
        smiles_to_iupac("[Se]=C1CCC(CC1)[C@@H](C)CC")
        == "4-[(2S)-butan-2-yl]cyclohexane-1-selone"
    )


def test_selone_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(C)=[Se]") == "3-methylpentane-2-selone"


def test_phenyl_chain_selone():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#283's
    # carboxylic-acid/.../thione chains): the ring is cited as a "phenyl"
    # substituent prefix, mirroring `_thione.py`'s
    # '1-phenylpropane-2-thione'.
    assert smiles_to_iupac("c1ccccc1CC(=[Se])C") == "1-phenylpropane-2-selone"


def test_phenyl_chain_selone_longer_chain():
    assert smiles_to_iupac("c1ccccc1CC(=[Se])CC") == "1-phenylbutane-2-selone"


def test_phenyl_directly_attached_selone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=[Se])C")


def test_phenyl_substituted_benzene_ring_selone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=[Se])C")
