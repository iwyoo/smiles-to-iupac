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
    ],
)
def test_selone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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
