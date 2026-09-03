import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_worked_example_3_oxohexanal():
    # CH3CH2CH2COCH2CHO -> 3-oxohexanal: this specific worked example is
    # cited and confirmed in _aldehyde_ketone.py's own module docstring
    # (P-41/Table 3.3's 'al' > 'one' demotion).
    assert smiles_to_iupac("CCCC(=O)CC=O") == "3-oxohexanal"


def test_simplest_case():
    # A single-axis variant of the worked example above (shorter chain,
    # same demotion mechanism) -- CH3COCH2CHO -> 3-oxobutanal.
    assert smiles_to_iupac("CC(=O)CC=O") == "3-oxobutanal"


def test_ketone_locant_from_other_end():
    # aldehyde carbon must be C1, so the ketone locant is counted from the
    # aldehyde end regardless of which end it's closer to.
    assert smiles_to_iupac("O=CCCC(=O)C") == "4-oxopentanal"


def test_two_ketones():
    assert smiles_to_iupac("CC(=O)CC(=O)CC=O") == "3,5-dioxohexanal"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(=O)C(Cl)C=O") == "2-chloro-3-oxobutanal"


def test_plain_aldehyde_still_works():
    assert smiles_to_iupac("CCC=O") == "propanal"


def test_plain_ketone_still_works():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_two_aldehydes_with_ketone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC(=O)CC=O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)CC=O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC1CCC(=O)C1")


def test_hydroxyl_coexistence_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=O)CC=O")


def test_phenyl_chain_aldehyde_ketone():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_aldehyde.py`'s
    # phenyl-chain path): the ring is cited as a "phenyl" substituent
    # prefix alongside the demoted ketone's "oxo" prefix. PubChem PUG
    # REST: "3-oxo-4-phenylbutanal"/"2-oxo-4-phenylbutanal".
    assert smiles_to_iupac("c1ccccc1CC(=O)CC=O") == "3-oxo-4-phenylbutanal"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C=O") == "2-oxo-4-phenylbutanal"


def test_phenyl_directly_attached_to_ketone_carbon():
    # The ring attaches directly to the ketone carbon, not the aldehyde --
    # still a valid chain-parent case (no separate ring-suffix
    # construction needed here, unlike a -CHO directly on the ring).
    # PubChem PUG REST: "3-oxo-3-phenylpropanal".
    assert smiles_to_iupac("c1ccccc1C(=O)CC=O") == "3-oxo-3-phenylpropanal"


def test_phenyl_substituted_benzene_ring_aldehyde_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)CC=O")


def test_phenyl_chain_aldehyde_ketone_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)CC=O")
