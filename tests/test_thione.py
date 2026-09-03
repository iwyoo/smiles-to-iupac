import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-64.6.1 worked example: "propane-2-thione (PIN) (not
        # thioacetone)". PubChem structure match confirms.
        ("CC(=S)C", "propane-2-thione"),
        # Blue Book P-64.6.1 worked example: "butane-2-thione (PIN)".
        ("CCC(=S)C", "butane-2-thione"),
        ("CC(=S)CC", "butane-2-thione"),
        # PubChem structure match: "hexane-3-thione".
        ("CCCC(=S)CC", "hexane-3-thione"),
        # Blue Book P-64.6.1 worked example: "pentane-2,4-dithione (PIN)".
        ("CC(=S)CC(=S)C", "pentane-2,4-dithione"),
        # PubChem structure match: "cyclohexanethione".
        ("C1CCC(=S)CC1", "cyclohexanethione"),
        # Halogen coexistence, PubChem structure match:
        # "1-chloropropane-2-thione".
        ("ClCC(=S)C", "1-chloropropane-2-thione"),
        # Unsaturated chain, PubChem structure match: "pent-4-ene-2-thione".
        ("C=CCC(=S)C", "pent-4-ene-2-thione"),
    ],
)
def test_thione_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_thial_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=S")


def test_aromatic_thione_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc(cc1)C(=S)C")


def test_polycyclic_thione_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S=C1CCC2(CCCCC2)CC1")


def test_thione_with_hydroxyl_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=S)C")


def test_thione_with_ketone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC(=S)C")


def test_thiol_not_confused_with_thione():
    assert smiles_to_iupac("CS") == "methanethiol"


def test_sulfide_not_confused_with_thione():
    assert smiles_to_iupac("CSC") == "methylsulfanylmethane"


def test_disulfide_not_confused_with_thione():
    assert smiles_to_iupac("CSSC") == "(methyldisulfanyl)methane"


def test_ketone_not_confused_with_thione():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_acyclic_thione_stereocenter():
    # A single specified tetrahedral stereocenter (P-92): a thione's C=S
    # carbon is double-bonded to sulfur exactly like a ketone's C=O
    # carbon, so this mirrors `_ketone.py` cleanly (CIP computed entirely
    # by RDKit's `rdCIPLabeler`). PubChem has no registered stereoisomer
    # for this molecule (thiones are sparsely covered there), so this is
    # a structural/regression check on the already-proven mechanism, not
    # an independent PubChem cross-check.
    assert smiles_to_iupac("CC[C@@H](C)C(C)=S") == "(3R)-3-methylpentane-2-thione"


def test_cyclic_thione_stereocenter():
    # A stereocenter on the ring itself (P-92), same pattern as
    # `_ketone.py`'s `_name_cyclic_ketone`. Same PubChem-sparsity caveat.
    assert smiles_to_iupac("S=C1CCCC[C@H]1C") == "(2R)-2-methylcyclohexane-1-thione"


def test_thione_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(C)=S") == "3-methylpentane-2-thione"


def test_phenyl_chain_thione():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#282's
    # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/acyl-
    # halide/sulfonic-acid/thiol/sulfinic-acid/thioic-acid/selenoic-acid/
    # telluroic-acid chains): the ring is cited as a "phenyl" substituent
    # prefix, mirroring `_ketone.py`'s '1-phenylpropan-2-one'.
    assert smiles_to_iupac("c1ccccc1CC(=S)C") == "1-phenylpropane-2-thione"


def test_phenyl_chain_thione_longer_chain():
    assert smiles_to_iupac("c1ccccc1CC(=S)CC") == "1-phenylbutane-2-thione"


def test_phenyl_directly_attached_thione_raises():
    # An aryl thione (C=S directly on the ring) is out of scope for this
    # acyclic-chain-parent module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=S)C")


def test_phenyl_substituted_benzene_ring_thione_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=S)C")


def test_phenyl_chain_thione_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=S)C")
