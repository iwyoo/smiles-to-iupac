import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-sulfinic acid' mirrors '-sulfonic acid' (see
        # test_sulfonic_acid.py) with one fewer oxygen; same locant rules.
        # 'methanesulfinic acid'/'ethanesulfinic acid' are also
        # independently verifiable real compound names.
        ("CS(=O)O", "methanesulfinic acid"),
        ("CCS(=O)O", "ethanesulfinic acid"),
        ("CCCS(=O)O", "propane-1-sulfinic acid"),
        ("CC(S(=O)O)C", "propane-2-sulfinic acid"),
        ("CCCCS(=O)O", "butane-1-sulfinic acid"),
    ],
)
def test_saturated_sulfinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_sulfinic_acid():
    assert smiles_to_iupac("C=CCS(=O)O") == "prop-2-ene-1-sulfinic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)S(=O)O") == "1-chloroethanesulfinic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C(S(=O)O)C")


def test_two_sulfinic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)CS(=O)O")


def test_sulfonic_acid_not_confused_with_sulfinic():
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_cyclohexanesulfinic_acid():
    # PubChem CID 3302189.
    assert smiles_to_iupac("OS(=O)C1CCCCC1") == "cyclohexanesulfinic acid"


def test_2_methylcyclohexane_1_sulfinic_acid():
    # PubChem CID 67183560.
    assert smiles_to_iupac("OS(=O)C1CCCCC1C") == "2-methylcyclohexane-1-sulfinic acid"


def test_cyclopentanesulfinic_acid():
    # PubChem CID 14138538.
    assert smiles_to_iupac("OS(=O)C1CCCC1") == "cyclopentanesulfinic acid"


def test_2_chlorocyclohexane_1_sulfinic_acid():
    # PubChem CID 67182750.
    assert smiles_to_iupac("OS(=O)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfinic acid"


def test_polycyclic_sulfinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfinic_acid():
    # Monocyclic ring, single -SO2H, single ring double bond (P-31.1.3):
    # the sulfinic acid always gets locant 1 (suffix priority), the ring
    # double bond's locant is minimized by choosing direction.
    # Cross-checked against PubChem (CID 12616418/69649903).
    assert smiles_to_iupac("OS(=O)C1CCCC=C1") == "cyclohex-2-ene-1-sulfinic acid"
    assert smiles_to_iupac("OS(=O)C1CC=CCC1") == "cyclohex-3-ene-1-sulfinic acid"


def test_unsaturated_ring_sulfinic_acid_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CCCC=C1C")


def test_unsaturated_ring_sulfinic_acid_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)C1CCCC#C1")


def test_ring_substituent_chain_sulfinic_acid():
    # A sulfinic acid entirely on a chain hanging off a plain saturated
    # ring (the ring itself bears no sulfinic acid) -- mirrors
    # `_sulfonic_acid.py`'s `test_ring_substituent_chain_sulfonic_acid`.
    # PubChem PUG REST-verified "cyclohexylmethanesulfinic acid"
    # (CID 57427118).
    assert smiles_to_iupac("OS(=O)CC1CCCCC1") == "cyclohexylmethanesulfinic acid"


def test_ring_substituent_chain_sulfinic_acid_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)CC1CCC(C)CC1")


def test_sulfinic_acid_unspecified_stereocenter_unaffected():
    # A genuine chain stereocenter left unspecified (no @/@@) is named
    # exactly as before -- no error, matching this project's long-standing
    # convention for unspecified stereochemistry.
    assert smiles_to_iupac("CCC(C)S(=O)O") == "butane-2-sulfinic acid"


def test_sulfinic_acid_specified_chain_stereocenter_raises():
    # The sulfinic sulfur (-R, =O, -OH) is itself a potential stereocenter
    # in this molecule too (unlike `_sulfonic_acid.py`'s sulfur, whose two
    # identical =O make it never stereogenic) -- so a specified chain
    # stereocenter here always coexists with an unspecified sulfur one,
    # and `specified_stereocenters` correctly rejects the combination
    # (P-92) instead of the silent drop this project's stereodescriptor
    # safety net exists to fix.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)S(=O)O")


def test_sulfinic_acid_specified_sulfur_stereocenter_cited():
    # A specified sulfur configuration, with no chain stereocenter
    # alongside it, gets a bare (R)-/(S)- prefix per P-93.3.4.1 (the
    # sulfinic sulfur is a trigonal pyramidal stereocenter, cited the same
    # way as a tetrahedral one). PubChem itself doesn't distinguish the
    # two sulfur configurations for this molecule (CCC[S@](=O)O and
    # CCC[S@@](=O)O both resolve to CID 643586 with the same, unstereo
    # name) -- confirmed instead against RDKit's own CIP assignment.
    assert smiles_to_iupac("CCC[S@](=O)O") == "(R)-propane-1-sulfinic acid"
    assert smiles_to_iupac("CCC[S@@](=O)O") == "(S)-propane-1-sulfinic acid"


def test_ring_sulfinic_acid_non_stereocenter_marker_unaffected():
    # The ring carbon bearing -SO2H here isn't a genuine stereocenter (both
    # ring neighbors are identical -CH2- groups), so RDKit never reports it
    # as a stereo element at all -- this stray '@' marker is silently
    # ignored exactly as before, same policy as every other module's
    # "unspecified/non-stereogenic marker" handling.
    assert smiles_to_iupac("O=S(O)[C@H]1CCCCC1") == "cyclohexanesulfinic acid"


def test_phenyl_chain_sulfinic_acid():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#278's
    # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/acyl-
    # halide/sulfonic-acid/thiol chains): the ring is cited as a "phenyl"
    # substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCCS(=O)O") == "3-phenylpropane-1-sulfinic acid"


def test_benzenesulfinic_acid():
    # -SO2H directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1S(=O)O") == "benzenesulfinic acid"  # CID 12057


def test_substituted_benzenesulfinic_acid():
    # A substituent on a different ring atom than the -SO2H: the
    # mancude-ring numbering is free to start at the -SO2H carbon, so its
    # own locant is never cited, unlike the cycloalkane case.
    assert smiles_to_iupac("Cc1ccccc1S(=O)O") == "2-methylbenzenesulfinic acid"  # CID 12661295


def test_phenyl_substituted_benzene_ring_sulfinic_acid_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-2.md).
    assert smiles_to_iupac("Cc1ccccc1CCS(=O)O") == "2-(2-methylphenyl)ethanesulfinic acid"


def test_phenyl_chain_sulfinic_acid_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propane-1-sulfinic
    # acid" for this structure.
    assert smiles_to_iupac("Clc1ccc(CCCS(=O)O)cc1") == "3-(4-chlorophenyl)propane-1-sulfinic acid"


def test_phenyl_chain_sulfinic_acid_ring_methyl():
    # PubChem PUG REST-verified "3-(4-methylphenyl)propane-1-sulfinic acid".
    assert smiles_to_iupac("Cc1ccc(CCCS(=O)O)cc1") == "3-(4-methylphenyl)propane-1-sulfinic acid"


def test_phenyl_chain_sulfinic_acid_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CS(=O)O") == "(4-ethylphenyl)methanesulfinic acid"


def test_phenyl_chain_sulfinic_acid_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCS(=O)O")
