import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methaneselenonic_acid():
    # PubChem structure match: "methaneselenonic acid".
    assert smiles_to_iupac("C[Se](=O)(=O)O") == "methaneselenonic acid"


def test_ethaneselenonic_acid():
    # PubChem structure match: "ethaneselenonic acid".
    assert smiles_to_iupac("CC[Se](=O)(=O)O") == "ethaneselenonic acid"


def test_propane_1_selenonic_acid():
    # PubChem structure match: "propane-1-selenonic acid".
    assert smiles_to_iupac("CCC[Se](=O)(=O)O") == "propane-1-selenonic acid"


def test_propane_2_selenonic_acid():
    assert smiles_to_iupac("CC([Se](=O)(=O)O)C") == "propane-2-selenonic acid"


def test_chlorobutaneselenonic_acid():
    assert smiles_to_iupac("ClCCCC[Se](=O)(=O)O") == "4-chlorobutane-1-selenonic acid"


def test_pent_4_ene_1_selenonic_acid():
    assert smiles_to_iupac("C=CCCC[Se](=O)(=O)O") == "pent-4-ene-1-selenonic acid"


def test_diselenonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)(=O)C[Se](=O)(=O)O")


def test_sulfonic_acid_not_confused_with_selenonic_acid():
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_ring_selenonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)(=O)C1CCCCC1")


def test_selenonic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)(=O)CCO")


def test_benzeneselenonic_acid():
    # -Se(=O)(=O)OH directly on a benzene ring carbon, cross-checked
    # against PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1[Se](=O)(=O)O") == "benzeneselenonic acid"


def test_substituted_benzeneselenonic_acid():
    # The mancude-ring numbering is free to start at the -Se(=O)(=O)OH
    # carbon, so its own locant is never cited, mirroring
    # benzenesulfonic acid.
    assert smiles_to_iupac("Cc1ccccc1[Se](=O)(=O)O") == "2-methylbenzeneselenonic acid"
    assert smiles_to_iupac("Cc1ccc(cc1)[Se](=O)(=O)O") == "4-methylbenzeneselenonic acid"  # PubChem PUG REST


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_sulfonic_acid.py`'s acyclic stereocenter task -- selenium
        # itself is not a potential stereocenter (its two double-bonded
        # oxygens are identical), confirmed via RDKit `FindPotentialStereo`.
        # PubChem has no registered selenonic acid with a branch (even
        # unspecified), so only structure/CIP-label consistency with
        # `_sulfonic_acid.py`'s sulfur analogue is checked here.
        ("CC[C@@H](C)[Se](=O)(=O)O", "(2R)-butane-2-selenonic acid"),
        ("CC[C@H](C)[Se](=O)(=O)O", "(2S)-butane-2-selenonic acid"),
    ],
)
def test_acyclic_selenonic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_selenonic_acid_stereocenter_with_coexisting_substituent():
    # A stereocenter that also bears a halogen substituent: the suffix's
    # own locant is still forced to be C1 by P-44.4.1.8, so it's omitted
    # here too (P-14.3.4.2(b)) even though the halogen shares that
    # position, mirroring `_sulfonic_acid.py`'s identical rule.
    assert smiles_to_iupac("C[C@@H](Cl)[Se](=O)(=O)O") == "(1S)-1-chloroethaneselenonic acid"


def test_selenonic_acid_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(Cl)[Se](=O)(=O)O") == "1-chloropropane-1-selenonic acid"
