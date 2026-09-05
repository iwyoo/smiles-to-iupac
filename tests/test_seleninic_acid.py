import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-seleninic acid' mirrors '-selenonic acid' (see
        # test_selenonic_acid.py) with one fewer oxygen; same locant rules
        # as '-sulfinic acid'.
        ("C[Se](=O)O", "methaneseleninic acid"),
        ("CC[Se](=O)O", "ethaneseleninic acid"),
        ("CCC[Se](=O)O", "propane-1-seleninic acid"),
        ("CC([Se](=O)O)C", "propane-2-seleninic acid"),
        ("CCCC[Se](=O)O", "butane-1-seleninic acid"),
    ],
)
def test_saturated_seleninic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_seleninic_acid():
    assert smiles_to_iupac("C=CC[Se](=O)O") == "prop-2-ene-1-seleninic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)[Se](=O)O") == "1-chloroethane-1-seleninic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C([Se](=O)O)C")


def test_two_seleninic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)C[Se](=O)O")


def test_sulfinic_acid_not_confused_with_seleninic_acid():
    assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"


def test_selenonic_acid_not_confused_with_seleninic_acid():
    assert smiles_to_iupac("C[Se](=O)(=O)O") == "methaneselenonic acid"


def test_ring_seleninic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)C1CCCCC1")


def test_seleninic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Se](=O)CCO")


def test_seleninic_acid_unspecified_stereocenter_unaffected():
    # A genuine chain stereocenter left unspecified (no @/@@) is named
    # exactly as before -- no error, matching this project's long-standing
    # convention for unspecified stereochemistry.
    assert smiles_to_iupac("CCC(C)[Se](=O)O") == "butane-2-seleninic acid"


def test_seleninic_acid_specified_chain_stereocenter_raises():
    # The seleninic selenium (-R, =O, -OH) is itself a potential
    # stereocenter in this molecule too (unlike `_sulfonic_acid.py`'s
    # sulfur or `_tellurinic_acid.py`'s tellurium), so a specified chain
    # stereocenter here always coexists with an unspecified selenium one,
    # and `specified_stereocenters` correctly rejects the combination
    # (P-92) instead of the silent drop this project's stereodescriptor
    # safety net exists to fix (same reasoning as `_sulfinic_acid.py`).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)[Se](=O)O")


def test_phenyl_chain_seleninic_acid():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring
    # `_sulfinic_acid.py`'s identical PR): the ring is cited as a
    # "phenyl" substituent prefix. PubChem CID 125616
    # ("phenylmethaneseleninic acid") confirms the mononuclear (no
    # locant) case.
    assert smiles_to_iupac("O=[Se](O)Cc1ccccc1") == "phenylmethaneseleninic acid"
    assert smiles_to_iupac("c1ccccc1CCC[Se](=O)O") == "3-phenylpropane-1-seleninic acid"


def test_benzeneseleninic_acid():
    # -Se(=O)OH directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1[Se](=O)O") == "benzeneseleninic acid"


def test_substituted_benzeneseleninic_acid():
    # The mancude-ring numbering is free to start at the -Se(=O)OH
    # carbon, so its own locant is never cited, mirroring
    # benzenesulfonic acid.
    assert smiles_to_iupac("Cc1ccccc1[Se](=O)O") == "2-methylbenzeneseleninic acid"
    assert smiles_to_iupac("Cc1ccc(cc1)[Se](=O)O") == "4-methylbenzeneseleninic acid"  # PubChem PUG REST


def test_benzeneseleninic_acid_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[Se@](=O)O")


def test_phenyl_substituted_benzene_ring_seleninic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[Se](=O)O")


def test_phenyl_chain_seleninic_acid_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[Se](=O)O")
