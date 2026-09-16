import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-tellurinic acid' mirrors '-tellurionic acid'/'-seleninic acid'
        # (see test_telluronic_acid.py/test_seleninic_acid.py) with one
        # fewer oxygen; same locant rules as '-sulfinic acid'. Only the
        # methane case is registered on PubChem (ethane/propane come back
        # as CID 0) -- the same sparse-Te-data gap test_telluronic_acid.py
        # already documented.
        ("C[Te](=O)O", "methanetellurinic acid"),
        ("CC[Te](=O)O", "ethanetellurinic acid"),
        ("CCC[Te](=O)O", "propane-1-tellurinic acid"),
        ("CC([Te](=O)O)C", "propane-2-tellurinic acid"),
        ("CCCC[Te](=O)O", "butane-1-tellurinic acid"),
    ],
)
def test_saturated_tellurinic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_tellurinic_acid():
    assert smiles_to_iupac("C=CC[Te](=O)O") == "prop-2-ene-1-tellurinic acid"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(Cl)[Te](=O)O") == "1-chloroethanetellurinic acid"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C([Te](=O)O)C")


def test_two_tellurinic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)C[Te](=O)O")


def test_telluronic_acid_not_confused_with_tellurinic_acid():
    assert smiles_to_iupac("C[Te](=O)(=O)O") == "methanetelluronic acid"


def test_ring_tellurinic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)C1CCCCC1")


def test_tellurinic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O[Te](=O)CCO")


def test_benzenetellurinic_acid():
    # -Te(=O)OH directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST IUPACName.
    assert smiles_to_iupac("c1ccccc1[Te](=O)O") == "benzenetellurinic acid"


def test_substituted_benzenetellurinic_acid():
    # The mancude-ring numbering is free to start at the -Te(=O)OH
    # carbon, so its own locant is never cited, mirroring
    # benzeneseleninic/benzenesulfinic acid.
    assert smiles_to_iupac("Cc1ccccc1[Te](=O)O") == "2-methylbenzenetellurinic acid"
    assert smiles_to_iupac("Cc1ccc(cc1)[Te](=O)O") == "4-methylbenzenetellurinic acid"  # PubChem PUG REST


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): unlike
        # `_seleninic_acid.py`'s selenium, RDKit's `Chem.FindPotentialStereo`
        # never flags this module's tellurinic tellurium as a potential
        # stereocenter (module docstring) -- so this mirrors
        # `_sulfonic_acid.py`/`_carboxylic_acid.py` cleanly (CIP computed
        # entirely by RDKit's `rdCIPLabeler`). PubChem has no registered
        # tellurinic acid stereoisomer (or even most non-stereo ones --
        # module docstring), so this is a structural/regression check on
        # the already-proven mechanism, not an independent PubChem
        # cross-check.
        ("CC[C@@H](C)[Te](=O)O", "(2R)-butane-2-tellurinic acid"),
        ("CC[C@H](C)[Te](=O)O", "(2S)-butane-2-tellurinic acid"),
    ],
)
def test_tellurinic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_tellurinic_acid_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)[Te](=O)O") == "butane-2-tellurinic acid"
