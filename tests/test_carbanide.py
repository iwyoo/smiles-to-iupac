import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-72.2.2.1's own worked examples: 'methanide (PIN)'
        # for H3C- (line 889), and P-71.1's own text confirms 'ethanide'
        # (not 'ethyl anion') for CH3-CH2- (line 131) -- the suffix
        # locant is omitted at one and two carbons, mirroring
        # `_alcohol.py`'s '-ol' citation.
        ("[CH3-]", "methanide"),
        ("C[CH2-]", "ethanide"),
        # A 3+-carbon chain always cites the suffix locant (P-72.2.2.1's
        # own worked '2-(...)propan-2-ide' examples, e.g. line 2872).
        ("CC[CH2-]", "propan-1-ide"),
        ("C[CH-]C", "propan-2-ide"),
        ("CCCC[CH2-]", "pentan-1-ide"),
    ],
)
def test_carbanide_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_carbanide():
    # A branched skeleton is named via the same longest-chain-through-
    # the-anion-carbon search `_alkoxide.py` uses for its own oxygen
    # substituent -- e.g. isobutyl anion, "2-methylpropan-1-ide".
    assert smiles_to_iupac("CC([CH2-])C") == "2-methylpropan-1-ide"


def test_isocyanide_carbon_not_misnamed_as_carbanide():
    # An isocyanide carbon (-N+#C-, P-61.9) also carries a formal charge
    # -1, but it's not a plain hydrocarbon carbanion -- `_isocyanide.py`
    # already handles this shape, so it must not be misrouted here.
    assert smiles_to_iupac("C[N+]#[C-]") == "isocyanomethane"


def test_multiple_carbanion_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2-][CH2-]")


def test_carbanide_with_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH-]1CCCCC1")


def test_carbanide_with_halogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC[CH2-]")


def test_carbanide_with_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[CH2-]")


def test_carbanide_with_heteroatom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC[CH2-]")


def test_acetyl_anion_with_oxo_on_anion_carbon():
    # P-72.2.2.1's own worked example: 'acetyl anion' -> '1-oxoethan-1-ide
    # (PIN)' -- unlike the plain two-carbon 'ethanide' above, the '-ide'
    # suffix locant is cited here since it coincides with the cited 'oxo'
    # substituent locant (tmp/bluebook/P7.txt lines 887-891/1218).
    assert smiles_to_iupac("C[C-]=O") == "1-oxoethan-1-ide"


def test_oxo_substituent_elsewhere_on_chain():
    assert smiles_to_iupac("CCC(=O)[CH-]C") == "3-oxopentan-2-ide"


def test_carbanide_with_second_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC(=O)[CH-]C")


def test_carbanide_with_aldehyde_shaped_carbonyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C[CH-]C")
