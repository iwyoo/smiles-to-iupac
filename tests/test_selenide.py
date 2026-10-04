import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_selenide():
    # PubChem PUG REST CID 11648, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Se]C") == "methylselanylmethane"


def test_diethyl_selenide():
    # PubChem PUG REST CID 61173, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Se]CC") == "ethylselanylethane"


def test_ethyl_methyl_selenide():
    # The shorter (methyl) side becomes the substituent, the longer
    # (ethyl) side the parent. PubChem PUG REST CID 12248622,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("C[Se]CC") == "methylselanylethane"


def test_methyl_propyl_selenide():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 15932888,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[Se]C") == "1-methylselanylpropane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'selanyl' outside. Structure verified against PubChem: CID 11321295
    # ("1-propan-2-ylselanylbutane" -- PubChem elides the parentheses this
    # project's usual compound-substituent formatting keeps).
    assert smiles_to_iupac("CCCC[Se]C(C)C") == "1-(propan-2-yl)selanylbutane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Same tied-parent-selection cases as `test_ether.py`, selanyl in
        # place of oxy -- see that module's tests for the full reasoning.
        # PubChem can't compute a name for the branched/branched cases here
        # (CID 0, no IUPACName returned for either structure), so these are
        # verified structurally instead: `_selenide.py` reuses the exact
        # same graph algorithm already PubChem-confirmed for `_ether.py`/
        # `_sulfide.py`, with only the prefix word ('selanyl') differing.
        ("CC(C)[Se]C(C)C", "2-(propan-2-yl)selanylpropane"),
        ("CC(C)C[Se]C(C)(C)C", "1-tert-butylselanyl-2-methylpropane"),
        ("CC(C)(C)C[Se]CCC(C)C", "1-(2,2-dimethylpropyl)selanyl-3-methylbutane"),
        ("CC(C)C[Se]C(C)CC", "2-(2-methylpropyl)selanylbutane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-44.1.2.2 rule (1): 'selanyl' has no suffix form, so a plain
        # benzene ring is always the parent, mirroring `_sulfide.py`'s
        # identical benzene-ring path.
        ("c1ccccc1[Se]CC", "ethylselanylbenzene"),  # PubChem CID 140285
        # Direct ring-selenium bond with a branched R': the same
        # '(...)selanyl' parenthesization as the plain two-chain path
        # above, even though PubChem's own auto-generated name for the
        # same structure omits the parentheses
        # ('c1ccccc1[Se]C(C)C' -> 'propan-2-ylselanylbenzene', CID 140894).
        ("c1ccccc1[Se]C(C)C", "(propan-2-yl)selanylbenzene"),
        # Chain spacer between the ring and the selenide selenium --
        # matching this project's existing benzene-ring-chain convention
        # rather than PubChem's own un-parenthesized auto-name
        # ('ethylselanylmethylbenzene', CID 12975596).
        ("c1ccccc1C[Se]CC", "(ethylselanylmethyl)benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents_named_with_prefix():
    assert smiles_to_iupac("C[Se]c1ccccc1[Se]C") == "1,2-bis(methylselanyl)benzene"


def test_benzene_ring_stereocenter_named_with_prefix():
    assert smiles_to_iupac("c1ccccc1[Se][C@H](C)CC") == "{[(2R)-butan-2-yl]selanyl}benzene"
