import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_dimethyl_telluride():
    # PubChem PUG REST CID 68977, auto-generated name matches exactly.
    assert smiles_to_iupac("C[Te]C") == "methyltellanylmethane"


def test_diethyl_telluride():
    # PubChem PUG REST CID 69394, auto-generated name matches exactly.
    assert smiles_to_iupac("CC[Te]CC") == "ethyltellanylethane"


def test_methyl_ethyl_telluride():
    # The shorter (methyl) side becomes the substituent, the longer
    # (ethyl) side the parent. PubChem PUG REST CID 13981584,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("C[Te]CC") == "methyltellanylethane"


def test_methyl_propyl_telluride():
    # A 3-carbon parent needs the locant. PubChem PUG REST CID 15932889,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[Te]C") == "1-methyltellanylpropane"


def test_branched_prefix_side_is_enclosed():
    # P-63.2.2.1.1: a branched R' encloses only R' in parentheses, with
    # 'tellanyl' outside. Structure verified against PubChem: CID 13975014
    # ("1-propan-2-yltellanylbutane" -- PubChem elides the parentheses this
    # project's usual compound-substituent formatting keeps).
    assert smiles_to_iupac("CCCC[Te]C(C)C") == "1-(propan-2-yl)tellanylbutane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Same tied-parent-selection cases as `test_ether.py`, tellanyl in
        # place of oxy -- see that module's tests for the full reasoning.
        # PubChem can't compute a name for the branched/branched cases here
        # (CID 0, no IUPACName returned for either structure), so these are
        # verified structurally instead: `_telluride.py` reuses the exact
        # same graph algorithm already PubChem-confirmed for `_ether.py`/
        # `_sulfide.py`, with only the prefix word ('tellanyl') differing.
        ("CC(C)[Te]C(C)C", "2-(propan-2-yl)tellanylpropane"),
        ("CC(C)C[Te]C(C)(C)C", "1-tert-butyltellanyl-2-methylpropane"),
        ("CC(C)(C)C[Te]CCC(C)C", "1-(2,2-dimethylpropyl)tellanyl-3-methylbutane"),
        ("CC(C)C[Te]C(C)CC", "2-(2-methylpropyl)tellanylbutane"),
    ],
)
def test_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-44.1.2.2 rule (1): 'tellanyl' has no suffix form, so a plain
        # benzene ring is always the parent, mirroring `_selenide.py`'s/
        # `_sulfide.py`'s identical benzene-ring path.
        ("c1ccccc1[Te]CC", "ethyltellanylbenzene"),  # PubChem CID 5325650
        # Direct ring-tellurium bond with a branched R': the same
        # '(...)tellanyl' parenthesization as the plain two-chain path
        # above, even though PubChem's own auto-generated name for the
        # same structure omits the parentheses.
        ("c1ccccc1[Te]C(C)C", "(propan-2-yl)tellanylbenzene"),
        # Chain spacer between the ring and the telluride tellurium --
        # matching this project's existing benzene-ring-chain convention
        # rather than PubChem's own un-parenthesized auto-name
        # ('ethyltellanylmethylbenzene', CID 23077050).
        ("c1ccccc1C[Te]CC", "(ethyltellanylmethyl)benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents_named_with_prefix():
    assert smiles_to_iupac("C[Te]c1ccccc1[Te]C") == "1,2-bis(methyltellanyl)benzene"


def test_benzene_ring_stereocenter_named_with_prefix():
    assert smiles_to_iupac("c1ccccc1[Te][C@H](C)CC") == "{[(2R)-butan-2-yl]tellanyl}benzene"
