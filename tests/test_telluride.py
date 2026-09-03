import pytest

from smiles_to_iupac import smiles_to_iupac


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
