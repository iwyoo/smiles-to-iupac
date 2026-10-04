import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem D-2-ketohexopyranose structures, spanning all 4
        # retained D-2-ketohexose stem names and both anomers where a
        # registered structure exists (#1039 M2 step 3).
        ("C1[C@H]([C@H]([C@@H]([C@@](O1)(CO)O)O)O)O", "alpha-D-fructopyranose"),  # CID 439709
        ("C1[C@H]([C@H]([C@@H]([C@](O1)(CO)O)O)O)O", "beta-D-fructopyranose"),  # CID 24310
        ("C1[C@H]([C@@H]([C@@H]([C@@](O1)(CO)O)O)O)O", "alpha-D-tagatopyranose"),
        ("C1[C@H]([C@@H]([C@@H]([C@](O1)(CO)O)O)O)O", "beta-D-tagatopyranose"),
        ("C1[C@H]([C@@H]([C@H]([C@](O1)(CO)O)O)O)O", "beta-D-sorbopyranose"),
        ("C1[C@H]([C@H]([C@H]([C@](O1)(CO)O)O)O)O", "beta-D-psicopyranose"),
    ],
)
def test_cyclic_ketohexopyranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aldopyranose_still_resolves():
    # A plain aldohexopyranose (anomeric carbon bears only -OH, not also a
    # -CH2OH branch) must still route to the aldose path unchanged.
    assert smiles_to_iupac("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O") == "alpha-D-glucopyranose"


def test_hetero_ring_ketone_still_resolves():
    # A plain saturated single-heteroatom ring bearing an actual ring
    # ketone (not a hemiketal) must still route to `_ketone.py` unchanged.
    assert smiles_to_iupac("O=C1CCCCO1") == "oxan-2-one"


def test_ketofuranose_cites_the_specified_elements():
    # A 5-membered cyclic 2-ketose hemiketal (furanose) is a different
    # ring size, out of scope here (later M2 step) -- still collides with
    # `_ketone.py`'s hetero-ring-ketone path unchanged.
    assert smiles_to_iupac("C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O") == '(2R,3S,4S,5R)-2,5-bis(hydroxymethyl)oxolane-2,3,4-triol'