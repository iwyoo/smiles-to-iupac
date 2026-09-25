import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem D-aldohexofuranose structures, both anomers,
        # spanning distinct stem names (#1039 M2 step 2).
        ("C([C@H]([C@@H]1[C@@H]([C@H]([C@H](O1)O)O)O)O)O", "alpha-D-glucofuranose"),  # CID 11137711
        ("C([C@H]([C@@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O)O", "beta-D-glucofuranose"),  # CID 11309871
        (
            "C([C@H]([C@@H]1[C@@H]([C@@H]([C@H](O1)O)O)O)O)O",
            "alpha-D-mannofuranose",
        ),  # CID 9920533
        (
            "C([C@H]([C@H]1[C@@H]([C@H]([C@H](O1)O)O)O)O)O",
            "alpha-D-galactofuranose",
        ),  # CID 21627868
    ],
)
def test_cyclic_aldofuranose_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_pyranose_still_resolves_unchanged():
    # Regression guard: the 6-membered ring form must still dispatch to
    # `name_cyclic_aldopyranose`, not accidentally to the new furanose path.
    assert smiles_to_iupac("C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O") == "alpha-D-glucopyranose"
