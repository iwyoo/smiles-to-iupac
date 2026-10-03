import warnings

import pytest

from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac


@pytest.mark.parametrize(
    "smiles",
    ["C[Ti](Cl)(Cl)Cl", "C[Li]", "C[Zn]C", "CC1C[Pt](Cl)(Cl)C1"],
)
def test_organometallic_without_a_pin_warns(smiles):
    with pytest.warns(NonPreferredNameWarning, match="no PIN"):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize("smiles", ["CC[Al](CC)CC", "CCO", "[Fe+2].C[c-]1cccc1.[cH-]1cccc1"])
def test_names_with_a_pin_do_not_warn(smiles):
    with warnings.catch_warnings():
        warnings.simplefilter("error", NonPreferredNameWarning)
        smiles_to_iupac(smiles)


def test_warning_is_emitted_once_per_call():
    with pytest.warns(NonPreferredNameWarning) as record:
        smiles_to_iupac("COC(=O)C(C)CC1C[Pt](C)(I)(P(CC)(CC)CC)(P(CC)(CC)CC)C1")
    assert len(record) == 1
