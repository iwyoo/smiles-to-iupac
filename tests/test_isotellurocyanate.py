import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC[N]=C=[Te]", "isotellurocyanatoethane"),
        ("CCC[N]=C=[Te]", "1-isotellurocyanatopropane"),
        ("C[N]=C=[Te]", "isotellurocyanatomethane"),
    ],
)
def test_isotellurocyanate(smiles, expected):
    # Unlike `_isocyanate.py`/`_isothiocyanate.py`/`_isoselenocyanate.py`,
    # no PubChem-listed compound was found for any of these three
    # structures (all CID 0), so these are reviewed (eyeballed), not
    # independently verified, results -- a mechanical extension of the
    # identical pattern already confirmed twice over for each of the
    # oxygen/sulfur/selenium chalcogen analogues (see module docstring).
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC([N]=C=[Te])CC1")


def test_phenyl_isotellurocyanate_direct_bond():
    # P-44.1.2.2 rule (1): 'isotellurocyanato' has no suffix form, so the
    # ring is always senior to a chain of the same class -- not
    # independently PubChem-verified (CID 0, same sparse-tellurium gap as
    # the acyclic path), inherited unchanged from the identical mechanism
    # already confirmed for oxygen/sulfur/selenium.
    assert smiles_to_iupac("c1ccccc1[N]=C=[Te]") == "isotellurocyanatobenzene"


def test_phenyl_isotellurocyanate_chain():
    assert smiles_to_iupac("c1ccccc1C[N]=C=[Te]") == "(isotellurocyanatomethyl)benzene"


def test_phenyl_isotellurocyanate_substituted_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1[N]=C=[Te]")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC[N]=C=[Te]")


def test_two_isotellurocyanate_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]=C=NCN=C=[Te]")


def test_isoselenocyanate_still_works():
    # Sanity check: the selenium analogue must not be misrouted here.
    assert smiles_to_iupac("CC[N]=C=[Se]") == "isoselenocyanatoethane"
