import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Same skeletons/locants as `test_spiro_amine.py`'s
        # PubChem-confirmed cases, with the ring carbon itself carrying
        # the +1 charge instead of an exocyclic -NH2 substituent -- see
        # `test_von_baeyer_carbenium.py` for why verification here
        # cross-checks the neutral parent-hydride skeleton name instead
        # of the charged species directly.
        ("[CH+]1CCC2(CC1)CCCC2", "spiro[4.5]decan-8-ylium"),
        ("[CH+]1CCC2(CC1)CCCCC2", "spiro[5.5]undecan-3-ylium"),
        ("[CH+]1CCC2(CCC2)CC1", "spiro[3.5]nonan-7-ylium"),
    ],
)
def test_spiro_carbenium_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_spiro_carbenium_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CCC2(CC1)CCCC2[CH2+]")


def test_spiro_carbenium_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CCC2(CC1)CCCC2")


def test_unsaturated_spiro_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH+]1CCC2(C=CC2)CC1")


def test_stereo_spiro_carbenium_resolves():
    # CIP label cross-checked independently via rdCIPLabeler.
    assert smiles_to_iupac("[CH+]1CCC2(CC1)CCC[C@@H](C)C2") == "(8R)-8-methylspiro[5.5]undecan-3-ylium"
