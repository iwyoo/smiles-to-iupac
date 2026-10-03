import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Same skeletons/locants as `test_von_baeyer_amine.py`'s
        # PubChem-confirmed cases, with the ring carbon itself carrying
        # the +1 charge (losing one hydrogen) instead of an exocyclic
        # -NH2 substituent -- these charged structures aren't reliably
        # registered/named by PubChem (see `_carbenium.py`'s module
        # docstring), so verification here cross-checks the neutral
        # parent-hydride skeleton name (already PubChem-confirmed via the
        # sibling amine/alcohol modules on the identical skeletons) plus
        # the shared suffix-locant mechanism reused unchanged from those
        # modules.
        ("[CH+]1CC2CCC1CC2", "bicyclo[2.2.2]octan-2-ylium"),
        ("[CH+]1CCC2CCC1C2", "bicyclo[3.2.1]octan-2-ylium"),
        ("[CH+]1CC2CCC1C2", "bicyclo[2.2.1]heptan-2-ylium"),
        ("[CH+]1CC2CCC(C1)C2", "bicyclo[3.2.1]octan-3-ylium"),
        ("[CH+]1CC2CCCC(C1)C2", "bicyclo[3.3.1]nonan-3-ylium"),
        ("[CH+]1CC2CCC2C1", "bicyclo[3.2.0]heptan-3-ylium"),
        # A bridgehead cation (degree 3, no hydrogen) rather than a
        # secondary ring position -- same tricyclic skeleton and
        # systematic-name convention as `test_von_baeyer_amine.py`'s
        # 1-adamantanamine analog.
        ("[C+]12CC3CC(CC(C3)C1)C2", "adamantan-1-ylium"),
    ],
)
def test_von_baeyer_carbenium_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_carbenium_centers_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CC2CCC1C2[CH2+]")


def test_carbenium_on_substituent_branch_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2+]C1CC2CCC1C2")


def test_unsaturated_von_baeyer_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH+]1CC2C=CC1C2")


def test_stereo_von_baeyer_carbenium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH+]1C[C@@H]2CCC1C2")
