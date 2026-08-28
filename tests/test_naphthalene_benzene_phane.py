import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_naphthalene_benzene_unequal_bridge_phane():
    # Blue Book P-26.4.1.4's own worked example:
    # "1(2,7)-naphthalena-4(1,4)-benzenacycloheptaphane (PIN)" -- a
    # naphthalene superatom attached at its 2,7 positions and a benzene
    # superatom attached at its para (1,4) positions, joined by bridges of
    # length 2 and 3. No PubChem-listed compound exists for this exotic
    # macrocycle (checked), so this is settled by the primary-source
    # worked example itself.
    assert (
        smiles_to_iupac("c1cc2ccc1CCCc1ccc3ccc(cc3c1)CC2")
        == "1(2,7)-naphthalena-4(1,4)-benzenacycloheptaphane"
    )


def test_substituted_naphthalene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cc2ccc3cc2cc1CCc1ccc(cc1)CCC3")


def test_heteroaromatic_second_ring_raises():
    # pyridine, not benzene, as the second superatom -- out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1cc2ncc1CCCc1ccc3ccc(cc3c1)CC2")
