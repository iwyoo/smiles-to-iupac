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


def test_pyridine_benzene_equal_bridge_phane():
    # A second registered senior amplificant kind (#1031): pyridine
    # attached at its own 2,5 positions, joined to a para-attached benzene
    # by two equal-length bridges. No literal Blue Book worked example or
    # PubChem structure exists for this exact pairing -- derived by hand
    # from P-26.4.1.1/4.1.2/4.1.4 and P-44.2 (pyridine outranks benzene,
    # same as naphthalene). The '(2,5)' citation order (not '(5,2)')
    # exercises P-26.4.1.4's tie-break explicitly: pyridine has less
    # symmetry than naphthalene (only one mirror automorphism), so unlike
    # every naphthalene case tested here, the two bridge directions don't
    # trivially agree on which attachment locant is lower.
    assert (
        smiles_to_iupac("c1cc2ccc1CCc1ccc(nc1)CC2")
        == "1(2,5)-pyridina-4(1,4)-benzenacyclohexaphane"
    )


def test_pyridine_benzene_unequal_bridge_phane():
    # Same pyridine/benzene attachment pattern, bridges 2 and 3 instead of
    # 2 and 2 -- the shorter bridge goes forward (P-26.4.1.1), same rule
    # already verified for naphthalene.
    assert (
        smiles_to_iupac("c1cc2ccc1CCCc1ccc(nc1)CC2")
        == "1(2,5)-pyridina-4(1,4)-benzenacycloheptaphane"
    )


def test_pyridinium_nitrogen_attachment_raises():
    # Pyridine's own ring nitrogen (locant 1) can't be a bridge-attachment
    # point while neutral -- attaching a bridge there requires a charged
    # pyridinium-type nitrogen, already out of scope (charged atoms are
    # rejected everywhere in this project).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1Cc2ccc(cc2)CC[n+]2ccccc21")
