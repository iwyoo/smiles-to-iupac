import pytest
from rdkit import Chem

from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._fusion_numbering import derive_letter_by_pair, letter_by_pair_candidates
from smiles_to_iupac._tetracene_fusion import _LETTER_BY_PAIR as _TETRACENE_LETTER_BY_PAIR
from smiles_to_iupac._tetracene_fusion import _TETRACENE_REF

_PHENANTHRENE_REF = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
# PubChem CID 41601 (1-methylphenanthrene), CID 68154 (9-methylphenanthrene),
# CID 68162 (4-methylphenanthrene).
_PHEN_1_METHYL = Chem.MolFromSmiles("CC1=C2C=CC3=CC=CC=C3C2=CC=C1")
_PHEN_9_METHYL = Chem.MolFromSmiles("CC1=CC2=CC=CC=C2C3=CC=CC=C13")
_PHEN_4_METHYL = Chem.MolFromSmiles("CC1=C2C(=CC=C1)C=CC3=CC=CC=C32")

# PubChem CID 23505253/15660752/23505242 (1-/2-/5-methyltetracene), the same
# anchors `_tetracene_fusion.py`'s own docstring describes using by hand.
_TET_1_METHYL = Chem.MolFromSmiles("CC1=CC=CC2=CC3=CC4=CC=CC=C4C=C3C=C12")
_TET_2_METHYL = Chem.MolFromSmiles("CC1=CC2=CC3=CC4=CC=CC=C4C=C3C=C2C=C1")
_TET_5_METHYL = Chem.MolFromSmiles("CC1=C2C=CC=CC2=CC3=CC4=CC=CC=C4C=C13")


def test_single_anchor_leaves_many_candidates():
    with pytest.raises(UnsupportedStructure, match="numbering candidates undecided"):
        derive_letter_by_pair(_PHENANTHRENE_REF, [(_PHEN_1_METHYL, 1)])


def test_symmetric_base_never_converges_past_its_automorphism_group():
    # Phenanthrene's only nontrivial automorphism is its own mirror
    # symmetry (order 2) -- no combination of real anchors can ever
    # narrow this below 2 candidates (see module docstring), and this
    # must raise a clear error rather than silently pick one.
    with pytest.raises(UnsupportedStructure, match="numbering candidates undecided"):
        derive_letter_by_pair(
            _PHENANTHRENE_REF, [(_PHEN_1_METHYL, 1), (_PHEN_9_METHYL, 9), (_PHEN_4_METHYL, 4)]
        )
    candidates = letter_by_pair_candidates(
        _PHENANTHRENE_REF, [(_PHEN_1_METHYL, 1), (_PHEN_9_METHYL, 9), (_PHEN_4_METHYL, 4)]
    )
    assert len(candidates) == 2
    # The historically-hardcoded numbering must still be one of the two
    # mirror-symmetric candidates the algorithm correctly narrows down to.
    from smiles_to_iupac._phenanthrene_fusion import _LETTER_BY_PAIR as _PHEN_LETTER_BY_PAIR

    assert _PHEN_LETTER_BY_PAIR in candidates


def test_straight_acene_family_hardcoded_numbering_is_among_the_candidates():
    # Tetracene has the same D2h (order-4) symmetry every plain straight
    # acene does -- even with all three of the real anchors
    # `_tetracene_fusion.py`'s own docstring used by hand, 4 candidates
    # remain (one per symmetry-group element), and the historically-
    # hardcoded one must be among them.
    candidates = letter_by_pair_candidates(
        _TETRACENE_REF, [(_TET_1_METHYL, 1), (_TET_2_METHYL, 2), (_TET_5_METHYL, 5)]
    )
    assert len(candidates) == 4
    assert _TETRACENE_LETTER_BY_PAIR in candidates


def test_inconsistent_anchors_raise():
    # A locant claim that cannot be satisfied by any valid periphery
    # numbering at all (not even a symmetric alternative).
    with pytest.raises(UnsupportedStructure, match="not consistent with any valid numbering"):
        derive_letter_by_pair(_PHENANTHRENE_REF, [(_PHEN_1_METHYL, 1), (_PHEN_1_METHYL, 2)])


def test_anchor_not_matching_reference_raises():
    naphthalene_1_methyl = Chem.MolFromSmiles("Cc1cccc2ccccc12")
    with pytest.raises(UnsupportedStructure, match="does not match the reference molecule"):
        derive_letter_by_pair(_PHENANTHRENE_REF, [(naphthalene_1_methyl, 1)])
