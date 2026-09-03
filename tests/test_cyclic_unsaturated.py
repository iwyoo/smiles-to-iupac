import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Plain unsubstituted cycloalkenes/dienes are well-known, unambiguous
        # compound names (cyclohexene, cyclohexadiene, cyclopentene,
        # cycloheptene); the substituted variants below are covered by
        # in-line locant-rule citations instead of a per-case lookup.
        ("C1=CCCCC1", "cyclohexene"),
        ("C1=CC=CCC1", "cyclohexa-1,3-diene"),
        ("C1=CCC=CC1", "cyclohexa-1,4-diene"),
        ("C1=CCCC1", "cyclopentene"),
        ("C1=CCCCCC1", "cycloheptene"),
    ],
)
def test_unsubstituted_cyclic_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_single_ene_locant_always_omitted_with_substituent():
    # methyl on the ring carbon three bonds from the double bond: numbering
    # is chosen so the ene is at the (always omitted) position '1' and the
    # substituent gets the lowest achievable locant, '3'.
    assert smiles_to_iupac("C1=CC(C)CCC1") == "3-methylcyclohexene"


def test_substituent_on_ene_carbon_gets_lowest_locant():
    # methyl on a double-bond carbon itself: numbering starts there so both
    # the ene (implicit '1') and methyl land on locant '1', not '2'.
    assert smiles_to_iupac("C1=C(C)CCCC1") == "1-methylcyclohexene"


def test_diene_with_substituent_cites_all_locants():
    assert smiles_to_iupac("C1=CC=C(C)CC1") == "1-methylcyclohexa-1,3-diene"


def test_halogen_substituent():
    assert smiles_to_iupac("C1=CC(Cl)CCC1") == "3-chlorocyclohexene"


def test_unsubstituted_cycloalkyne():
    # cyclooctyne: well-known unambiguous compound name, single ring
    # triple bond so its locant is always omitted (P-14.3.3), same as the
    # single-ene case above.
    assert smiles_to_iupac("C1#CCCCCCC1") == "cyclooctyne"


def test_cycloalkyne_with_substituent():
    # Same locant-tiebreak logic already verified for the ene case
    # (test_single_ene_locant_always_omitted_with_substituent) applied to
    # a triple bond: the yne stays at the omitted '1', methyl gets the
    # lowest achievable locant.
    assert smiles_to_iupac("CC1CCCCCC#C1") == "3-methylcyclooctyne"


def test_mixed_enyne_ring():
    # P-31.1.3.1's own worked example (Blue Book Chapter P-3): the double
    # bond is allocated locant '1' and the triple bond '4' -- lower
    # locants go to the double bond specifically once the combined
    # {1,4}/{1,2}... locant-set choice is tied. Stereodescriptor omitted
    # (E/Z on ring multiple bonds is out of this module's scope).
    assert smiles_to_iupac("C1=CCC#CCCCCCCCCCC1") == "cyclopentadec-1-en-4-yne"


def test_polycyclic_ring_triple_bond_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1#CC2(CCCC1)CCCCC2")


def test_exocyclic_double_bond_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1(=C)CCCCC1")


def test_exocyclic_triple_bond_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1(C#CC)CCCCC1")
