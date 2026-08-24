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


def test_ring_triple_bond_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1#CCCCC1")


def test_exocyclic_double_bond_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1(=C)CCCCC1")
