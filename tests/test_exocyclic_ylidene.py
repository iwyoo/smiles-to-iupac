import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 14502, matches its own IUPACName exactly -- the
        # motivating case (#1090): a saturated ring bearing a single
        # exocyclic '=CH2' substituent was previously misclassified as
        # ring unsaturation.
        ("C1CCC(=C)CC1", "methylidenecyclohexane"),
        ("C1CCC(=CC)CC1", "ethylidenecyclohexane"),
        ("C1CCC(=C(C)C)CC1", "(propan-2-ylidene)cyclohexane"),
        ("C1CCCC1=C", "methylidenecyclopentane"),
    ],
)
def test_exocyclic_ylidene_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_double_bond_alongside_ethenyl_prefix():
    assert smiles_to_iupac("C1=CCCCC1C=C") == "3-ethenylcyclohexene"


def test_two_exocyclic_double_bonds():
    assert smiles_to_iupac("C1CCC(=C)C(=C)C1") == "1,2-dimethylidenecyclohexane"


def test_ring_halogen_alongside_exocyclic_ylidene():
    assert smiles_to_iupac("ClC1CCC(=C)CC1") == "1-chloro-4-methylidenecyclohexane"
