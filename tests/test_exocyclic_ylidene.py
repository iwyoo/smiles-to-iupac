import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


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


def test_genuinely_unsaturated_ring_still_raises():
    # Ring-internal unsaturation combined with a separate exocyclic
    # double bond -- neither this module's own minimal shape nor
    # `find_cyclic_unsaturated_core`'s ring-only one, so the ring's own
    # unsaturation must still be rejected exactly as before.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCCC1C=C")


def test_two_exocyclic_double_bonds_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(=C)C(=C)C1")


def test_ring_halogen_alongside_exocyclic_ylidene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CCC(=C)CC1")
