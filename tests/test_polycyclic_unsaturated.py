import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A real registered tricyclic ene structure (PubChem CID 12478950),
        # whose own IUPACName confirms this exact consecutive-locant form.
        ("C1=CC2CC3CC1CC(C2)C3", "tricyclo[4.3.1.1^3,8]undec-4-ene"),
        # Same skeleton, triple bond instead.
        ("C1#CC2CC3CC1CC(C2)C3", "tricyclo[4.3.1.1^3,8]undec-4-yne"),
        # A halogen substituent coexisting with the ring double bond -- the
        # bond locant still wins the numbering tie-break, mirroring
        # `test_bicyclic_unsaturated.py`'s identical case.
        ("ClC1=CC2CC3CC1CC(C2)C3", "4-chlorotricyclo[4.3.1.1^3,8]undec-4-ene"),
        # A bridgehead alkene (double bond touches locant 1 itself).
        ("C1CC2CC3CC1CC(=C2)C3", "tricyclo[4.3.1.1^3,8]undec-1-ene"),
        # Two double bonds (a tricyclic diene), multiplying prefix 'di'.
        ("C1=CC2CC3CC1C=C(C2)C3", "tricyclo[4.3.1.1^3,8]undeca-1,4-diene"),
    ],
)
def test_polycyclic_unsaturated_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_unsaturated_saturated_case_unaffected():
    assert smiles_to_iupac("C1CC2CC3CC1CC(C2)C3") == "tricyclo[4.3.1.1^3,8]undecane"
