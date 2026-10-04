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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A literally-aromatic-flagged ring bridged into a larger von
        # Baeyer polycyclic skeleton is treated as a Kekule cyclohexatriene
        # (P-31.1.4.2). Three real, PubChem-registered structures confirm
        # this exact von Baeyer (not fusion-name) treatment: CIDs 138272,
        # 139650, 582301.
        ("C1CC2CC1c1ccccc12", "tricyclo[6.2.1.0^2,7]undeca-2,4,6-triene"),
        ("C1C2CC1c1ccccc12", "tricyclo[6.1.1.0^2,7]deca-2,4,6-triene"),
        ("C1CC2CCC1c1ccccc12", "tricyclo[6.2.2.0^2,7]dodeca-2,4,6-triene"),
    ],
)
def test_polycyclic_kekule_aromatic_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_two_aromatic_rings_still_is_named():
    # Two independently-aromatic rings bridged by saturated chains (a
    # cyclophane) is a distinct shape this project's dedicated cyclophane
    # modules own -- not a candidate for this single-aromatic-ring
    # mechanism (`_cyclophane.py`'s own tests cover its actual scope).
    assert smiles_to_iupac("c1cc2cc(c1)CCc1ccc(cc1)CC2") == '1(1,3),4(1,4)-dibenzenacyclohexaphane'


def test_polycyclic_extra_unsaturation_outside_aromatic_ring_still_raises():
    # A real double bond coexisting with the aromatic ring, outside the
    # ring itself, is out of this step's scope (#834) -- deferred rather
    # than guessed at.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC2C(C)C1c1ccccc12")
