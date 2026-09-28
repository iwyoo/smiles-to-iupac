from smiles_to_iupac import smiles_to_iupac


def test_pyrrolo_ij_quinoline():
    # PubChem CID 12311231, C11H9N; name traced against the Blue Book's
    # pyrrolo[3,2,1-de]acridine worked example (P-25.3.1.3).
    assert smiles_to_iupac("C1C=CC2=CC=CC3=C2N1C=C3") == "4H-pyrrolo[3,2,1-ij]quinoline"
