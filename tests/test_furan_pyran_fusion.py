from smiles_to_iupac import smiles_to_iupac


def test_furo_3_2_b_pyran():
    # Blue Book's own worked example, P-25.3.2.4(c), `tmp/bluebook/P2.pdf`
    # page 75: "2H-furo[3,2-b]pyran (PIN)".
    assert smiles_to_iupac("C1=COC2=CCOC2=C1") == "2H-furo[3,2-b]pyran"
