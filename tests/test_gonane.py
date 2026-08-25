from smiles_to_iupac import smiles_to_iupac


def test_gonane():
    # Gonane's PIN per the 1989 IUPAC steroid nomenclature (Rule 2.1,
    # https://iupac.qmul.ac.uk/steroid/3S02a.html): "the parent tetracyclic
    # hydrocarbon without methyl groups at C-10 and C-13 and without a side
    # chain at C-17". Structure (C17H28) cross-checked against PubChem CID
    # 6857523's ConnectivitySMILES, which canonicalizes identically.
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("C1CCC2C(C1)CCC3C2CCC4C3CCC4") == "gonane"


def test_stereo_specified_gonane_falls_through_to_von_baeyer():
    # Gonane's own definition (Rule 2.1) doesn't specify ring-fusion
    # stereochemistry, so this module deliberately leaves stereo-marked
    # input unmatched, falling through to the general von Baeyer engine
    # (`_polycyclic.py`), which already names this exact skeleton.
    assert (
        smiles_to_iupac("C1CCC2[C@H](C1)CCC3C2CCC4C3CCC4")
        == "tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_methylated_gonane_falls_through_to_von_baeyer():
    # Angular methyls at C-10/C-13 (present from androstane up) are outside
    # this module's scope by definition -- falls through to the general
    # von Baeyer engine instead of being (mis)recognized as gonane.
    assert (
        smiles_to_iupac("CC12CCCC1(C)CCC1C2CCC2C1CCCC2")
        == "11,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
