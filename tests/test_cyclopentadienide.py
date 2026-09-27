from smiles_to_iupac import smiles_to_iupac


def test_cyclopentadienide_name():
    # Blue Book P-72.2.2.1's own worked example: the cyclopentadienyl
    # anion names as 'cyclopenta-2,4-dien-1-ide (PIN)', with the anion
    # center fixed at locant 1.
    assert smiles_to_iupac("[CH-]1C=CC=C1") == "cyclopenta-2,4-dien-1-ide"


def test_cyclopentadiene_neutral_parent_unaffected():
    assert smiles_to_iupac("C1=CC=CC1") == "cyclopenta-1,3-diene"
