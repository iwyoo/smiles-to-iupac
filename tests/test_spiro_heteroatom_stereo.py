from smiles_to_iupac import smiles_to_iupac


def test_oxaspiro_single_off_spiro_stereocenter():
    assert smiles_to_iupac("C[C@H]1CCCC12CCOCC2") == "(1S)-1-methyl-8-oxaspiro[4.5]decane"


def test_oxaspiro_other_configuration():
    assert smiles_to_iupac("C[C@@H]1CCCC12CCOCC2") == "(1R)-1-methyl-8-oxaspiro[4.5]decane"


def test_azaspiro_single_off_spiro_stereocenter():
    assert smiles_to_iupac("C[C@H]1CCCC12CCNCC2") == "(1S)-1-methyl-8-azaspiro[4.5]decane"


def test_stereocenter_at_the_spiro_atom_itself():
    assert smiles_to_iupac("C[C@H]1CCC[C@]2(C1)CCOC2") == "(5S,7S)-7-methyl-2-oxaspiro[4.5]decane"


def test_plain_oxaspiro_without_stereo_still_resolves():
    assert smiles_to_iupac("C1CCCC12CCOCC2") == "8-oxaspiro[4.5]decane"
