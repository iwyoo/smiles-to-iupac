from smiles_to_iupac import smiles_to_iupac


def test_single_off_spiro_stereocenter():
    assert smiles_to_iupac("C[C@H]1CCC[C@]2(C1)CCCC2") == "(7S)-7-methylspiro[4.5]decane"


def test_single_off_spiro_stereocenter_other_configuration():
    assert smiles_to_iupac("C[C@@H]1CCC[C@]2(C1)CCCC2") == "(7R)-7-methylspiro[4.5]decane"


def test_two_off_spiro_stereocenters_different_rings():
    assert smiles_to_iupac("C[C@H]1CCC[C@]2(C1)C[C@H](C)CCC2") == "(2S,6R,8R)-2,8-dimethylspiro[5.5]undecane"


def test_pseudoasymmetric_spiro_atom_is_cited_in_lower_case():
    assert smiles_to_iupac("C[C@H]1CC[C@]2(C1)CC[C@@H](C)CC2") == "(2S,5s,8S)-2,8-dimethylspiro[4.5]decane"


def test_plain_spiro_without_stereo_still_resolves():
    assert smiles_to_iupac("C1CCCC12CCCCC2") == "spiro[4.5]decane"


def test_substituted_spiro_without_stereo_still_resolves():
    assert smiles_to_iupac("CC1CCCC12CCCCC2") == "1-methylspiro[4.5]decane"
