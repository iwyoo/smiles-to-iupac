from smiles_to_iupac import smiles_to_iupac


def test_oxadispiro_off_spiro_stereocenter():
    # 6-oxadispiro[4.1.4^7.2^5]tridecane (tests/test_polyspiro_heteroatom.py's
    # own halogen-substituent case, moved one position over so it's the sole
    # potential stereocenter), stereo-tagged at the chlorine-bearing ring
    # carbon.
    assert (
        smiles_to_iupac("Cl[C@@H]1CC2(CCCC2)OC12CCCC2")
        == "(12R)-12-chloro-6-oxadispiro[4.1.4^7.2^5]tridecane"
    )


def test_oxadispiro_off_spiro_stereocenter_other_configuration():
    assert (
        smiles_to_iupac("Cl[C@H]1CC2(CCCC2)OC12CCCC2")
        == "(12S)-12-chloro-6-oxadispiro[4.1.4^7.2^5]tridecane"
    )


def test_oxadispiro_stereocenter_at_the_spiro_atom_itself():
    # The oxygen's position (adjacent to one spiro atom but not the other)
    # breaks the local symmetry that would otherwise make that spiro atom
    # achiral -- same dependent-stereocenter shape as _spiro_heteroatom.py's
    # own "stereocenter at the spiro atom itself" test.
    assert (
        smiles_to_iupac("C[C@H]1CCC[C@]2(COC3(CCC3)C2)C1")
        == "(6S,8S)-8-methyl-13-oxadispiro[3.1.5^6.2^4]tridecane"
    )


def test_plain_oxadispiro_without_stereo_still_resolves():
    assert smiles_to_iupac("C1CCC2(C1)CCC3(O2)CCCC3") == "6-oxadispiro[4.1.4^7.2^5]tridecane"
