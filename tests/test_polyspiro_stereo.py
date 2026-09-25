import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_linear_polyspiro_off_spiro_stereocenter():
    # Same dispiro[3.2.3^7.2^4]dodecane skeleton as tests/test_polyspiro.py's
    # own halogen-substituent case, stereo-tagged at the chlorine-bearing
    # ring carbon (P-92, mirroring _spiro.py's name_monospiro).
    assert smiles_to_iupac("Cl[C@H]1CCC12CCC1(CCC1)CC2") == "(1S)-1-chlorodispiro[3.2.3^7.2^4]dodecane"


def test_linear_polyspiro_off_spiro_stereocenter_other_configuration():
    assert smiles_to_iupac("Cl[C@@H]1CCC12CCC1(CCC1)CC2") == "(1R)-1-chlorodispiro[3.2.3^7.2^4]dodecane"


def test_branched_polyspiro_off_spiro_stereocenter():
    # trispiro[2.2.2^6.2.3^11.2^3]hexadecane (tests/test_polyspiro.py's own
    # asymmetric-terminal-ring case), stereo-tagged on the cyclobutane
    # terminal ring.
    assert (
        smiles_to_iupac("Cl[C@H]1CCC12CCC1(CC1)CCC1(CC1)CC2")
        == "(12S)-12-chlorotrispiro[2.2.2^6.2.3^11.2^3]hexadecane"
    )


def test_branched_polyspiro_off_spiro_stereocenter_other_configuration():
    assert (
        smiles_to_iupac("Cl[C@@H]1CCC12CCC1(CC1)CCC1(CC1)CC2")
        == "(12R)-12-chlorotrispiro[2.2.2^6.2.3^11.2^3]hexadecane"
    )


def test_plain_linear_polyspiro_without_stereo_still_resolves():
    assert smiles_to_iupac("C1CCC12CCC3(CC2)CCC3") == "dispiro[3.2.3^7.2^4]dodecane"


def test_substituent_branch_stereocenter_raises():
    # A stereocenter on the ring's substituent branch rather than the
    # polyspiro ring skeleton itself (P-92) -- matches _spiro.py's own
    # existing behavior for the identical shape (confirmed directly:
    # smiles_to_iupac("C1CCCC12CCCC2[C@@H](C)CC") also raises).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("F[C@@H](Cl)C1CCC12CCC1(CCC1)CC2")
