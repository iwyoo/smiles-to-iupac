import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanaminiumyl():
    # C[NH2+], the aminiumyl radical cation (P-75.3.1): methanaminium
    # (already-supported ammonium cation) minus one radical hydrogen.
    assert smiles_to_iupac("C[NH2+]") == "methanaminiumyl"


def test_ethanaminiumyl():
    assert smiles_to_iupac("CC[NH2+]") == "ethanaminiumyl"


def test_aniliniumyl():
    # aromatic ring case, mirrors the Blue Book's own literal worked
    # example `benzenaminiumyl (PIN)`, tmp/bluebook/P7.txt ~3478-3496 --
    # this project's own ammonium naming uses the retained "anilinium"
    # instead of the systematic "benzenaminium", consistent with
    # name_ammonium's own established output.
    assert smiles_to_iupac("c1ccccc1[NH2+]") == "aniliniumyl"


def test_secondary_aminiumyl():
    assert smiles_to_iupac("C[NH+]C") == "N-methylmethanaminiumyl"


def test_ammonium_still_resolves():
    # a sanity check that the new radical-ion dispatch doesn't misfire on
    # a plain, non-radical ammonium cation.
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"


def test_plain_radical_still_resolves():
    # a sanity check that the new radical-ion dispatch doesn't misfire on
    # a plain, uncharged radical.
    assert smiles_to_iupac("C[CH2]") == "ethyl"


def test_coexisting_second_charge_unsupported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[NH2+]CC[O-]")
