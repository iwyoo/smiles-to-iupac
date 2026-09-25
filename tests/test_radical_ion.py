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


def test_methyloxidaniumyl():
    # C[OH+], the oxidaniumyl radical cation (P-75.3.2): methyloxidanium
    # (already-supported oxonium cation) minus one radical hydrogen.
    assert smiles_to_iupac("C[OH+]") == "methyloxidaniumyl"


def test_ethyloxidaniumyl():
    assert smiles_to_iupac("CC[OH+]") == "ethyloxidaniumyl"


def test_ethylsulfaniumyl():
    # note: _sulfonium.py's own established output is "sulfanium", not
    # "sulfonium".
    assert smiles_to_iupac("CC[SH+]") == "ethylsulfaniumyl"


def test_methanaminyliumyl():
    # C[N+], the aminyliumyl radical cation (P-73.2.3.2's 'ylium' cation
    # plus P-75.3.1's radical 'yl', stacked): needs a double
    # reconstruction all the way back to neutral methanamine, since the
    # intermediate 'methanaminylium' cation is itself still radical-
    # carrying.
    assert smiles_to_iupac("C[N+]") == "methanaminyliumyl"


def test_ethaniminyliumyl():
    assert smiles_to_iupac("CC=[N+]") == "ethaniminyliumyl"


def test_ethanamidyliumyl():
    # matches the Blue Book's own literal worked example shape
    # `acetamidyliumyl (PIN)`, tmp/bluebook/P7.txt ~3478-3496 -- this
    # project's own established never-special-case-acetic/formic policy
    # gives the systematic "ethanamidyliumyl" instead of the retained
    # "acetamidyliumyl".
    assert smiles_to_iupac("CC(=O)[N+]") == "ethanamidyliumyl"


def test_oxonium_still_resolves():
    # a sanity check that the new radical-ion dispatch doesn't misfire on
    # a plain, non-radical oxonium cation.
    assert smiles_to_iupac("C[OH2+]") == "methyloxidanium"


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
