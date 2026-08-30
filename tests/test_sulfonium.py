import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_sulfanium():
    # PubChem structure match: "sulfanium" (SH3+).
    assert smiles_to_iupac("[SH3+]") == "sulfanium"


def test_methylsulfanium():
    # PubChem structure match: "methylsulfanium".
    assert smiles_to_iupac("C[SH2+]") == "methylsulfanium"


def test_ethylsulfanium():
    # PubChem structure match: "ethylsulfanium".
    assert smiles_to_iupac("CC[SH2+]") == "ethylsulfanium"


def test_dimethylsulfanium():
    # PubChem structure match: "dimethylsulfanium".
    assert smiles_to_iupac("C[SH+]C") == "dimethylsulfanium"


def test_trimethylsulfanium():
    # PubChem structure match: "trimethylsulfanium".
    assert smiles_to_iupac("C[S+](C)C") == "trimethylsulfanium"


def test_branched_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)[SH2+]")


def test_ring_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SH2+]C1CCCCC1")


def test_thiol_not_confused_with_sulfonium():
    assert smiles_to_iupac("CS") == "methanethiol"


def test_phosphonium_not_confused_with_sulfonium():
    assert smiles_to_iupac("C[PH3+]") == "methylphosphanium"


def test_sulfonium_specified_stereocenter_raises():
    # The sulfonium sulfur (three distinct alkyl substituents) is itself a
    # genuine, configurationally stable stereocenter -- this project has no
    # established way to cite a heteroatom-centered stereodescriptor, so a
    # specified one must raise rather than silently collapsing both
    # configurations to the same name (P-92).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[S@+](CC)CCC")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[S@@+](CC)CCC")


def test_sulfonium_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("C[S+](CC)CCC") == "ethyl(methyl)(propyl)sulfanium"
