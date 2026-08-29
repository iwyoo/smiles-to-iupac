import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-sulfinamide' mirrors '-sulfonamide' (see test_sulfonamide.py)
        # with one fewer oxygen; same locant rules as '-sulfinic acid'.
        ("CS(=O)N", "methanesulfinamide"),
        ("CCS(=O)N", "ethanesulfinamide"),
        ("CCCS(=O)N", "propane-1-sulfinamide"),
        ("CC(S(=O)N)C", "propane-2-sulfinamide"),
        ("CCCCS(=O)N", "butane-1-sulfinamide"),
    ],
)
def test_saturated_sulfinamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_sulfinamide():
    assert smiles_to_iupac("C=CCS(=O)N") == "prop-2-ene-1-sulfinamide"


def test_halogen_substituent():
    # PubChem's own generated name ("1-chloroethanesulfinamide") omits the
    # locant; this project cites it once a substituent is present on a
    # 2-carbon chain, the same accepted divergence as
    # '2-chloroethane-1-selenol' (see test_selenol.py).
    assert smiles_to_iupac("CC(Cl)S(=O)N") == "1-chloroethane-1-sulfinamide"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C(S(=O)N)C")


def test_two_sulfinamides_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CS(=O)N")


def test_sulfonamide_not_confused_with_sulfinamide():
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_sulfinic_acid_not_confused_with_sulfinamide():
    assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"


def test_cyclohexanesulfinamide():
    # PubChem structure match: "cyclohexanesulfinamide".
    assert smiles_to_iupac("O=S(N)C1CCCCC1") == "cyclohexanesulfinamide"


def test_2_methylcyclohexane_1_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCCC1C") == "2-methylcyclohexane-1-sulfinamide"


def test_cyclopentanesulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCC1") == "cyclopentanesulfinamide"


def test_2_chlorocyclohexane_1_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfinamide"


def test_polycyclic_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(N)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(N)C1CCCC=C1")


def test_sulfinamide_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CC1CCCCC1")


def test_sulfinamide_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CCO")


def test_n_substituted_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)NC")
