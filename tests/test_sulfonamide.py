import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanesulfonamide():
    # PubChem structure match: "methanesulfonamide".
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_ethanesulfonamide():
    # PubChem structure match: "ethanesulfonamide".
    assert smiles_to_iupac("CCS(=O)(=O)N") == "ethanesulfonamide"


def test_propane_1_sulfonamide():
    # PubChem structure match: "propane-1-sulfonamide".
    assert smiles_to_iupac("CCCS(=O)(=O)N") == "propane-1-sulfonamide"


def test_propane_2_sulfonamide():
    assert smiles_to_iupac("CC(S(=O)(=O)N)C") == "propane-2-sulfonamide"


def test_chlorobutanesulfonamide():
    assert smiles_to_iupac("ClCCCCS(=O)(=O)N") == "4-chlorobutane-1-sulfonamide"


def test_2_chloroethane_1_sulfonamide():
    # Locant cited even on a 2-carbon chain once a substituent (here, the
    # halogen) is present -- same project-wide convention as
    # '2-chloroethane-1-selenol'/'2-chloroethane-1-thiol' (PubChem's own
    # generated name omits the locant here; that divergence is accepted
    # project-wide, see test_selenol.py's identical case).
    assert smiles_to_iupac("ClCCS(=O)(=O)N") == "2-chloroethane-1-sulfonamide"


def test_pent_4_ene_1_sulfonamide():
    assert smiles_to_iupac("C=CCCCS(=O)(=O)N") == "pent-4-ene-1-sulfonamide"


def test_disulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)CCS(=O)(=O)N")


def test_cyclohexanesulfonamide():
    # PubChem structure match: "cyclohexanesulfonamide".
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1") == "cyclohexanesulfonamide"


def test_2_methylcyclohexane_1_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1C") == "2-methylcyclohexane-1-sulfonamide"


def test_cyclopentanesulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC1") == "cyclopentanesulfonamide"


def test_2_chlorocyclohexane_1_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfonamide"


def test_polycyclic_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CCCC=C1")


def test_sulfonamide_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)CC1CCCCC1")


def test_sulfonamide_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)CCO")


def test_n_substituted_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)(=O)NC")
