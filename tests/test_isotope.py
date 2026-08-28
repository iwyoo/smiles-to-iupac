import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_deuteromethane_name():
    # methane-d1, PubChem CID 12669 (structure-verified via
    # ConnectivitySMILES "[2H]C"; PubChem's own auto-generated name
    # "deuteriomethane" isn't in PIN format, so the name itself is
    # justified by the Blue Book's own worked example, P-82.2.1:
    # "CH3-2H -> (2H1)methane (PIN)").
    assert smiles_to_iupac("[2H]C") == "(2H1)methane"


def test_plain_methane_unaffected():
    assert smiles_to_iupac("C") == "methane"


def test_ethane_deuterium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]CC")


def test_tritium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[3H]C")


def test_dideuteromethane_name():
    # CD2H2, PubChem structure-verified via ConnectivitySMILES (CID
    # 259772, "[2H]C([2H])" core matches); name is the Blue Book's own
    # P-82.2.1 subscript-counting rule applied to the confirmed single-
    # deuterium worked example ("(2H1)methane (PIN)") -- the count is
    # "always specified as a right subscript ... even in case of
    # monosubstitution", so 2 atoms substituted is '(2H2)'.
    assert smiles_to_iupac("[2H]C([2H])") == "(2H2)methane"


def test_tetradeuteromethane_name():
    # CD4, PubChem CID 24246 (structure-verified via ConnectivitySMILES
    # "[2H]C([2H])([2H])[2H]"); name derived the same way as the dideuterio
    # case above.
    assert smiles_to_iupac("[2H]C([2H])([2H])[2H]") == "(2H4)methane"


def test_dichlorodideuteromethane_name():
    # CD2Cl2, Blue Book P-82.2.1 worked example: "dichloro(2H2)methane
    # (PIN)". Structure cross-checked against PubChem CID 160586
    # (ConnectivitySMILES "[2H]C([2H])(Cl)Cl"; PubChem's own auto-generated
    # name "dichloro(dideuterio)methane" isn't PIN format).
    assert smiles_to_iupac("[2H]C([2H])(Cl)Cl") == "dichloro(2H2)methane"


def test_carbon_14_methane_name():
    # (14C)methane, Blue Book P-82.2.1 worked example: "(14C)methane
    # (PIN)". Structure cross-checked against PubChem CID 26873
    # (MolecularFormula CH4, isotopically labeled 14C).
    assert smiles_to_iupac("[14CH4]") == "(14C)methane"


def test_trichloro_carbon_12_methane_name():
    # CHCl3 with the carbon explicitly 12C, Blue Book P-82.2.1 worked
    # example: "trichloro(12C)methane (PIN) / (12C)chloroform". 12C is
    # carbon's most abundant natural isotope, so it has no distinct
    # PubChem CID to structure-verify against -- the name itself is the
    # Blue Book's own cited worked example.
    assert smiles_to_iupac("[12CH](Cl)(Cl)Cl") == "trichloro(12C)methane"


def test_deuterium_and_carbon_isotope_together_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H][13CH3]")


def test_isotopically_labeled_halogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]C([37Cl])")
