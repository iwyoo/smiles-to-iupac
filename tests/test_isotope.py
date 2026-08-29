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


def test_carbon_14_butane_name():
    # Blue Book P-82.5.2 worked example (numbering-priority rule): "(2-14C)butane
    # (PIN) [not (3-14C)butane]" -- the internal chain carbon adjacent to a
    # terminus must be numbered '2' (not '3'), confirming a locant IS needed
    # (and how it's chosen) for an unhalogenated multi-carbon chain, unlike
    # methane's own never-locanted case.
    assert smiles_to_iupac("C[14CH2]CC") == "(2-14C)butane"


def test_trifluoro_deuterio_ethane_name():
    # Blue Book P-82.6.2 worked example: "1,1,1-trifluoro(2-2H1)ethane (PIN)".
    assert smiles_to_iupac("FC(F)(F)C[2H]") == "1,1,1-trifluoro(2-2H1)ethane"


def test_carbon_isotope_propane_name():
    # Generalizes the confirmed butane rule (P-82.5.2) to a shorter
    # unhalogenated chain: propane's central carbon is likewise not
    # numbering-direction-symmetric with a terminal carbon, so its locant
    # is cited the same way.
    assert smiles_to_iupac("C[14CH2]C") == "(2-14C)propane"


def test_plain_ethane_deuterium_raises():
    # An unhalogenated 2-carbon chain with a single isotopic modification
    # isn't covered by any confirmed Blue Book worked example this module
    # has verified -- P-82.6.1.1's omission rule for this exact shape is
    # plausible but unconfirmed, so this stays a deliberate rejection
    # rather than a guessed name.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[2H]CC")


def test_single_halogen_ethane_deuterium_raises():
    # A 2-carbon chain with exactly one halogen substituent alongside an
    # isotopic modification is likewise unconfirmed.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("FCC[2H]")


def test_deuterium_locant_set_two_positions():
    # Deuterium spread across two chain positions (a locant *set*, mirroring
    # Blue Book's own '(1,1,1,3,3-2H5)pentan-2-one' shape). PubChem
    # structure match: "1,2-dideuterioethane" (its own systematic
    # 'deuterio' prefix style, not the Blue Book nuclide descriptor).
    assert smiles_to_iupac("[2H]CC[2H]") == "(1,2-2H2)ethane"


def test_deuterium_locant_set_four_atoms_two_positions():
    # PubChem structure match: "1,1,2,2-tetradeuterioethane".
    assert smiles_to_iupac("[2H]C([2H])C([2H])[2H]") == "(1,1,2,2-2H4)ethane"


def test_deuterium_locant_set_three_carbon_chain():
    # PubChem structure match: "1,1,3,3-tetradeuteriopropane".
    assert smiles_to_iupac("[2H]C([2H])CC([2H])[2H]") == "(1,1,3,3-2H4)propane"


def test_deuterium_locant_set_with_halogen():
    # Isotope and halogen locants are minimized together as one combined
    # series (P-82.5.2) when picking chain-numbering direction; here that
    # unambiguously favors citing the chlorine at position 3, not 1.
    assert smiles_to_iupac("[2H]C([2H])C([2H])CCl") == "3-chloro(1,1,2-2H3)propane"


def test_branched_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C[2H]")
