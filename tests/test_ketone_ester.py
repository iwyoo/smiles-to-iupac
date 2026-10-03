import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_worked_example_methyl_acetoacetate():
    # methyl acetoacetate (PubChem CID 7757, CAS 105-45-3): IUPAC name
    # 'methyl 3-oxobutanoate', confirmed via PubChem/NIST WebBook.
    assert smiles_to_iupac("CC(=O)CC(=O)OC") == "methyl 3-oxobutanoate"


def test_ethyl_acetoacetate():
    # A single-axis variant of the worked example above (methyl -> ethyl
    # ester), same demotion mechanism.
    assert smiles_to_iupac("CCOC(=O)CC(=O)C") == "ethyl 3-oxobutanoate"


def test_two_ketones():
    assert smiles_to_iupac("CC(=O)CC(=O)CC(=O)OC") == "methyl 3,5-dioxohexanoate"


def test_halogen_substituent():
    assert smiles_to_iupac("CC(=O)C(Cl)C(=O)OC") == "methyl 2-chloro-3-oxobutanoate"


def test_plain_ester_still_works():
    assert smiles_to_iupac("CCC(=O)OC") == "methyl propanoate"


def test_plain_ketone_still_works():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


def test_unsaturated_chain():
    assert smiles_to_iupac("C=CC(=O)CC(=O)OC") == "methyl 3-oxopent-4-enoate"


def test_ring():
    assert smiles_to_iupac("COC(=O)C1CCC1=O") == "methyl 2-oxocyclobutane-1-carboxylate"


def test_hydroxyl_coexistence():
    assert smiles_to_iupac("OCC(=O)CC(=O)OC") == "methyl 4-hydroxy-3-oxobutanoate"


def test_acyl_carbon_off_longest_chain():
    assert (
        smiles_to_iupac("CCCCC(=O)C(CC(C)C)C(CC(C)C)C(=O)OC")
        == "methyl 2,3-bis(2-methylpropyl)-4-oxooctanoate"
    )


def test_phenyl_chain_ketone_ester():
    # A plain, unsubstituted benzene ring on the acyl chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring
    # `_aldehyde_ketone.py`'s phenyl-chain path): the ring is cited as a
    # "phenyl" substituent prefix alongside the demoted ketone's "oxo"
    # prefix. PubChem PUG REST: "methyl 3-oxo-4-phenylbutanoate"/
    # "methyl 2-oxo-4-phenylbutanoate".
    assert smiles_to_iupac("c1ccccc1CC(=O)CC(=O)OC") == "methyl 3-oxo-4-phenylbutanoate"
    assert smiles_to_iupac("c1ccccc1CCC(=O)C(=O)OC") == "methyl 2-oxo-4-phenylbutanoate"


def test_phenyl_chain_ketone_ester_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(=O)CC(=O)OC") == "methyl 4-(2-ethenylphenyl)-3-oxobutanoate"
