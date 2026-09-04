import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanediazonium():
    # PubChem structure match: "methanediazonium".
    assert smiles_to_iupac("C[N+]#N") == "methanediazonium"


def test_ethanediazonium():
    # PubChem structure match: "ethanediazonium".
    assert smiles_to_iupac("CC[N+]#N") == "ethanediazonium"


def test_propane_1_diazonium():
    # PubChem structure match: "propane-1-diazonium".
    assert smiles_to_iupac("CCC[N+]#N") == "propane-1-diazonium"


def test_propane_2_diazonium():
    # PubChem structure match: "propane-2-diazonium" (a branched chain).
    assert smiles_to_iupac("CC(C)[N+]#N") == "propane-2-diazonium"


def test_chloroethane_diazonium():
    # Locant cited once a substituent is present on a 2-carbon chain, the
    # same project-wide convention as '2-chloroethane-1-selenol'
    # (PubChem's own generated name, "2-chloroethanediazonium", omits it).
    assert smiles_to_iupac("ClCC[N+]#N") == "2-chloroethane-1-diazonium"


def test_pent_4_ene_1_diazonium():
    assert smiles_to_iupac("C=CCCC[N+]#N") == "pent-4-ene-1-diazonium"


def test_two_diazonium_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#[N+]CC[N+]#N")


def test_cyclohexanediazonium():
    # PubChem structure match: "cyclohexanediazonium".
    assert smiles_to_iupac("C1CCCCC1[N+]#N") == "cyclohexanediazonium"


def test_cyclopentanediazonium():
    # PubChem structure match: "cyclopentanediazonium".
    assert smiles_to_iupac("C1CCCC1[N+]#N") == "cyclopentanediazonium"


def test_cyclopropanediazonium():
    assert smiles_to_iupac("C1CC1[N+]#N") == "cyclopropanediazonium"


def test_substituted_ring_diazonium_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCCCC1[N+]#N")


def test_halogen_substituted_ring_diazonium_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CCCCC1[N+]#N")


def test_benzenediazonium():
    # -N#N+ directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1[N+]#N") == "benzenediazonium"  # CID 9718


def test_polycyclic_diazonium_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12(CCC(CC1)CC2)[N+]#N")


def test_diazonium_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[N+]#N")


def test_ammonium_not_confused_with_diazonium():
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"


def test_phenyl_chain_diazonium():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_sulfonic_acid.py`'s
    # phenyl-chain path): the ring is cited as a "phenyl" substituent
    # prefix. PubChem PUG REST: "2-phenylethanediazonium" (no locant) /
    # "3-phenylpropane-1-diazonium" -- the two-carbon case keeps its
    # locant here, same known, out-of-scope-here limitation already
    # established elsewhere in this project (e.g. `_thiol.py`'s own
    # phenyl-chain two-carbon case).
    assert smiles_to_iupac("c1ccccc1CC[N+]#N") == "2-phenylethane-1-diazonium"
    assert smiles_to_iupac("c1ccccc1CCC[N+]#N") == "3-phenylpropane-1-diazonium"


def test_substituted_benzenediazonium():
    # A substituent on a different ring atom than the diazonium group:
    # the ring is renumbered to give the -N#N+ carbon locant 1 (never
    # cited), then the other substituent gets the lowest remaining
    # locant, cross-checked against PubChem PUG REST.
    assert smiles_to_iupac("Cc1ccccc1[N+]#N") == "2-methylbenzenediazonium"  # CID 192837


def test_dimethyl_benzenediazonium():
    # Two substituents on the ring, lowest-locant-set + alphabetical
    # citation order (P-14.5.2), structurally verified (no matching
    # PubChem CID for this exact regiochemistry).
    assert smiles_to_iupac("Cc1ccc(C)cc1[N+]#N") == "2,5-dimethylbenzenediazonium"


def test_phenyl_substituted_benzene_ring_diazonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[N+]#N")


def test_phenyl_chain_diazonium_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[N+]#N")
