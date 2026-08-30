import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanethiol():
    # '-thiol' suffix construction mirrors '-ol' (_alcohol.py, see
    # test_alcohol.py, already cross-checked against PubChem) with the same
    # locant-citation rules; 'methanethiol'/'ethanethiol' are also the
    # well-known common names of these compounds.
    assert smiles_to_iupac("SC") == "methanethiol"


def test_ethanethiol():
    assert smiles_to_iupac("SCC") == "ethanethiol"


def test_propane_1_thiol():
    assert smiles_to_iupac("SCCC") == "propane-1-thiol"


def test_propane_2_thiol():
    assert smiles_to_iupac("CC(S)C") == "propane-2-thiol"


def test_butane_1_thiol_with_chloro_substituent():
    assert smiles_to_iupac("SCCCCCl") == "4-chlorobutane-1-thiol"


def test_pent_4_ene_1_thiol():
    assert smiles_to_iupac("SCCCC=C") == "pent-4-ene-1-thiol"


def test_ethane_1_2_dithiol():
    # PubChem CID 10902.
    assert smiles_to_iupac("SCCS") == "ethane-1,2-dithiol"


def test_propane_1_3_dithiol():
    # PubChem CID 8013.
    assert smiles_to_iupac("SCCCS") == "propane-1,3-dithiol"


def test_propane_1_2_dithiol():
    # PubChem CID 61217.
    assert smiles_to_iupac("CC(S)CS") == "propane-1,2-dithiol"


def test_butane_1_4_dithiol():
    # PubChem CID 79148.
    assert smiles_to_iupac("SCCCCS") == "butane-1,4-dithiol"


def test_cyclohexanethiol():
    # PubChem CID 15290.
    assert smiles_to_iupac("SC1CCCCC1") == "cyclohexanethiol"


def test_cyclohexane_1_2_dithiol():
    # PubChem CID 5242646.
    assert smiles_to_iupac("SC1CCCCC1S") == "cyclohexane-1,2-dithiol"


def test_2_methylcyclohexane_1_thiol():
    # PubChem CID 519947.
    assert smiles_to_iupac("CC1CCCCC1S") == "2-methylcyclohexane-1-thiol"


def test_cyclopentanethiol():
    # PubChem CID 15510.
    assert smiles_to_iupac("SC1CCCC1") == "cyclopentanethiol"


def test_polycyclic_thiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CC2CCC1CC2")


def test_unsaturated_ring_thiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC=C1")


def test_thiol_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCC1CCCCC1")


def test_thiol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCCO")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): unlike
        # `_sulfinic_acid.py`, a thiol's -SH sulfur is monovalent (bonded
        # only to carbon and H) and can never itself be a stereocenter, so
        # this mirrors `_carboxylic_acid.py`/`_aldehyde.py`/`_ketone.py`/
        # `_sulfonic_acid.py` cleanly (CIP computed entirely by RDKit's
        # `rdCIPLabeler`). PubChem CID 444090.
        ("CC[C@@H](C)S", "(2R)-butane-2-thiol"),
        ("CC[C@H](C)S", "(2S)-butane-2-thiol"),
    ],
)
def test_acyclic_thiol_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_thiol_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_alcohol.py`'s `_name_cyclic_alcohol`. PubChem CID 22211639 (name
    # carries a redundant 'cis-' relative descriptor this project drops
    # once full R/S is given, same policy as the existing ring-alcohol
    # stereocenter task).
    assert smiles_to_iupac("S[C@H]1CCCC[C@H]1C") == "(1S,2R)-2-methylcyclohexane-1-thiol"


def test_thiol_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)S") == "butane-2-thiol"


def test_thiol_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S[C@H]1CCCCC1Cl")
