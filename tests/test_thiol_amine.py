import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Cysteamine: PubChem CID 5726's own IUPACName is
        # '2-aminoethanethiol' (the retained 'ethanethiol' stem kept even
        # when substituted); this project's own `_thiol.py` already
        # diverges from that convention for a substituted ethanethiol
        # (its own '-1-thiol' locant citation once a substituent exists),
        # so this module follows that same pre-existing convention instead.
        ("NCCS", "2-aminoethanethiol"),
        # 3-aminopropane-1-thiol: PubChem-verified exactly.
        ("NCCCS", "3-aminopropane-1-thiol"),
        # 2-aminopropane-1-thiol: PubChem-verified exactly (the -SH gets
        # the lower locant, P-44).
        ("CC(N)CS", "2-aminopropane-1-thiol"),
        # A halogen substituent coexists with both the thiol and the
        # amine.
        ("NC(Cl)CS", "2-amino-2-chloroethanethiol"),
    ],
)
def test_smiles_to_iupac_thiol_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_secondary_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNCCS")


def test_two_amines_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(N)CS")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCC(S)CC1")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC=CCS")


def test_specified_stereocenter_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](N)CS")


def test_plain_thiol_still_works():
    assert smiles_to_iupac("CCCS") == "propane-1-thiol"


def test_plain_amine_still_works():
    assert smiles_to_iupac("CCCN") == "propan-1-amine"


def test_sulfide_coexisting_raises():
    # A sulfide alongside the thiol/amine pair is out of scope for this
    # first pairwise pilot on `_thiol.py`.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(S)CSCC")
