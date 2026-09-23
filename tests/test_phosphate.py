import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-67.1.3.2), symmetric trialkyl/aryl
        # phosphate esters -- verified against each structure's own
        # PubChem-registered IUPAC name.
        ("COP(=O)(OC)OC", "trimethyl phosphate"),  # CID 10541
        ("CCOP(=O)(OCC)OCC", "triethyl phosphate"),  # CID 6579
        ("CCCOP(=O)(OCCC)OCCC", "tripropyl phosphate"),  # CID 12659
        ("CCCCOP(=O)(OCCCC)OCCCC", "tributyl phosphate"),  # CID 6503
        ("c1ccc(OP(=O)(Oc2ccccc2)Oc2ccccc2)cc1", "triphenyl phosphate"),  # CID 6987
        ("CC(C)OP(=O)(OC(C)C)OC(C)C", "tripropan-2-yl phosphate"),  # CID 31289
        ("CCCCCCCCOP(=O)(OCCCCCCCC)OCCCCCCCC", "trioctyl phosphate"),  # CID 31261
        # A branch name starting with its own locant digit needs 'tris'
        # plus enclosing marks, unlike the internally-locanted
        # 'tripropan-2-yl' case above (confirmed against both real
        # PubChem structures, see the PR this was introduced in).
        ("ClCCOP(=O)(OCCCl)OCCCl", "tris(2-chloroethyl) phosphate"),  # CID 6588
        (
            "CCCCC(CC)COP(=O)(OCC(CC)CCCC)OCC(CC)CCCC",
            "tris(2-ethylhexyl) phosphate",
        ),  # CID 31284
    ],
)
def test_phosphate_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_phosphate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOP(=O)(OC)OC")


def test_partial_hydrogen_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COP(=O)(OC)O")


def test_phosphite_ester_unaffected():
    # No P=O bond -- routed to a different module entirely
    # (`_phosphane.py`), not this one; still unsupported, but for a
    # different reason (a future milestone step, not this pilot).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOP(OCC)OCC")


def test_phosphonic_acid_unaffected():
    assert smiles_to_iupac("CP(=O)(O)O") == "methylphosphonic acid"


def test_phosphanone_unaffected():
    assert smiles_to_iupac("CCP(=O)(CC)CC") == "triethyl-λ5-phosphanone"
