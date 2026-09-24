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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # All three R groups distinct: three separate un-prefixed words,
        # alphabetic order. PubChem CID 12494385/12494386.
        ("CCOP(=O)(OC)Oc1ccccc1", "ethyl methyl phenyl phosphate"),
        # Two identical R groups + one distinct: one multiplied word plus
        # one plain word, ordered by each word's own base name (ignoring
        # the multiplying prefix). PubChem CID 120420.
        ("CCOP(=O)(OCC)OC", "diethyl methyl phosphate"),
        # Same shape, roles reversed (two methyls, one ethyl) -- exercises
        # the same grouping/ordering mechanism as CID 120420 above with a
        # different element assignment.
        ("CCOP(=O)(OC)OC", "ethyl dimethyl phosphate"),
    ],
)
def test_mixed_alkyl_phosphate_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Real PubChem structures (P-67.1.3.2), partial ("hydrogen")
        # esters -- one R group + two remaining P-OH ("dihydrogen") or
        # two R groups + one remaining P-OH ("hydrogen").
        ("COP(=O)(O)O", "methyl dihydrogen phosphate"),  # CID 13130
        ("CCOP(=O)(O)O", "ethyl dihydrogen phosphate"),  # CID 74190
        ("COP(=O)(O)OC", "dimethyl hydrogen phosphate"),  # CID 13134
        ("CCOP(=O)(O)OCC", "diethyl hydrogen phosphate"),  # CID 654
    ],
)
def test_partial_hydrogen_phosphate_ester_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_salt_of_partial_ester():
    assert smiles_to_iupac("COP(=O)(O)[O-].[Na+]") == "sodium methyl hydrogen phosphate"
    assert smiles_to_iupac("CCOP(=O)(O)[O-].[Na+]") == "sodium ethyl hydrogen phosphate"
    assert smiles_to_iupac("COP(=O)(O)[O-].[K+]") == "potassium methyl hydrogen phosphate"


def test_salt_of_partial_ester_multivalent_cation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COP(=O)(O)[O-].[Ca+2]")


def test_phosphoric_acid_itself_unaffected():
    # All three P-O positions are plain hydroxyl (no R group at all) --
    # the parent acid itself, not an ester, still out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OP(=O)(O)O")


def test_phosphonic_acid_unaffected():
    assert smiles_to_iupac("CP(=O)(O)O") == "methylphosphonic acid"


def test_phosphanone_unaffected():
    assert smiles_to_iupac("CCP(=O)(CC)CC") == "triethyl-λ5-phosphanone"
