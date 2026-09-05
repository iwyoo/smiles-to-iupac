import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_acetylpiperidine():
    # Blue Book P-64.1.2.1(b) worked example directly: "1-(piperidin-1-yl)
    # ethan-1-one (PIN)" / "1-acetylpiperidine (a 'hidden amide')".
    assert smiles_to_iupac("CC(=O)N1CCCCC1") == "1-(piperidin-1-yl)ethan-1-one"


def test_acetylpyrrolidine():
    # PubChem CID 77650 "1-pyrrolidin-1-ylethanone".
    assert smiles_to_iupac("CC(=O)N1CCCC1") == "1-(pyrrolidin-1-yl)ethan-1-one"


def test_propanoylpiperidine():
    # Blue Book P-64.3.2 worked example: "1-(piperidin-1-yl)propan-1-one
    # (PIN)" / "1-propanoylpiperidine".
    assert smiles_to_iupac("CCC(=O)N1CCCCC1") == "1-(piperidin-1-yl)propan-1-one"


def test_propanoylpyrrolidine():
    # PubChem CID 551043 "1-pyrrolidin-1-ylpropan-1-one".
    assert smiles_to_iupac("CCC(=O)N1CCCC1") == "1-(pyrrolidin-1-yl)propan-1-one"


def test_acetylazepane():
    # PubChem CID 22059 "1-(azepan-1-yl)ethanone".
    assert smiles_to_iupac("CC(=O)N1CCCCCC1") == "1-(azepan-1-yl)ethan-1-one"


def test_chloroacetylpiperidine():
    # PubChem CID 222312 "2-chloro-1-piperidin-1-ylethanone" -- this
    # project's usual explicit-locant/parenthesization convention applied.
    assert smiles_to_iupac("ClCC(=O)N1CCCCC1") == "2-chloro-1-(piperidin-1-yl)ethan-1-one"


def test_formylpiperidine_raises():
    # R = H (a formyl group) is a distinct 'carbaldehyde' suffix
    # construction (PubChem CID 17429 "piperidine-1-carbaldehyde"), not a
    # pseudoketone -- out of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CN1CCCCC1")


def test_branched_acyl_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)N1CCCCC1")


def test_substituted_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N1CCC(C)CC1")


def test_two_heteroatom_ring_raises():
    # Morpholine (1,4-N,O) attaching via N is a separate axis (locant 4,
    # not 1) -- out of scope for this first slice.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N1CCOCC1")


def test_plain_ring_lactam_unaffected():
    # A genuine lactam (the ketone carbon itself in the ring) is a
    # different shape entirely and must still route correctly.
    assert smiles_to_iupac("O=C1CCCCN1") == "piperidin-2-one"
