import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName) and matching the
        # Blue Book's own worked examples (P6a.pdf, P-66.4.1.1). The
        # amidine carbon is always C1, and its own locant is never cited
        # (P-14.3.3), same as `_amide.py`.
        ("C(=N)N", "methanimidamide"),
        ("CC(=N)N", "ethanimidamide"),
        ("CCC(=N)N", "propanimidamide"),
        ("CCCCCC(=N)N", "hexanimidamide"),
        # A real positional choice for a substituent: the amidine carbon
        # fixes C1 regardless.
        ("CC(C)C(=N)N", "2-methylpropanimidamide"),
        # -imidamide combined with existing unsaturation support.
        ("C=CCC(=N)N", "but-3-enimidamide"),
        # -imidamide + halogen substituent prefix.
        ("ClCC(=N)N", "2-chloroethanimidamide"),
    ],
)
def test_amidine_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_amino_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=N)NC")


def test_n_substituted_imino_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=NC)N")


def test_ring_amidine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C(=N)N")


def test_diamidine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=N)CC(=N)N")


def test_guanidine_raises():
    # H2N-C(=NH)-NH2 has two amino nitrogens on the same carbon -- not a
    # plain primary amidine (P-66.4.1.2.1, a separate retained name), out
    # of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=N)N")


def test_imine_not_misnamed_as_amidine():
    # A plain imine (`_imine.py`) has only one nitrogen and must not be
    # routed here.
    assert smiles_to_iupac("CC=N") == "ethanimine"


def test_amide_not_misnamed_as_amidine():
    # A plain amide (`_amide.py`) has an oxygen, not a second nitrogen,
    # and must not be routed here.
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"
