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


def test_guanidine_not_confused_with_amidine():
    # H2N-C(=NH)-NH2 has two amino nitrogens on the same carbon -- not a
    # plain primary amidine (P-66.4.1.2.1), but its own separate retained
    # name, "guanidine" (see `_guanidine.py`/test_guanidine.py).
    assert smiles_to_iupac("NC(=N)N") == "guanidine"


def test_imine_not_misnamed_as_amidine():
    # A plain imine (`_imine.py`) has only one nitrogen and must not be
    # routed here.
    assert smiles_to_iupac("CC=N") == "ethanimine"


def test_amide_not_misnamed_as_amidine():
    # A plain amide (`_amide.py`) has an oxygen, not a second nitrogen,
    # and must not be routed here.
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"


def test_amidine_unspecified_stereocenter_unaffected():
    # A genuine chain stereocenter left unspecified (no @/@@) is named
    # exactly as before -- no error, matching this project's long-standing
    # convention for unspecified stereochemistry.
    assert smiles_to_iupac("CCC(C)C(=N)N") == "2-methylbutanimidamide"


def test_amidine_specified_chain_stereocenter_raises():
    # This amidine's own C=NH imine bond is always an unspecified
    # potential Bond_Double stereo element to RDKit, regardless of
    # substituents (module docstring) -- so a specified chain
    # stereocenter here always coexists with it, and
    # `specified_stereocenters` correctly rejects the combination (P-92/
    # P-93) instead of the silent drop this project's stereodescriptor
    # safety net exists to fix.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(=N)N")


def test_phenyl_chain_amidine():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_amide.py`'s
    # `test_phenyl_chain_amide`): the ring is cited as a "phenyl"
    # substituent prefix. PubChem PUG REST: "2-phenylethanimidamide"/
    # "3-phenylpropanimidamide".
    assert smiles_to_iupac("c1ccccc1CC(=N)N") == "2-phenylethanimidamide"
    assert smiles_to_iupac("c1ccccc1CCC(=N)N") == "3-phenylpropanimidamide"


def test_phenyl_directly_attached_amidine_raises():
    # Benzamidine-type naming (the amidine carbon directly on the ring) is
    # a separate construction, out of scope for this acyclic-chain-parent
    # module, mirroring `_amide.py`'s benzamide-type rejection.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=N)N")


def test_phenyl_substituted_benzene_ring_amidine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=N)N")


def test_phenyl_chain_amidine_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propanimidamide" for
    # this structure.
    assert smiles_to_iupac("Clc1ccc(CCC(=N)N)cc1") == "3-(4-chlorophenyl)propanimidamide"


def test_phenyl_chain_amidine_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=N)N")
