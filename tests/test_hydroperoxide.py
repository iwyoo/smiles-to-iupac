import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-56.1's own worked example: 'ethaneperoxol (PIN)'.
        ("CCOO", "ethaneperoxol"),
        # A mononuclear parent: P-14.3.4.2(a) locant always omitted.
        # Structure verified against PubChem: CID 18199 ("hydroperoxymethane").
        ("COO", "methaneperoxol"),
        # Chain length 3+: 'peroxol' begins with a consonant, so no 'e'
        # elision even with a locant (mirrors _thiol.py's
        # 'propane-2-thiol', not _alcohol.py's vowel-elided 'propan-2-ol').
        # Structure verified against PubChem: CID 123234
        # ("1-hydroperoxypropane").
        ("CCCOO", "propane-1-peroxol"),
        # A branched chain: PubChem CID 18200 ("2-hydroperoxypropane").
        ("CC(C)OO", "propane-2-peroxol"),
        # A halogen substituent coexists freely (P-35.2.1); the -OOH
        # suffix locant is minimized ahead of the substituent-prefix
        # locant (P-44.4.1), so it's cited as C1 even though PubChem's own
        # auto-generated name (CID 55301917, "1-chloro-2-hydroperoxyethane")
        # numbers from the other end.
        ("ClCCOO", "2-chloroethane-1-peroxol"),
    ],
)
def test_hydroperoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCOO")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOC1CCCCC1")


def test_two_hydroperoxide_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OOCCOO")


def test_phenyl_chain_hydroperoxide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#287's
    # carboxylic-acid/.../sulfinamide chains): the ring is cited as a
    # "phenyl" substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCCOO") == "3-phenylpropane-1-peroxol"


def test_phenyl_chain_hydroperoxide_internal_locant():
    assert smiles_to_iupac("C(c1ccccc1)C(C)OO") == "1-phenylpropane-2-peroxol"


def test_hydroperoxybenzene():
    # PubChem IUPACName confirmed directly: "hydroperoxybenzene" (no
    # retained ring-plus-suffix name the way 'phenol' has for -OH; benzene
    # stays the parent with 'hydroperoxy' as a plain prefix instead).
    assert smiles_to_iupac("c1ccccc1OO") == "hydroperoxybenzene"


def test_phenyl_substituted_benzene_ring_hydroperoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CCOO")


def test_phenyl_chain_hydroperoxide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCOO")
