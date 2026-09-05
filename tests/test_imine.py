import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName).
        ("C=N", "methanimine"),
        ("CC=N", "ethanimine"),
        # A ketimine (interior C=N) always cites its locant, even on the
        # shortest possible chain -- like 'propan-2-one', 'propan-2-imine'
        # (acetone imine) still cites '2'.
        ("CC(C)=N", "propan-2-imine"),
        ("CCC(C)=N", "butan-2-imine"),
        # An aldimine (terminal C=N) on a 3+-carbon chain *does* cite its
        # locant, unlike '-al' -- the Blue Book's own P-62.3.1.1 worked
        # example is 'hexan-1-imine (PIN)'.
        ("CCC=N", "propan-1-imine"),
        ("CCCCCC=N", "hexan-1-imine"),
        # A two-carbon aldimine chain omits the locant (P-14.3.4.2(b)) --
        # and, unlike this project's own `_alcohol.py`/`_ketone.py`, that
        # omission isn't gated on having zero other substituents (see
        # module docstring).
        ("ClCC=N", "2-chloroethanimine"),
        # An N-substituent is cited as an "N-" prefix with no locant.
        ("C=NC", "N-methylmethanimine"),
        ("CC=NC", "N-methylethanimine"),
        ("CC(C)=NC", "N-methylpropan-2-imine"),
        # Oximes (P-68.3.1.1.2): the PIN is the N-hydroxy derivative of the
        # imine named by this module -- the Blue Book's own worked example
        # is 'N-hydroxypentan-2-imine (PIN)' for pentan-2-one oxime. Note
        # this deliberately does NOT match PubChem's own auto-generated
        # name for this exact structure (CID 136433,
        # "N-pentan-2-ylidenehydroxylamine", a different hydroxylamine-
        # parent pattern) -- implemented per the Blue Book's direct PIN
        # citation instead (see module docstring).
        ("CC(=NO)CCC", "N-hydroxypentan-2-imine"),
        # An O-alkyl oxime ether: PubChem CID 54150571 matches the Blue
        # Book's own 'N-ethoxypropan-1-imine (PIN)' worked example exactly
        # (unlike the plain -OH case above).
        ("CCC=NOCC", "N-ethoxypropan-1-imine"),
        ("C=NO", "N-hydroxymethanimine"),
        ("CC=NO", "N-hydroxyethanimine"),
    ],
)
def test_smiles_to_iupac_simple_imine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_raises():
    # A cyclic/aromatic imine (e.g. 'thiolan-2-imine') is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N=C1CCCC1")


def test_multiple_imine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N=CC=N")


def test_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NC(C)C")


def test_unsaturated_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC=N")


def test_amine_hetero_mix_raises():
    # A structure with both an imine and a coexisting -OH is rejected
    # outright -- this module doesn't attempt Table 3.3 seniority
    # competition between suffixes.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC=N")


def test_branched_oxime_o_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NOC(C)C")


def test_two_oxygens_on_oxime_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NOO")


def test_imine_unspecified_stereocenter_unaffected():
    # A genuine chain stereocenter left unspecified (no @/@@) is named
    # exactly as before -- no error, matching this project's long-standing
    # convention for unspecified stereochemistry.
    assert smiles_to_iupac("CCC(C)C(C)=N") == "3-methylpentan-2-imine"


def test_imine_specified_chain_stereocenter_raises():
    # This module's own C=N bond is always an unspecified potential
    # Bond_Double stereo element to RDKit, regardless of substituents
    # (module docstring) -- so a specified chain stereocenter here always
    # coexists with it, and `specified_stereocenters` correctly rejects
    # the combination (P-92/P-93) instead of the silent drop this
    # project's stereodescriptor safety net exists to fix.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)C(C)=N")


def test_phenyl_chain_imine():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring
    # `_sulfonic_acid.py`'s phenyl-chain path): the ring is cited as a
    # "phenyl" substituent prefix. PubChem PUG REST: "2-phenylethanimine"/
    # "3-phenylpropan-1-imine"/"1-phenylpropan-2-imine" (ketimine, internal
    # locant) -- all exact matches, including the two-carbon and
    # mononuclear locant-omission rules this module already gets right.
    assert smiles_to_iupac("c1ccccc1CC=N") == "2-phenylethanimine"
    assert smiles_to_iupac("c1ccccc1CCC=N") == "3-phenylpropan-1-imine"
    assert smiles_to_iupac("c1ccccc1CC(C)=N") == "1-phenylpropan-2-imine"


def test_phenyl_directly_attached_imine():
    # The imine carbon directly on the ring (chain length 1, e.g.
    # benzaldimine) is still a valid mononuclear-chain case for this
    # module -- P-14.3.4.2(a) omits the phenyl substituent's own locant
    # too. PubChem PUG REST: "phenylmethanimine".
    assert smiles_to_iupac("c1ccccc1C=N") == "phenylmethanimine"


def test_phenyl_substituted_benzene_ring_imine_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-5.md).
    assert smiles_to_iupac("Cc1ccccc1CC=N") == "2-(2-methylphenyl)ethanimine"


def test_phenyl_chain_imine_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propan-1-imine" for
    # this structure.
    assert smiles_to_iupac("Clc1ccc(CCC=N)cc1") == "3-(4-chlorophenyl)propan-1-imine"


def test_phenyl_chain_imine_ring_methyl():
    assert smiles_to_iupac("Cc1ccc(CCC=N)cc1") == "3-(4-methylphenyl)propan-1-imine"


def test_phenyl_chain_imine_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC=N")


def test_phenyl_chain_n_substituted_imine():
    # An N-alkyl-substituted imine alongside a benzene-ring-substituent
    # chain. PubChem PUG REST: "N-methyl-3-phenylpropan-1-imine"/
    # "N-ethyl-2-phenylethanimine" -- both exact matches, including the
    # phenyl-chain locant rules already established for the plain case.
    assert smiles_to_iupac("c1ccccc1CCC=NC") == "N-methyl-3-phenylpropan-1-imine"
    assert smiles_to_iupac("c1ccccc1CC=NCC") == "N-ethyl-2-phenylethanimine"


def test_phenyl_chain_oxime():
    # A plain -OH oxime alongside a benzene-ring-substituent chain. No
    # PubChem-registered example exists for this combination (PubChem's
    # own auto-generated oxime name uses a different "N-...ylidene-
    # hydroxylamine" pattern this project already diverges from in the
    # non-benzene case, per this module's own docstring) -- this follows
    # the same "N-hydroxy..." pattern this module already established for
    # the plain-chain oxime case, combined with the phenyl-chain locant
    # rules already established for the plain imine case.
    assert smiles_to_iupac("c1ccccc1CCC=NO") == "N-hydroxy-3-phenylpropan-1-imine"


def test_phenyl_chain_oxime_ether_raises():
    # An O-alkyl oxime ether (=N-O-R) alongside a benzene-ring-substituent
    # chain remains out of scope -- no PubChem-registered example exists
    # to verify the combination against.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCC=NOCC")
