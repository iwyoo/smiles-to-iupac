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
