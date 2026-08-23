import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). Note the 'al'
        # locant is never cited (see module docstring): the -CHO carbon is
        # always the chain terminus, so 'propan-1-al' is never written.
        ("CC=O", "ethanal"),
        ("CCC=O", "propanal"),
        ("CC(C)C=O", "2-methylpropanal"),
        # Dialdehyde: multiplying prefix, still no locants at all (both ends
        # of the chain, unambiguous), cross-checked against PubChem.
        ("O=CCCCC=O", "pentanedial"),
        # -al combined with existing unsaturation support: the 'ene' locant
        # is still cited (counted from the -CHO end, suffix priority,
        # P-44.4.1.8), cross-checked against PubChem.
        ("C=CCCC=O", "pent-4-enal"),
        # -al combined with a halogen substituent prefix, cross-checked
        # against PubChem.
        ("ClCCC=O", "3-chloropropanal"),
    ],
)
def test_aldehyde_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formaldehyde_raises():
    # A carbonyl carbon with zero carbon neighbors (formaldehyde) is out of
    # scope for this module (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=O")


def test_ketone_aldehyde_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)CC=O")


def test_carboxylic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)O")


def test_aryl_aldehyde_raises():
    # Benzaldehyde: an aromatic ring elsewhere in the molecule is out of
    # scope for this module (separate, in-progress aromatic-ring module's
    # territory).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=Cc1ccccc1")


def test_alcohol_aldehyde_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC=O")


def test_ring_aldehyde_raises():
    # -CHO on a ring is the 'carbaldehyde' suffix (P-33.3.1.2), a different
    # naming pattern this module deliberately excludes (see module
    # docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC1CCCCC1")


def test_bicyclic_carbon_skeleton_with_stray_aldehyde_raises():
    # A -CHO group whose carbon is not on any single longest chain of the
    # molecule (here, a branch off a longer chain) is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC(C=O)CCCC")
