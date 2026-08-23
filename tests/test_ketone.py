import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). Note the locant is
        # cited even though C2 is the only chemically possible position
        # (see module docstring): P-14.3.4.2's omission only covers
        # mononuclear/two-carbon parents, neither of which a ketone can be.
        ("CC(=O)C", "propan-2-one"),
        ("CCC(=O)C", "butan-2-one"),
        # P-14.3.4.2 style locant choice for a real positional choice on a
        # longer chain, cross-checked against PubChem.
        ("CCCC(=O)CC", "hexan-3-one"),
        # Diketone: multiplying prefix + full locant set, cross-checked
        # against PubChem.
        ("CC(=O)CC(=O)C", "pentane-2,4-dione"),
        # Monocyclic saturated ring, single ketone: P-14.3.3 locant omission
        # (like 'methylcyclohexane') applies to the suffix too. Cross-checked
        # against PubChem.
        ("O=C1CCCCC1", "cyclohexanone"),
        # -one combined with existing unsaturation support, on carbons that
        # don't touch the double bond (avoiding the enone guard).
        # Cross-checked against PubChem: the ketone gets locant 2 (suffix
        # priority, P-44.4.1.8) even though numbering from the other end
        # would give the double bond a lower locant instead.
        ("CC(=O)CCC=C", "hex-5-en-2-one"),
        # -one combined with a halogen substituent prefix. Cross-checked
        # against PubChem.
        ("CC(=O)CCCl", "4-chlorobutan-2-one"),
        # -one + halogen together on a ring: suffix locant priority
        # (P-44.4.1.8) fixes C1 at the ketone carbon, then the halogen gets
        # the lowest remaining locant. Cross-checked against PubChem.
        ("O=C1CCCCC1Cl", "2-chlorocyclohexan-1-one"),
        # An alpha,beta-unsaturated ketone (the C=C adjacent to, but not on,
        # the carbonyl carbon): unlike an enol/enamine, a ketone carbon can
        # never itself also be an alkene carbon (see module docstring), so
        # this is a perfectly ordinary combined suffix name. Cross-checked
        # against PubChem.
        ("CC(=O)C=CC", "pent-3-en-2-one"),
    ],
)
def test_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aldehyde_raises():
    # A carbonyl carbon with only one carbon neighbor is an aldehyde, a more
    # senior characteristic group (Table 3.3) than a plain ketone; must not
    # be misread as one.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=O")


def test_carboxylic_acid_not_misread_as_ketone():
    # A carboxylic acid is routed to the dedicated carboxylic-acid module
    # (see test_carboxylic_acid.py) instead of falling through here.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_aryl_ketone_raises():
    # A carbonyl on an aromatic ring (an aryl ketone) is out of scope for
    # this module (separate, in-progress aromatic-ring module's territory).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)c1ccccc1")


def test_bicyclic_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2CCC1CC2")


def test_ether_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOCC")


def test_alcohol_hetero_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC(=O)C")
