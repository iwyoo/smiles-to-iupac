import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName).
        ("CN", "methanamine"),
        ("CCN", "ethanamine"),
        # P-14.3.4.2 style locant choice for a real positional choice:
        # 'propan-1-amine' vs 'propan-2-amine', both cross-checked against
        # PubChem.
        ("CCCN", "propan-1-amine"),
        ("CC(N)C", "propan-2-amine"),
        # Diamine: multiplying prefix + full locant set, cross-checked
        # against PubChem.
        ("NCCN", "ethane-1,2-diamine"),
        # Monocyclic saturated ring, single -NH2: P-14.3.3 locant omission
        # (like 'methylcyclohexane') applies to the suffix too. Cross-checked
        # against PubChem.
        ("NC1CCCCC1", "cyclohexanamine"),
        # -NH2 combined with existing unsaturation support, on carbons that
        # don't touch the double bond (avoiding the enamine guard).
        # Cross-checked against PubChem: the -NH2 gets locant 1 (suffix
        # priority, P-44.4.1.8) even though numbering from the other end
        # would give the double bond locant 1 instead.
        ("NCCCC=C", "pent-4-en-1-amine"),
        # -NH2 combined with a halogen substituent prefix, on a chain long
        # enough (3 carbons) that the P-14.3.4.2(b) short-chain omission
        # never applies. Cross-checked against PubChem.
        ("NCCCCl", "3-chloropropan-1-amine"),
        # -NH2 + halogen together on a ring: suffix locant priority
        # (P-44.4.1.8) fixes C1 at the -NH2 carbon, then the halogen gets the
        # lowest remaining locant. Cross-checked against PubChem.
        ("NC1CCCCC1Cl", "2-chlorocyclohexan-1-amine"),
    ],
)
def test_amine_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_secondary_amine_raises():
    # A nitrogen bonded to two carbons (secondary amine) is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNCC")


def test_tertiary_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCN(CC)CC")


def test_nitrile_routes_to_nitrile_module():
    # A nitrogen triple-bonded to carbon is not a plain primary amine; it is
    # routed to the dedicated nitrile module instead (see test_nitrile.py).
    assert smiles_to_iupac("CCC#N") == "propanenitrile"


def test_aniline_raises():
    # An aromatic ring bearing -NH2 is out of scope for this module
    # (separate, in-progress aromatic-ring module's territory).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Nc1ccccc1")


def test_bicyclic_amine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CC2CCC1CC2")


def test_enamine_raises():
    # -NH2 on a carbon that is also part of a C=C bond: deliberately
    # narrowed out of scope (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC=CC")
