import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName).
        ("CO", "methanol"),
        ("CCO", "ethanol"),
        # P-14.3.4.2 style locant choice for a real positional choice:
        # 'propan-1-ol' vs 'propan-2-ol', both cross-checked against PubChem.
        ("CCCO", "propan-1-ol"),
        ("CC(O)C", "propan-2-ol"),
        # Diol / triol: multiplying prefix + full locant set, both
        # cross-checked against PubChem (glycol/glycerol's systematic names).
        ("OCCO", "ethane-1,2-diol"),
        ("OCC(O)CO", "propane-1,2,3-triol"),
        # Monocyclic saturated ring, single -OH: P-14.3.3 locant omission
        # (like 'methylcyclohexane') applies to the suffix too. Cross-checked
        # against PubChem.
        ("OC1CCCCC1", "cyclohexanol"),
        # -OH combined with existing unsaturation support, on carbons that
        # don't touch the double bond (avoiding the enol guard). Cross-
        # checked against PubChem: the -OH gets locant 1 (suffix priority,
        # P-44.4.1.8) even though numbering from the other end would give
        # the double bond locant 1 instead.
        ("OCCCC=C", "pent-4-en-1-ol"),
        # -OH combined with a halogen substituent prefix, on a chain long
        # enough (3 carbons) that the P-14.3.4.2(b) short-chain omission
        # never applies, sidestepping the ambiguous two-carbon case (see
        # module docstring). Cross-checked against PubChem.
        ("OCCCCl", "3-chloropropan-1-ol"),
        ("OCC(F)CC", "2-fluorobutan-1-ol"),
        # -OH + halogen together on a ring: suffix locant priority
        # (P-44.4.1.8) fixes C1 at the -OH carbon, then the halogen gets the
        # lowest remaining locant. Cross-checked against PubChem.
        ("OC1CCCCC1Cl", "2-chlorocyclohexan-1-ol"),
    ],
)
def test_alcohol_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_carboxylic_acid_raises():
    # -C(=O)-OH: a carboxylic acid is a more senior characteristic group
    # (Table 3.3) than a plain alcohol; must not be misread as one.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)O")


def test_ketone_raises():
    # C=O not part of -COOH: a ketone is a more senior characteristic group
    # (Table 3.3) than a plain alcohol.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)C")


def test_aldehyde_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC=O")


def test_phenol_raises():
    # An aromatic ring bearing -OH is out of scope for this module (separate,
    # in-progress aromatic-ring module's territory); must not be picked up
    # here before reaching that module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Oc1ccccc1")


def test_bicyclic_alcohol_raises():
    # Same skeleton as test_halogens.py's 'ClC1CC2CCC1CC2' ->
    # '2-chlorobicyclo[2.2.2]octane', with the halogen swapped for -OH: a
    # real, nameable compound ('bicyclo[2.2.2]octan-2-ol', confirmed via
    # PubChem), but von Baeyer polycyclic alcohols are explicitly deferred.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CC2CCC1CC2")


def test_ether_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOCC")


def test_enol_raises():
    # -OH on a carbon that is also part of a C=C bond: deliberately
    # narrowed out of scope (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC")


def test_amine_hetero_mix_raises():
    # A structure with both -OH and a non-halogen heteroatom (N) is rejected
    # outright by the alcohol module; this module never attempts
    # suffix-vs-suffix seniority competition (Table 3.3) since only
    # C/O(-OH)/halogen atoms are accepted at all. (A pure amine, no -OH, is
    # dispatched to the separate amine module instead — see test_amine.py.)
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCCN")
