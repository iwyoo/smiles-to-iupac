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
        # P-16.3.3 multiplying-prefix elision: 'tetra' elides its terminal
        # 'a' before the vowel-initial 'amine' suffix ('tetramine', not
        # 'tetraamine'). Cross-checked against PubChem PUG REST IUPACName,
        # CID 6395580.
        ("NCC(N)C(N)CN", "butane-1,2,3,4-tetramine"),
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


def test_secondary_amine():
    # PubChem structure match: diethylamine -> "N-ethylethanamine".
    assert smiles_to_iupac("CCNCC") == "N-ethylethanamine"


def test_tertiary_amine_symmetric():
    # PubChem structure match: triethylamine -> "N,N-diethylethanamine".
    assert smiles_to_iupac("CCN(CC)CC") == "N,N-diethylethanamine"


def test_tertiary_amine_all_same():
    # PubChem structure match: trimethylamine -> "N,N-dimethylmethanamine".
    assert smiles_to_iupac("CN(C)C") == "N,N-dimethylmethanamine"


def test_tertiary_amine_asymmetric():
    # PubChem structure match: "N-ethyl-N-methylpropan-1-amine" -- the
    # longest N-linked chain (propyl) becomes the parent, the other two
    # (ethyl, methyl) are cited as alphabetized N-prefixes.
    assert smiles_to_iupac("CCN(C)CCC") == "N-ethyl-N-methylpropan-1-amine"


def test_secondary_amine_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCNC(C)C")


def test_secondary_amine_unsaturated_n_substituent_raises():
    # The allyl group is the smaller N-linked chain here (pentyl is longer
    # and becomes the parent), so it's cited as the N-substituent -- and an
    # unsaturated N-substituent is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCNCCCCC")


def test_secondary_amine_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(C)C1CCCCC1")


def test_diamine_with_tertiary_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCCN(C)C")


def test_secondary_amine_with_halogen_on_n_substituent_raises():
    # A halogen on the smaller (N-substituent) chain would otherwise be
    # silently dropped, since only that chain's length is used to name it.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClCCNCCC")


def test_secondary_amine_with_halogen_on_parent_chain_raises():
    # P-14.5.2's alphanumerical interleaving of the 'N-' prefix with other
    # substituent prefixes isn't implemented yet (PubChem:
    # "2-chloro-N-ethylethanamine") -- reject rather than silently mis-order
    # ("N-ethyl-2-chloro...").
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClCCNCC")


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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92) on a primary
        # amine chain -- the amine nitrogen itself is never a potential
        # stereocenter (pyramidal inversion), confirmed via RDKit
        # `FindPotentialStereo`. PubChem CID 2724537.
        ("CC[C@@H](C)N", "(2R)-butan-2-amine"),
        ("CC[C@H](C)N", "(2S)-butan-2-amine"),
    ],
)
def test_primary_amine_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_secondary_tertiary_amine_stereocenter():
    # The stereodescriptor sits outermost, ahead of both the N,N- and N-
    # prefixes (P-91.3: always at the very front of the complete name),
    # same pattern as `_amide.py`.
    assert smiles_to_iupac("CC[C@@H](C)N(C)C") == "(2R)-N,N-dimethylbutan-2-amine"
    assert smiles_to_iupac("CC[C@@H](C)NCC") == "(2R)-N-ethylbutan-2-amine"


def test_cyclic_amine_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`. PubChem CID
    # 92244014 (name carries a redundant 'trans-' relative descriptor this
    # project drops once full R/S is given, same policy as the existing
    # ring-alcohol stereocenter task).
    assert smiles_to_iupac("N[C@H]1CCCC[C@@H]1Cl") == "(1S,2S)-2-chlorocyclohexan-1-amine"


def test_amine_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(Cl)N") == "1-chloropropan-1-amine"


def test_amine_partially_specified_stereocenters_raises():
    # Only one of the ring's two genuine stereocenters is marked -- must
    # raise rather than silently dropping the marker.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@H]1CCCCC1Cl")
