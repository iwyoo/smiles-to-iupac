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
        # Monocyclic ring, single -NH2, single ring double bond (P-31.1.3):
        # the amine always gets locant 1 (suffix priority), the ring
        # double bond's locant is minimized by choosing direction.
        # Cross-checked against PubChem (CID 10931371/13040867).
        ("NC1CCCC=C1", "cyclohex-2-en-1-amine"),
        ("NC1CC=CCC1", "cyclohex-3-en-1-amine"),
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


def test_n_locant_sorts_before_numeric_when_names_coincide():
    # PubChem PUG REST IUPACName match: "N,2-dimethylpropan-1-amine", not
    # "2,N-dimethyl..." -- a non-numeric ('N') locant is cited before a
    # numeric one when both share the same substituent name
    # (`_substituents.py`'s `_locant_sort_key`).
    assert smiles_to_iupac("CNCC(C)C") == "N,2-dimethylpropan-1-amine"


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


def test_secondary_amine_with_halogen_on_parent_chain():
    # P-14.5.2: the 'N-' prefix interleaves alphabetically with the
    # halogen prefix rather than always citing first (PubChem:
    # "2-chloro-N-ethylethanamine"; this project keeps the amine's own
    # locant, same established style as "2-phenylethan-1-amine" above).
    assert smiles_to_iupac("ClCCNCC") == "2-chloro-N-ethylethan-1-amine"


def test_secondary_amine_with_dihalogen_on_parent_chain():
    assert smiles_to_iupac("ClCCC(Cl)NCC") == "1,3-dichloro-N-ethylpropan-1-amine"


def test_nitrile_routes_to_nitrile_module():
    # A nitrogen triple-bonded to carbon is not a plain primary amine; it is
    # routed to the dedicated nitrile module instead (see test_nitrile.py).
    assert smiles_to_iupac("CCC#N") == "propanenitrile"


def test_aniline():
    # -NH2 attached directly to a benzene ring carbon (P-62.2.1.1.1). The
    # retained name 'aniline' stands for the whole ring+NH2 system (like
    # 'phenol'), so the amine's own ring locant is never cited, only other
    # substituents'. PubChem PUG REST: "aniline"/"4-methylaniline"/
    # "2-chloroaniline" -- all exact matches.
    assert smiles_to_iupac("Nc1ccccc1") == "aniline"
    assert smiles_to_iupac("Nc1ccc(C)cc1") == "4-methylaniline"
    assert smiles_to_iupac("Nc1ccccc1Cl") == "2-chloroaniline"


def test_n_substituted_aniline_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNc1ccccc1")


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


def test_phenyl_chain_amine():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_thiol.py`'s
    # `test_phenyl_chain_thiol`): the ring is cited as a "phenyl"
    # substituent prefix. PubChem PUG REST: "phenylmethanamine"/
    # "3-phenylpropan-1-amine" (both exact matches; the two-carbon case
    # "2-phenylethan-1-amine" keeps its locant unlike PubChem's
    # locant-omitted "2-phenylethanamine" -- same known, out-of-scope-here
    # limitation as `_thiol.py`'s own two-carbon phenyl-chain case).
    assert smiles_to_iupac("c1ccccc1CN") == "phenylmethanamine"
    assert smiles_to_iupac("c1ccccc1CCN") == "2-phenylethan-1-amine"
    assert smiles_to_iupac("c1ccccc1CCCN") == "3-phenylpropan-1-amine"


def test_phenyl_chain_amine_internal_locant():
    # The -NH2 locant is a genuine choice on the chain, same as the base
    # acyclic module. PubChem PUG REST: "1-phenylpropan-2-amine".
    assert smiles_to_iupac("C(c1ccccc1)C(C)N") == "1-phenylpropan-2-amine"


def test_phenyl_directly_attached_amine():
    # Same aniline-type case as test_aniline, reached via a different
    # ring-atom ordering in the SMILES.
    assert smiles_to_iupac("c1ccccc1N") == "aniline"


def test_phenyl_substituted_benzene_ring_amine_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-4.md) -- PubChem PUG
    # REST-verified "2-(2-methylphenyl)ethanamine".
    assert smiles_to_iupac("Cc1ccccc1CCN") == "2-(2-methylphenyl)ethan-1-amine"


def test_phenyl_chain_amine_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propan-1-amine" for
    # this structure.
    assert smiles_to_iupac("Clc1ccc(CCCN)cc1") == "3-(4-chlorophenyl)propan-1-amine"


def test_phenyl_chain_amine_ring_methyl():
    # PubChem PUG REST-verified "3-(4-methylphenyl)propan-1-amine".
    assert smiles_to_iupac("Cc1ccc(CCCN)cc1") == "3-(4-methylphenyl)propan-1-amine"


def test_phenyl_chain_amine_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCN")


def test_phenyl_chain_diamine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(N)CCN")


def test_phenyl_chain_secondary_amine_raises():
    # A secondary/tertiary amine nitrogen alongside a ring is out of
    # scope, same as the existing (non-phenyl-chain) ring guard in
    # `name_amine`.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCNC")


def test_unsaturated_ring_amine_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC=C1C")


def test_unsaturated_ring_amine_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC#C1")
