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


def test_secondary_amine_branched_n_substituent():
    # PubChem CID 89119: 'N-propan-2-ylpropan-1-amine' -- corrected to this
    # project's established compound-prefix parenthesization (see
    # `_amide.py`'s identical 'N-(propan-2-yl)acetamide' correction).
    assert smiles_to_iupac("CCCNC(C)C") == "N-(propan-2-yl)propan-1-amine"


def test_tertiary_amine_two_branched_n_substituents():
    # PubChem CID 13616122: 'N,N-di(propan-2-yl)propan-1-amine'.
    assert smiles_to_iupac("CCCN(C(C)C)C(C)C") == "N,N-di(propan-2-yl)propan-1-amine"


def test_secondary_amine_halogenated_n_substituent():
    # PubChem CID 3045065: 'N-2-chloroethylpropan-1-amine' -- corrected to
    # this project's established compound-prefix parenthesization.
    assert smiles_to_iupac("ClCCNCCC") == "N-(2-chloroethyl)propan-1-amine"


def test_secondary_amine_allyl_n_substituent():
    # PubChem CID 12442641 lists 'N-prop-2-enylpentan-1-amine'; the PIN cites
    # the free-valence locant (P-32.1.1: 'prop-2-en-1-yl').
    assert smiles_to_iupac("C=CCNCCCCC") == "N-(prop-2-en-1-yl)pentan-1-amine"


def test_secondary_amine_butenyl_n_substituent():
    # PubChem CID 11018968 lists 'N-but-3-enylpentan-1-amine'; PIN style cites '1'.
    assert smiles_to_iupac("C=CCCNCCCCC") == "N-(but-3-en-1-yl)pentan-1-amine"


def test_secondary_amine_propargyl_n_substituent():
    # PubChem CID 3465926 lists 'N-prop-2-ynylbutan-1-amine'; PIN style cites '1'.
    assert smiles_to_iupac("C#CCNCCCC") == "N-(prop-2-yn-1-yl)butan-1-amine"


def test_secondary_amine_vinyl_n_substituent_raises():
    # N directly on the alkene carbon (not one bond further away, as the
    # allyl/butenyl/propargyl cases above are) is still out of scope --
    # a real enamine-shaped nitrogen, rejected by `_reject_enamine_carbons`
    # regardless of which side of the amine it sits on.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNCCCC")


def test_secondary_amine_branched_unsaturated_n_substituent():
    assert smiles_to_iupac("C=CC(C)NCCCCC") == "N-(but-3-en-2-yl)pentan-1-amine"


def test_secondary_amine_two_multiple_bonds_on_n_substituent():
    assert smiles_to_iupac("C=CC=CCNCCCCCCC") == "N-(penta-2,4-dien-1-yl)heptan-1-amine"


def test_secondary_amine_on_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(C)C1CCCCC1")


def test_diamine_with_tertiary_nitrogen():
    # P-16.9.2: two coexisting amine nitrogens use superscript parent-chain
    # locants ('N1', 'N2', ...) rather than the older N/N'-prime
    # convention to tell them apart -- confirmed against the Blue Book's
    # own worked PIN example for this exact rule, see
    # `test_two_amines_use_superscript_locants_not_primes` below.
    assert smiles_to_iupac("NCCN(C)C") == "N1,N1-dimethylethane-1,2-diamine"


def test_disconnected_amine_nitrogens_named():
    assert smiles_to_iupac("NCCN(C)CCN") == "N1-(2-aminoethyl)-N1-methylethane-1,2-diamine"


def test_triamine_with_tertiary_nitrogen():
    # Three coexisting amines on one shared chain generalizes cleanly from
    # the two-amine case (P-16.9.2's superscript-locant mechanism doesn't
    # special-case a count) -- cross-checked structurally against PubChem
    # CID 14396522 (same skeleton, C5H15N3; PubChem's own autoname uses
    # yet another, non-PIN locant style, '3-N,3-N-dimethyl...', so it's
    # not usable for a name-string comparison here, only formula/skeleton).
    assert smiles_to_iupac("CN(C)CC(N)CN") == "N1,N1-dimethylpropane-1,2,3-triamine"


def test_tetramine_with_tertiary_nitrogen():
    assert smiles_to_iupac("NCC(N)C(N)CN(C)C") == "N1,N1-dimethylbutane-1,2,3,4-tetramine"


def test_two_amines_use_superscript_locants_not_primes():
    # Blue Book P-16.9.2's own worked PIN example, verbatim.
    assert smiles_to_iupac("CCNCCNC") == "N1-ethyl-N2-methylethane-1,2-diamine"


def test_geminal_diamine_with_substituent_raises():
    # P-16.9.2 explicitly excludes geminal amines from the superscript-
    # locant convention; both nitrogens sharing one carbon needs a
    # different (unimplemented) naming path.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(N)N(C)C")


def test_two_amines_with_unsaturated_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCCN(C)CC=C")


def test_secondary_amine_with_halogen_on_parent_chain():
    # P-14.5.2: the 'N-' prefix interleaves alphabetically with the
    # halogen prefix rather than always citing first (PubChem:
    # "2-chloro-N-ethylethanamine"; this project keeps the amine's own
    # locant, same established style as "2-phenylethanamine" above).
    assert smiles_to_iupac("ClCCNCC") == "2-chloro-N-ethylethanamine"


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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CNc1ccccc1", "N-methylaniline"),
        ("CN(C)c1ccccc1", "N,N-dimethylaniline"),
        ("CCNc1ccc(Cl)cc1", "4-chloro-N-ethylaniline"),
        ("c1ccccc1Nc1ccccc1", "N-phenylaniline"),
        ("c1ccccc1N(c1ccccc1)c1ccccc1", "N,N-diphenylaniline"),
        ("c1ccccc1CNc1ccccc1", "N-benzylaniline"),
        ("NC(c1ccccc1)c1ccccc1", "diphenylmethanamine"),
    ],
)
def test_n_substituted_aniline_and_multi_ring_chain_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


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


def test_cyclic_amine_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_alcohol.py`'s
    # `_name_cyclic_alcohol`/`_ketone.py`/`_thiol.py`/`_sulfonic_acid.py` --
    # the branch and the -NH2 share the same ring carbon (C1), same shape
    # as the alcohol/thiol/sulfonic-acid precedent. PubChem has no cached
    # record for this exact structure (CID 0, sparse data gap), but the
    # unstereo parent ('1-ethylcyclohexan-1-amine') matches PubChem
    # exactly, and the descriptor-construction code itself is identical to
    # the already-verified alcohol/ketone/thiol/sulfonic-acid cases.
    assert smiles_to_iupac("NC1(CCCCC1)[C@@H](C)CC") == "1-[(2S)-butan-2-yl]cyclohexan-1-amine"


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


def test_acyclic_amine_specified_ez_double_bond():
    # Standalone specified C=C E/Z stereo, no stereocenter -- mirrors
    # `_ketone.py`/`_aldehyde.py`/`_carboxylic_acid.py`/`_nitrile.py`'s
    # identical PR #710/#713/#715/#717 case. PubChem CID 6433772 confirms
    # the structure (SMILES "C/C=C/CN"), though PubChem's own
    # auto-generated name is "(E)-but-2-en-1-amine" (no locant); this
    # project's existing convention is to always cite the locant
    # explicitly.
    assert smiles_to_iupac("NC/C=C/C") == "(2E)-but-2-en-1-amine"


def test_acyclic_primary_amine_stereocenter_with_ez_double_bond_coexistence():
    # P-91.3: a specified tetrahedral stereocenter and a specified C=C
    # double-bond E/Z element on the same acyclic primary-amine chain are
    # cited together in one ascending-locant group, mirroring
    # `_alcohol.py`'s/`_ketone.py`'s/`_aldehyde.py`'s/
    # `_carboxylic_acid.py`'s/`_nitrile.py`'s identical case. PubChem CIDs
    # 58781915/93518073 confirm both structures (SMILES
    # "C/C=C/[C@@H](C)N" / "C/C=C/[C@H](C)N"), though PubChem's own names
    # group by type ("(E,2R)-...") rather than by ascending locant -- the
    # expected locant-ascending form here is the primary source's own
    # P-91.3 rule.
    assert smiles_to_iupac("C/C=C/[C@@H](C)N") == "(2R,3E)-pent-3-en-2-amine"
    assert smiles_to_iupac("C/C=C/[C@H](C)N") == "(2S,3E)-pent-3-en-2-amine"


def test_acyclic_secondary_amine_stereocenter_with_ez_double_bond_coexistence():
    # Same combined-descriptor mechanism, reached via
    # `_name_acyclic_secondary_tertiary_amine` instead of the simple
    # primary-amine path -- a mechanical side effect of sharing
    # `_best_chain_name`/the stereo format with the primary-amine case,
    # not independently PubChem-verified (see issue #718).
    assert smiles_to_iupac("CN[C@H](C)/C=C/C") == "(2R,3E)-N-methylpent-3-en-2-amine"


def test_ammonium_stereocenter_with_ez_double_bond_coexistence():
    # `_ammonium.py` reuses `_name_acyclic_secondary_tertiary_amine`
    # directly (see issue #718's `_carboxylic_acid_amine.py`-style hidden
    # caller finding) -- confirms its own stereo computation was also
    # updated to the combined-capable helper.
    assert smiles_to_iupac("C[N+](C)(C)[C@H](C)/C=C/C") == "(2R,3E)-N,N,N-trimethylpent-3-en-2-aminium"


def test_phenyl_chain_amine():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_thiol.py`'s
    # `test_phenyl_chain_thiol`): the ring is cited as a "phenyl"
    # substituent prefix. PubChem PUG REST: "phenylmethanamine"/
    # "3-phenylpropan-1-amine" (both exact matches; the two-carbon case
    # "2-phenylethanamine" keeps its locant unlike PubChem's
    # locant-omitted "2-phenylethanamine" -- same known, out-of-scope-here
    # limitation as `_thiol.py`'s own two-carbon phenyl-chain case).
    assert smiles_to_iupac("c1ccccc1CN") == "phenylmethanamine"
    assert smiles_to_iupac("c1ccccc1CCN") == "2-phenylethanamine"
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
    assert smiles_to_iupac("Cc1ccccc1CCN") == "2-(2-methylphenyl)ethanamine"


def test_phenyl_chain_amine_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propan-1-amine" for
    # this structure.
    assert smiles_to_iupac("Clc1ccc(CCCN)cc1") == "3-(4-chlorophenyl)propan-1-amine"


def test_phenyl_chain_amine_ring_methyl():
    # PubChem PUG REST-verified "3-(4-methylphenyl)propan-1-amine".
    assert smiles_to_iupac("Cc1ccc(CCCN)cc1") == "3-(4-methylphenyl)propan-1-amine"


def test_phenyl_chain_amine_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CCN") == "2-(4-ethylphenyl)ethanamine"


def test_phenyl_chain_amine_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCN")


def test_phenyl_chain_diamine_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(N)CCN")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1CCNC", "N-methyl-2-phenylethanamine"),
        ("CN(C)Cc1ccccc1", "N,N-dimethyl-1-phenylmethanamine"),
        ("CC(N(C)C)c1ccccc1", "N,N-dimethyl-1-phenylethanamine"),
        ("CC(N(C)C)c1ccccc1P(c1ccccc1)c1ccccc1", "1-[2-(diphenylphosphanyl)phenyl]-N,N-dimethylethanamine"),
        ("NCCc1ccccc1P(C)C", "2-[2-(dimethylphosphanyl)phenyl]ethanamine"),
    ],
)
def test_phenyl_chain_secondary_tertiary_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_ring_amine_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC=C1C")


def test_unsaturated_ring_amine_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC1CCCC#C1")


def test_ring_substituent_chain_amine():
    # A primary amine entirely on a chain hanging off a plain saturated
    # ring (the ring itself bears no amine) -- the ring is cited as a
    # "cyclo..." substituent prefix on the chain, mirroring `_alcohol.py`'s
    # `_name_ring_substituent_chain_alcohol`. PubChem PUG REST-verified
    # "cyclohexylmethanamine" (CID 76688).
    assert smiles_to_iupac("NCC1CCCCC1") == "cyclohexylmethanamine"


def test_ring_substituent_chain_amine_internal_locant():
    # PubChem PUG REST: "1-cyclohexylethanamine" (CID 110733) -- the
    # locant is not omitted here, the same known limitation as
    # `test_phenyl_chain_amine`'s two-carbon case.
    assert smiles_to_iupac("NC(C)C1CCCCC1") == "1-cyclohexylethanamine"


def test_ring_substituent_chain_amine_ring_with_substituent():
    # A ring atom other than the chain attachment carrying its own
    # substituent is out of scope for this first pass -- mirrors
    # `_name_ring_substituent_chain_alcohol`'s identical restriction.
    assert smiles_to_iupac("NCC1CCC(C)CC1") == "(4-methylcyclohexyl)methanamine"


def test_ring_substituent_chain_amine_unsaturated_ring():
    assert smiles_to_iupac("NCC1CCCC=C1") == "(cyclohex-2-en-1-yl)methanamine"


def test_ring_with_amine_chain_amine_tie():
    # Ring and chain each carry exactly one amine (P-44.1.1 tie,
    # P-44.1.2.2 resolves it in the ring's favor) -- the ring becomes the
    # parent and the chain's amine is cited as an "aminomethyl"
    # substituent prefix, mirroring `_alcohol.py`'s `_name_ring_with_
    # hydroxy_chain_alcohol`. PubChem PUG REST-verified
    # "2-(aminomethyl)cyclohexan-1-amine" (CID 430326).
    assert smiles_to_iupac("NC1CCCCC1CN") == "2-(aminomethyl)cyclohexan-1-amine"


def test_ring_with_amine_chain_amine_ring_wins_outright():
    # The ring carries two amines against the chain's one -- P-44.1.1's
    # greater-count rule picks the ring outright, no tie-break needed.
    # PubChem PUG REST-verified "1-(aminomethyl)cyclohexane-1,2-diamine"
    # (CID 149881934).
    assert (
        smiles_to_iupac("NC1CCCCC1(N)CN") == "1-(aminomethyl)cyclohexane-1,2-diamine"
    )


def test_ring_with_amine_chain_amine_longer_chain():
    # A two-carbon chain amine, same tie-break, verifying the chain's own
    # locant math still works once it's longer than a single carbon.
    assert smiles_to_iupac("NC1CCCCC1CCN") == "2-(2-aminoethyl)cyclohexan-1-amine"


def test_ring_with_amine_chain_amine_chain_wins():
    # The chain has strictly more amines than the ring, so it's the
    # senior parent (P-44.1.1) and the ring is cited as an "amino"-
    # decorated cyclic substituent instead, mirroring `_alcohol.py`'s
    # identical "chain wins" branch (`name_branch`'s `_ring_substituent_
    # with_named_atoms` generalization). PubChem PUG REST-verified
    # "1-(2-aminocyclohexyl)ethane-1,2-diamine" (CID 77647753).
    assert (
        smiles_to_iupac("NC1CCCCC1C(N)CN") == "1-(2-aminocyclohexyl)ethane-1,2-diamine"
    )
