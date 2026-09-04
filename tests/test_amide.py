import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). The amide carbon is
        # always C1, and its own locant is never cited (P-14.3.3).
        ("C(N)=O", "methanamide"),
        ("CC(N)=O", "ethanamide"),
        ("CCC(N)=O", "propanamide"),
        # A real positional choice for a substituent: the amide carbon fixes
        # C1 regardless.
        ("CC(C)CC(N)=O", "3-methylbutanamide"),
        # -amide combined with existing unsaturation support, cross-checked
        # against PubChem (crotonamide's PIN).
        ("CC=CC(N)=O", "but-2-enamide"),
        # -amide + halogen substituent prefix, cross-checked against
        # PubChem.
        ("NC(=O)CCCl", "3-chloropropanamide"),
    ],
)
def test_amide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_methylethanamide():
    # PubChem structure match: "N-methylacetamide" (PubChem uses the
    # retained 'acetamide' stem; this module uses the systematic PIN stem
    # 'ethanamide' consistently with its own plain-amide tests above).
    assert smiles_to_iupac("CNC(C)=O") == "N-methylethanamide"


def test_n_n_dimethylethanamide():
    # PubChem structure match: "N,N-dimethylacetamide".
    assert smiles_to_iupac("CC(=O)N(C)C") == "N,N-dimethylethanamide"


def test_n_ethyl_n_methylethanamide():
    # PubChem structure match: "N-ethyl-N-methylacetamide".
    assert smiles_to_iupac("CC(=O)N(C)CC") == "N-ethyl-N-methylethanamide"


def test_n_substituted_amide_with_longer_acyl_chain():
    # The N-substituent must not be mistaken for the acyl chain even when
    # it happens to be longer than it.
    assert smiles_to_iupac("CC(=O)NCCCC") == "N-butylethanamide"


def test_n_n_disubstituted_amide_with_locant_leading_parent_name():
    # The N,N-prefix must be hyphen-separated from a parent name that
    # itself starts with a numeric locant, not just concatenated. (This
    # SMILES previously used '-CH2CH2Cl' N-substituents, CCCl -- that
    # halogen was invisible to the old length-only carbon-chain check and
    # silently misnamed as plain 'diethyl'; now explicitly rejected (see
    # `test_halogenated_n_substituent_not_supported` below), so this test
    # uses plain ethyl instead to keep testing its actual point, the
    # hyphenation.)
    assert (
        smiles_to_iupac("O=C(C(C(F)(F)F)C(F)(F)F)N(CC)CC")
        == "N,N-diethyl-3,3,3-trifluoro-2-(trifluoromethyl)propanamide"
    )


def test_halogenated_n_substituent_not_supported():
    # A halogen on the N-substituent is out of scope (module docstring:
    # "plain, unsubstituted" N-substituent only) -- previously silently
    # misnamed instead of rejected (see comment above).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NCCCl")


def test_n_substituted_amide_with_locant_leading_parent_name_raises():
    # This SMILES's N-substituent isn't actually a plain propyl group -- it
    # carries two hydroxyls of its own ("CC(O)C(O)N..."), which this
    # module's docstring already scopes out ("plain, unbranched,
    # unsubstituted" N-substituent only). Before `specified_stereocenters`
    # was wired in, this was silently
    # accepted and misnamed as "N-propyl-2,3,4,5,6,7-hexahydroxyheptanamide"
    # -- a genuine pre-existing bug on two fronts: the substituted
    # N-substituent was never rejected, and the acyl chain's five specified
    # stereocenters were silently dropped. It now correctly raises (either
    # via the new stereo safety net noticing the N-substituent's own two
    # stereocenters are left unspecified while the acyl chain's are
    # specified, or via the new explicit N-substituent-hydroxyl check).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(O)C(O)NC(=O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 136874, 12985.
        ("CC(=O)NC(C)C", "N-(propan-2-yl)ethanamide"),
        ("CC(=O)NC(C)(C)C", "N-tert-butylethanamide"),
    ],
)
def test_branched_n_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_identical_branched_n_substituents_parenthesized_when_compound():
    # PubChem CID 69797.
    assert smiles_to_iupac("CC(=O)N(C(C)C)C(C)C") == "N,N-di(propan-2-yl)ethanamide"


def test_two_identical_branched_n_substituents_not_parenthesized_when_retained():
    # PubChem CID 18999412.
    assert smiles_to_iupac("CC(=O)N(C(C)(C)C)C(C)(C)C") == "N,N-ditert-butylethanamide"


def test_two_different_n_substituents_alphabetized_ignoring_italic_prefix():
    # PubChem CID 54197906.
    assert smiles_to_iupac("CC(=O)N(CC)C(C)(C)C") == "N-tert-butyl-N-ethylethanamide"


def test_unsaturated_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC=C")


def test_lactam_is_named_via_ketone_module():
    # A plain, unsubstituted lactam (ketone carbonyl directly bonded to the
    # ring's own N-H) is routed to `_ketone.py`'s hetero-ring ketone path
    # instead of being rejected here -- see test_ketone.py's own
    # coverage. PubChem-verified: CID 12025.
    assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"


def test_n_substituted_lactam_named_via_ketone_module():
    # An N-alkyl lactam is routed to `_ketone.py`'s hetero-ring ketone
    # path (same as the plain lactam above), which supports a single
    # plain alkyl substituent on the ring nitrogen. PubChem-verified:
    # 1-methylpyrrolidin-2-one.
    assert smiles_to_iupac("O=C1CCCN1C") == "1-methylpyrrolidin-2-one"


def test_diamide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)CC(N)=O")


def test_ester_not_misnamed_as_amide():
    # -COO- (ester, `_ester.py`) is not amide-shaped: its carbonyl carbon's
    # other oxygen neighbor is carbon-bonded, not a nitrogen, so this must
    # not be routed here and misnamed.
    assert smiles_to_iupac("CC(=O)OC") == "methyl ethanoate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/
        # #271/#272/#273's carboxylic-acid/ketone/alcohol/ester/aldehyde
        # chains): the ring is cited as a "phenyl" substituent prefix.
        # Cross-checked against PubChem CID 68140 ("3-phenylpropanamide").
        ("c1ccccc1CCC(N)=O", "3-phenylpropanamide"),
        # Two-carbon chain: PubChem's own name for this SMILES uses the
        # retained "phenylacetamide" stem, but this module always uses the
        # systematic 'ethanamide' stem (see test_amide_names above), so
        # this is an accepted, reviewed result rather than a
        # PubChem-confirmed one -- same policy as the aldehyde module's
        # '2-phenylethanal' (PR #273).
        ("c1ccccc1CC(N)=O", "2-phenylethanamide"),
    ],
)
def test_phenyl_chain_amide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_amide():
    # Same benzamide-type case as test_benzamide, reached via a different
    # ring-atom ordering in the SMILES.
    assert smiles_to_iupac("c1ccccc1C(N)=O") == "benzamide"


def test_phenyl_chain_amide_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCC(=O)NC")


def test_phenyl_substituted_benzene_ring_amide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(N)=O")


def test_phenyl_chain_amide_with_hydroxyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCc1ccccc1CC(N)=O")


def test_phenyl_chain_amide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(N)=O")


def test_carboxylic_acid_not_misnamed_as_amide():
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_benzamide():
    # -CONH2 attached directly to a benzene ring carbon (P-66.1.1.1.2.1,
    # one of only four retained amide names that are PINs and can be
    # substituted). The retained name 'benzamide' stands for the whole
    # ring+CONH2 system (like 'benzonitrile'), so the amide's own ring
    # locant is never cited, only other substituents'. PubChem PUG REST:
    # "benzamide"/"4-methylbenzamide"/"2-chlorobenzamide" -- all exact
    # matches.
    assert smiles_to_iupac("NC(=O)c1ccccc1") == "benzamide"
    assert smiles_to_iupac("NC(=O)c1ccc(C)cc1") == "4-methylbenzamide"
    assert smiles_to_iupac("NC(=O)c1ccccc1Cl") == "2-chlorobenzamide"


def test_ring_amide():
    # -CONH2 attached directly to a saturated monocyclic ring carbon
    # (P-66.1.1.1.1.3, the 'carboxamide' suffix) -- e.g.
    # 'cyclohexanecarboxamide'. PubChem PUG REST:
    # "cyclohexanecarboxamide"/"4-methylcyclohexane-1-carboxamide" -- both
    # exact matches.
    assert smiles_to_iupac("NC(=O)C1CCCCC1") == "cyclohexanecarboxamide"
    assert smiles_to_iupac("NC(=O)C1CCC(C)CC1") == "4-methylcyclohexane-1-carboxamide"


def test_alcohol_mix_names_hydroxy_prefix():
    # 'amide' outranks 'ol' in Table 3.3, so a coexisting standalone -OH is
    # cited as the 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("NC(=O)CCO") == "3-hydroxypropanamide"


def test_amide_enol_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(N)=O")


def test_n_substituent_with_hydroxyl_raises():
    # A hydroxyl on the N-substituent itself (not the acyl chain) is
    # invisible to the carbon-only chain walk that measures N-substituent
    # length, so it must be checked separately -- only a plain,
    # unsubstituted alkyl N-substituent is in scope (module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(O)NC(=O)C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): an amide
        # carbon is always chain-terminal (fixed C1), same pattern as
        # `_carboxylic_acid.py`/`_aldehyde.py` (CIP computed entirely by
        # RDKit's `rdCIPLabeler`). PubChem CID 41097729.
        ("CC[C@@H](C)C(N)=O", "(2R)-2-methylbutanamide"),
        ("CC[C@H](C)C(N)=O", "(2S)-2-methylbutanamide"),
    ],
)
def test_amide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_amide_stereocenter():
    # The stereodescriptor sits outermost, ahead of the N-alkyl prefix
    # (P-91.3: always at the very front of the complete name). PubChem
    # CID 89121954.
    assert smiles_to_iupac("CC[C@@H](C)C(=O)NCC") == "(2R)-N-ethyl-2-methylbutanamide"


def test_amide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(N)=O") == "2-methylbutanamide"
