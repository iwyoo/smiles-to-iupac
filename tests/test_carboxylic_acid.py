import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). The -COOH carbon is
        # always C1, and its own locant is never cited (P-14.3.3).
        ("C(=O)O", "methanoic acid"),
        ("CC(=O)O", "ethanoic acid"),
        ("CCC(=O)O", "propanoic acid"),
        # A real positional choice for a substituent, cross-checked against
        # PubChem: the -COOH carbon fixes C1 regardless.
        ("CC(C)CC(=O)O", "3-methylbutanoic acid"),
        # Dicarboxylic acid: both chain ends are -COOH carbons, multiplying
        # prefix + no locants cited, cross-checked against PubChem (adipic
        # acid's PIN).
        ("OC(=O)CCCCC(=O)O", "hexanedioic acid"),
        # -oic acid combined with existing unsaturation support, cross-checked
        # against PubChem (crotonic acid's PIN).
        ("CC=CC(=O)O", "but-2-enoic acid"),
        # Unsaturated dicarboxylic acid, cross-checked against PubChem
        # (fumaric/maleic acid's PIN).
        ("OC(=O)C=CC(=O)O", "but-2-enedioic acid"),
        # -oic acid + halogen substituent prefix, cross-checked against
        # PubChem.
        ("OC(=O)CCCl", "3-chloropropanoic acid"),
        # Diene acid: two multiplied 'ene' locants ahead of the acid suffix,
        # cross-checked against PubChem (sorbic acid's PIN).
        ("CC=CC=CC(=O)O", "hexa-2,4-dienoic acid"),
    ],
)
def test_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mononuclear_parent_omits_locant():
    # P-14.3.4.2(a): a mononuclear parent's own substituent locant is
    # always '1' and never cited, even for the acyclic single-carbon
    # -COOH chain itself (mirrors the ring path's identical P-14.3.3
    # rule, and `_carboxylic_acid_seleninic_acid.py`'s own precedent for
    # this same rule via a coexisting seleninic acid).
    assert smiles_to_iupac("OC(=O)Cl") == "chloromethanoic acid"


def test_ether_raises():
    # An ether coexisting with a carboxylic acid (as opposed to a plain
    # ether on its own, now handled by the separate ether module) is still
    # out of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOCC(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-65.1.2.2.2: a -COOH directly attached to an otherwise
        # unsubstituted saturated monocyclic ring is named with the
        # 'carboxylic acid' suffix on the ring parent hydride, cross-checked
        # against PubChem (e.g. cyclohexanecarboxylic acid's PIN). The sole
        # substituent's ring locant is omitted (P-14.3.3), same as
        # _sulfonic_acid.py's analogous single-substituent ring case.
        ("C1CCC(CC1)C(=O)O", "cyclohexanecarboxylic acid"),
        ("C1CCCC1C(=O)O", "cyclopentanecarboxylic acid"),
        ("C1CCCCCC1C(=O)O", "cycloheptanecarboxylic acid"),
        ("C1CC1C(=O)O", "cyclopropanecarboxylic acid"),
        # P-65.1.2.2.2 + substituent(s) on a *different* ring atom than the
        # carboxylic acid: the ring is renumbered to give the -COOH carbon
        # locant 1 (never omitted once there's another substituent to
        # locate), then the other substituent(s) get the lowest remaining
        # locant set, cross-checked against PubChem PUG REST.
        ("CC1CCC(CC1)C(=O)O", "4-methylcyclohexane-1-carboxylic acid"),  # CID 20330
        ("ClC1CCC(CC1)C(=O)O", "4-chlorocyclohexane-1-carboxylic acid"),  # CID 12603304
        ("CC1CCC(C)C(C1)C(=O)O", "2,5-dimethylcyclohexane-1-carboxylic acid"),  # CID 14048268
        # A substituent sharing the same ring atom as the -COOH (a
        # quaternary ring carbon) is named like any other ring position.
        ("CC1(CCCCC1)C(=O)O", "1-methylcyclohexane-1-carboxylic acid"),  # verified via PubChem PUG REST
        ("OC(=O)C1(Cl)CCCCC1", "1-chlorocyclohexane-1-carboxylic acid"),  # verified via PubChem PUG REST
    ],
)
def test_ring_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_carboxylic_acid_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCC1C(=O)O")


def test_ring_carboxylic_acid_standalone_hydroxyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCCC1C(=O)O")


def test_ring_carboxylic_acid_intervening_chain_carbon_raises():
    # A -COOH one chain carbon away from the ring (rather than directly
    # attached) is still out of scope for this first slice.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(CC1)CC(=O)O")


def test_ring_carboxylic_acid_two_carboxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC(C(=O)O)CC1")


def test_benzoic_acid():
    # -COOH directly on a benzene ring carbon: 'benzoic acid' is a fully
    # retained name (P-65.1.1), cross-checked against PubChem PUG REST.
    assert smiles_to_iupac("OC(=O)c1ccccc1") == "benzoic acid"  # CID 243


def test_substituted_benzoic_acid():
    # A substituent on a different ring atom than the -COOH: the retained
    # 'benzoic acid' name itself carries no locant, so only the other
    # substituent's position is cited.
    assert smiles_to_iupac("OC(=O)c1ccccc1C") == "2-methylbenzoic acid"  # CID 8373
    assert smiles_to_iupac("OC(=O)c1ccc(C)cc1") == "4-methylbenzoic acid"  # CID 7470


def test_phenyl_substituent_carboxylic_acid():
    # A plain, unsubstituted benzene ring hanging off a chain whose far end
    # carries the sole -COOH (P-65.1.1.2's territory doesn't apply here --
    # the -COOH is on the chain, not the ring itself). PubChem CID 999
    # structure-matches this to "phenylacetic acid", but this module's own
    # already-established convention uses the systematic stem once
    # substituted (see the halogen case just below) rather than the
    # retained 'acetic acid' name -- P-65.1.1 retains 'acetic acid' as PIN
    # only for unsubstituted CH3COOH itself.
    assert smiles_to_iupac("c1ccccc1CC(=O)O") == "2-phenylethanoic acid"


def test_phenyl_substituent_carboxylic_acid_matches_halogen_locant_pattern():
    # Same locant mechanism already verified for a halogen substituent
    # (test_carboxylic_acid_names' "3-chloropropanoic acid" case) applied
    # to a phenyl substituent instead -- one axis changed (substituent
    # identity), locant placement itself already proven correct.
    assert smiles_to_iupac("c1ccccc1CCC(=O)O") == "3-phenylpropanoic acid"


def test_phenyl_substituent_carboxylic_acid_directly_on_ring_is_benzoic_acid():
    # No intervening chain carbon between the ring and the -COOH carbon
    # is really the P-65.1.1 benzoic-acid case (test_benzoic_acid above),
    # reached here through the phenyl-chain module's own dispatch instead
    # of the ring module's.
    assert smiles_to_iupac("c1ccccc1C(=O)O") == "benzoic acid"


def test_phenyl_substituent_carboxylic_acid_substituted_ring_ortho_methyl():
    # A ring methyl substituent is now supported -- PubChem PUG
    # REST-verified "2-(2-methylphenyl)acetic acid".
    assert smiles_to_iupac("Cc1ccccc1CC(=O)O") == "2-(2-methylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_naphthalene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2ccccc2c1CC(=O)O")


def test_pyridine_substituent_carboxylic_acid():
    # A plain heteroaromatic monocycle (P-29.3.4.1) reuses the same
    # phenyl-chain path as a plain benzene ring -- structure-matches
    # PubChem CID 108 ("2-pyridin-3-ylacetic acid").
    assert smiles_to_iupac("OC(=O)Cc1cccnc1") == "2-(pyridin-3-yl)ethanoic acid"


def test_furan_substituent_carboxylic_acid():
    # PubChem CID 75974 ("2-(furan-2-yl)acetic acid").
    assert smiles_to_iupac("OC(=O)Cc1ccco1") == "2-(furan-2-yl)ethanoic acid"


def test_thiophene_substituent_carboxylic_acid():
    # PubChem CID 15970 ("2-thiophen-2-ylacetic acid").
    assert smiles_to_iupac("OC(=O)Cc1cccs1") == "2-(thiophen-2-yl)ethanoic acid"


def test_pyrrole_c_substituent_carboxylic_acid_cites_indicated_hydrogen():
    # Attachment at a ring carbon leaves pyrrole's own N-H tautomer intact,
    # so its indicated hydrogen must still be cited (P-25.7.1.3) -- PubChem
    # CID 4220146 ("2-(1H-pyrrol-2-yl)acetic acid").
    assert smiles_to_iupac("OC(=O)Cc1ccc[nH]1") == "2-(1H-pyrrol-2-yl)ethanoic acid"


def test_pyrrole_n_substituent_carboxylic_acid_no_indicated_hydrogen():
    # Attachment directly at the N-H position consumes that hydrogen
    # itself, so no indicated-hydrogen citation is needed -- PubChem
    # CID 242027 ("2-pyrrol-1-ylacetic acid").
    assert smiles_to_iupac("OC(=O)Cn1cccc1") == "2-(pyrrol-1-yl)ethanoic acid"


def test_heteroaromatic_substituent_carboxylic_acid_directly_on_ring_raises():
    # A -COOH directly on a heteroaromatic ring (e.g. nicotinic acid) needs
    # a different suffix construction than this chain-substituent path --
    # out of scope for this step.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)c1cccnc1")


def test_phenyl_substituent_carboxylic_acid_ring_halogen():
    # The ring's other atoms may each carry a single halogen alongside the
    # chain attachment (PubChem CID 15880 IUPACName "2-(4-chlorophenyl)
    # acetic acid" -- this module's own retained-name convention prefers
    # the systematic 'ethanoic acid' stem once substituted, same as the
    # unsubstituted-ring case above).
    assert smiles_to_iupac("Clc1ccc(CC(=O)O)cc1") == "2-(4-chlorophenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_dihalogen():
    # PubChem CID 88209 "2-(2,4-dichlorophenyl)acetic acid".
    assert smiles_to_iupac("Clc1cc(Cl)ccc1CC(=O)O") == "2-(2,4-dichlorophenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_methyl():
    # PubChem PUG REST-verified "2-(4-methylphenyl)acetic acid" -- a
    # plain methyl ring substituent alongside the chain, same
    # non-halogen-specific mechanism
    # `ring_chain_attachment_with_halogens`/`halogenated_phenyl_substituent`
    # already used for halogens.
    assert smiles_to_iupac("CC1=CC=C(CC(=O)O)C=C1") == "2-(4-methylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_halogen_and_methyl():
    # PubChem "2-(2-chloro-5-methylphenyl)acetic acid" -- halogen and
    # methyl ring substituents mixed on the same ring.
    assert smiles_to_iupac("ClC1=CC=C(C)C=C1CC(=O)O") == "2-(2-chloro-5-methylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_chain_terminal_halogens_not_absorbed():
    # Regression: `longest_branched_chain`'s BFS treated a chain-terminal
    # halogen as an ordinary carbon, silently absorbing it into the chain
    # and shortening the halogen substituent list by one. Correct answer
    # has a 3-carbon chain and 3 fluorines, not a 4-carbon chain with 2.
    assert (
        smiles_to_iupac("c1ccccc1C(F)C(F)(F)C(=O)O")
        == "2,2,3-trifluoro-3-phenylpropanoic acid"
    )


def test_phenyl_substituent_carboxylic_acid_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated acyclic
    # alkyl ring substituent (not just methyl) is now supported, reusing
    # `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CC(=O)O") == "2-(4-ethylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_branched_alkyl():
    # PubChem PUG REST IUPACName match: a branched ring substituent
    # (isopropyl) embeds with no extra inner parens of its own.
    assert smiles_to_iupac("CC(C)c1ccc(cc1)CC(=O)O") == "2-(4-propan-2-ylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_tert_butyl():
    # PubChem PUG REST IUPACName match.
    assert smiles_to_iupac("CC(C)(C)c1ccc(cc1)CC(=O)O") == "2-(4-tert-butylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_branched_acid_chain():
    # Ibuprofen: the acid chain itself branches at the ring-adjacent
    # carbon (`longest_branched_chain`, P-44.3.2) -- the alpha-methyl is
    # absorbed into the parent chain ('propanoic acid') rather than cited
    # as a separate substituent on a shorter 'ethanoic acid', matching
    # PubChem's own PIN exactly (CID 3672).
    assert (
        smiles_to_iupac("CC(C)Cc1ccc(C(C)C(=O)O)cc1")
        == "2-[4-(2-methylpropyl)phenyl]propanoic acid"
    )


def test_amine_coexisting_demotes_to_amino_prefix():
    # glycine (H2N-CH2-COOH): a coexisting primary amine is junior to -COOH
    # (Table 3.3) and demoted to the 'amino' prefix (see
    # _carboxylic_acid_amine.py) rather than rejected. This project's own
    # `_carboxylic_acid.py` always uses the systematic 'ethanoic acid' stem
    # rather than the retained 'acetic acid' PubChem uses for glycine's
    # actual PIN ('2-aminoacetic acid', CID 750) -- see
    # _carboxylic_acid_amine.py's module docstring for that pre-existing,
    # inherited divergence.
    assert smiles_to_iupac("NCC(=O)O") == "2-aminoethanoic acid"


def test_alcohol_mix_names_hydroxy_prefix():
    # '-oic acid' outranks 'ol' in Table 3.3, so a coexisting standalone -OH
    # is cited as the 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("OC(=O)CCO") == "3-hydroxypropanoic acid"


def test_carboxylic_acid_enol_mix_raises():
    # A hydroxyl on a C=C carbon (an enol) is a tautomer of a more senior
    # carbonyl form and out of scope, same as `_alcohol.py`'s own enol check.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_alcohol.py`'s `_name_acyclic_alcohol` (CIP computed
        # entirely by RDKit's `rdCIPLabeler`, not reimplemented here).
        # PubChem CID 2724540, name matches exactly (a rare case where
        # PubChem's auto-generated name equals the PIN).
        ("C[C@@H](Cl)C(=O)O", "(2R)-2-chloropropanoic acid"),
        ("C[C@H](Cl)C(=O)O", "(2S)-2-chloropropanoic acid"),
        # Two specified stereocenters, coexisting with a standalone
        # hydroxyl (P-91.3's own worked example, cited in
        # `_common.py`'s `specified_stereocenters` docstring). PubChem
        # CID 21586112, name matches exactly.
        ("C[C@H](O)[C@H](Cl)C(=O)O", "(2S,3S)-2-chloro-3-hydroxybutanoic acid"),
    ],
)
def test_carboxylic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_carboxylic_acid_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention (see `_common.py`'s `specified_stereocenters` docstring).
    assert smiles_to_iupac("CC(Cl)C(=O)O") == "2-chloropropanoic acid"


def test_carboxylic_acid_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](Cl)C(Cl)C(=O)O")


def test_acyclic_carboxylic_acid_specified_ez_double_bond():
    # Standalone specified C=C E/Z stereo, no stereocenter -- mirrors
    # `_ketone.py`/`_aldehyde.py`'s identical PR #710/#713 case. PubChem
    # CID 637090 confirms the structure (SMILES "C/C=C/C(=O)O"), though
    # PubChem's own auto-generated name is "(E)-but-2-enoic acid" (no
    # locant); this project's existing convention is to always cite the
    # locant explicitly.
    assert smiles_to_iupac("OC(=O)/C=C/C") == "(2E)-but-2-enoic acid"


def test_acyclic_carboxylic_acid_stereocenter_with_ez_double_bond_coexistence():
    # P-91.3: a specified tetrahedral stereocenter and a specified C=C
    # double-bond E/Z element on the same acyclic carboxylic acid chain
    # are cited together in one ascending-locant group, mirroring
    # `_alcohol.py`'s/`_ketone.py`'s/`_aldehyde.py`'s identical case.
    # PubChem CIDs 13383322/76965618 confirm both structures (SMILES
    # "C/C=C/[C@@H](C)C(=O)O" / "C/C=C/[C@H](C)C(=O)O"), though PubChem's
    # own names group by type ("(E,2R)-...") rather than by ascending
    # locant -- the expected locant-ascending form here is the primary
    # source's own P-91.3 rule.
    assert smiles_to_iupac("OC(=O)[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enoic acid"
    assert smiles_to_iupac("OC(=O)[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enoic acid"


def test_cyclic_carboxylic_acid_ring_stereocenter():
    # A specified stereocenter on the ring itself, alongside the ring's
    # sole -COOH substituent (P-92), same pattern as `_sulfonic_acid.py`'s
    # `_name_cyclic_sulfonic_acid`/`_ketone.py`'s `_name_cyclic_ketone`.
    # PubChem CID 92338088 (name carries a redundant 'trans-' relative
    # descriptor this project drops once full R/S is given, same policy
    # as the existing ring-alcohol/sulfonic-acid stereocenter tasks).
    assert (
        smiles_to_iupac("OC(=O)[C@H]1CCCC[C@@H]1Cl")
        == "(1R,2S)-2-chlorocyclohexane-1-carboxylic acid"
    )


def test_carboxylic_acid_ring_branch_stereocenter_raises():
    # A stereocenter on a substituent branch elsewhere on the ring (not
    # the -COOH carbon's own ring atom) is out of scope for now -- only
    # a stereocenter on the ring itself is supported (see
    # `tasks/carbo-suffix-ring-stereocenter.md` stage 1; stage 2 will
    # extend this to the ring's sole substituent branch, mirroring
    # `_sulfonic_acid.py`'s `_ring_branch_stereo_display`).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)[C@H]1CC[C@H](C[C@@H](C)CC)C1")


def test_carboxylic_acid_ring_stereocenter_unspecified_unaffected():
    # A genuine ring stereocenter left unspecified is named exactly as
    # before -- no stereo prefix.
    assert smiles_to_iupac("OC(=O)C1CCCCC1Cl") == "2-chlorocyclohexane-1-carboxylic acid"


def test_carboxylic_acid_ring_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(O)[C@H]1CCCCC1Cl")


def test_phenyl_substituent_carboxylic_acid_branch_tie_prefers_more_substituents():
    # Both the phenyl-bearing carbon and the methyl-bearing carbon sit one
    # bond from C2 -- tied for farthest from the -COOH carbon
    # (`longest_branched_chain`). P-44.3.2's own next tie-break after
    # chain length prefers the chain giving more substituents cited as
    # prefixes: 'methyl' + 'phenyl' (two) beats one compound 'benzyl'-
    # shaped substituent, confirmed against PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1CC(C)C(=O)O") == "2-methyl-3-phenylpropanoic acid"
