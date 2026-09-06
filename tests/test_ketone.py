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
        # Monocyclic ring, single ketone, single ring double bond
        # (P-31.1.3): the ketone always gets locant 1 (suffix priority),
        # the ring double bond's locant is minimized by choosing direction.
        # Cross-checked against PubChem (CID 13594/77727).
        ("O=C1CCCC=C1", "cyclohex-2-en-1-one"),
        ("O=C1CC=CCC1", "cyclohex-3-en-1-one"),
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


def test_carboxylic_acid_not_misread_as_ketone():
    # A carboxylic acid is routed to the dedicated carboxylic-acid module
    # (see test_carboxylic_acid.py) instead of falling through here.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_aryl_ketone_raises():
    # A carbonyl on an aromatic ring (an aryl ketone) is out of scope for
    # this module (separate, in-progress aromatic-ring module's territory).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)c1ccccc1")


def test_phenyl_substituent_ketone():
    # A plain, unsubstituted benzene ring hanging off a chain carrying the
    # sole ketone, mirroring `_carboxylic_acid.py`'s identical first slice
    # (test_phenyl_substituent_carboxylic_acid). PubChem CID 7678.
    assert smiles_to_iupac("CC(=O)Cc1ccccc1") == "1-phenylpropan-2-one"


def test_phenyl_substituent_ketone_matches_locant_mechanism():
    # Same locant mechanism, one carbon further from the ring. PubChem
    # CID 17355.
    assert smiles_to_iupac("CC(=O)CCc1ccccc1") == "4-phenylbutan-2-one"


def test_phenyl_substituent_ketone_directly_on_ring_raises():
    # No intervening chain carbon between the ring and the carbonyl carbon
    # is the aryl-ketone case, already covered by test_aryl_ketone_raises --
    # this is the same rejection reached through the new benzene-ring code
    # path instead of falling through to the old blanket ring rejection.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)C")


def test_phenyl_substituent_ketone_substituted_ring_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout.md) -- PubChem PUG
    # REST-verified "1-(2-methylphenyl)propan-2-one".
    assert smiles_to_iupac("Cc1ccccc1CC(C)=O") == "1-(2-methylphenyl)propan-2-one"


def test_phenyl_substituent_ketone_ring_halogen():
    # PubChem PUG REST computes "1-(4-chlorophenyl)propan-2-one" for this
    # structure, matching the systematic form directly (no retained-name
    # divergence here, unlike the acid/nitrile/aldehyde cases).
    assert smiles_to_iupac("CC(=O)Cc1ccc(Cl)cc1") == "1-(4-chlorophenyl)propan-2-one"


def test_phenyl_substituent_ketone_ring_dihalogen():
    assert smiles_to_iupac("Clc1cc(Cl)ccc1CC(C)=O") == "1-(2,4-dichlorophenyl)propan-2-one"


def test_phenyl_substituent_ketone_ring_methyl():
    # PubChem PUG REST-verified "1-(4-methylphenyl)propan-2-one".
    assert smiles_to_iupac("Cc1ccc(CC(C)=O)cc1") == "1-(4-methylphenyl)propan-2-one"


def test_phenyl_substituent_ketone_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CC(=O)C") == "1-(4-ethylphenyl)propan-2-one"


def test_phenyl_substituent_ketone_naphthalene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2ccccc2c1CC(C)=O")


def test_bicyclic_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2CCC1CC2")


def test_ether_now_supported_via_ether_ketone():
    # A coexisting ether is now handled by `_ether_ketone.py` (P-41: an
    # ether has no suffix at all, so it's always the 'alkoxy' prefix) --
    # see tests/test_ether_ketone.py for that module's own coverage.
    assert smiles_to_iupac("CCOCC(=O)C") == "1-ethoxypropan-2-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Single-heteroatom saturated ring, single ketone -- PubChem-verified
        # (see the module docstring). The heteroatom is always locant 1; numbering direction
        # is chosen to minimize the ketone locant(s), same P-44.4.1.8 rule
        # as the carbocyclic path.
        ("O=C1CCNCC1", "piperidin-4-one"),
        ("O=C1CCCNC1", "piperidin-3-one"),
        ("O=C1CCOCC1", "oxan-4-one"),
        ("O=C1CCSCC1", "thian-4-one"),
        # Lactams (the ketone carbonyl directly bonded to the ring N) fit
        # this same shape and are routed here ahead of `_amide.py`'s
        # ring-always-out-of-scope guard (see test_amide.py's
        # test_lactam_is_named_via_ketone_module). PubChem-verified: CID
        # 12665/12025/7768.
        ("O=C1CCCCN1", "piperidin-2-one"),
        ("O=C1CCCN1", "pyrrolidin-2-one"),
        ("O=C1CCCCCN1", "azepan-2-one"),
        # A symmetric, unsubstituted cyclic imide (succinimide) has both
        # ring carbonyls bonded to the same N and is routed here too (see
        # test_imide.py's test_cyclic_symmetric_imide_is_named_via_ketone_suffix) --
        # P-66.6.3's "N-acyl amide" imide construction only covers the
        # acyclic case. PubChem-verified: CID 11439.
        ("O=C1CCC(=O)N1", "pyrrolidine-2,5-dione"),
        # Lactones (the ketone carbonyl directly bonded to the ring O/S)
        # fit this same shape and are routed here ahead of `_ester.py`'s
        # acyclic-only construction.
        # PubChem-verified: CID 10953/10953473.
        ("O=C1CCCCO1", "oxan-2-one"),
        ("O=C1CCCCS1", "thian-2-one"),
    ],
)
def test_hetero_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # N-alkyl on the single-heteroatom ring nitrogen -- always locant
        # 1, so the N-substituent prefix carries locant 1 too. PubChem-
        # verified: 1-methylpyrrolidin-2-one (CID 12025 analog), and the
        # piperidin-2-one/azepan-2-one homologs (see the module's
        # unsubstituted lactam cases above for the parent ring names).
        ("O=C1CCCN1C", "1-methylpyrrolidin-2-one"),
        ("O=C1CCCCN1C", "1-methylpiperidin-2-one"),
        ("O=C1CCCCCN1CC", "1-ethylazepan-2-one"),
        # An N-alkyl heteroatom not adjacent to the ketone carbonyl is
        # equally in scope -- the fix isn't lactam-specific.
        ("O=C1CCN(C)CC1", "1-methylpiperidin-4-one"),
    ],
)
def test_hetero_ring_ketone_n_alkyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected




def test_hetero_ring_ketone_ring_substituent_raises():
    # A plain alkyl substituent elsewhere on the ring is out of scope for
    # this module's narrow first pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC(C)NCC1")


def test_hetero_ring_ketone_two_heteroatoms_wrong_ring_size():
    # A ring-fused urea (hydantoin, two N heteroatoms in a 5-membered
    # ring) doesn't fit the 6-membered 1,4-two-heteroatom shape, but is
    # named correctly via the separate five-membered 1,3-two-heteroatom
    # dione path (see test_five_membered_1_3_ring_dione_names) instead of
    # this module's 6-membered path.
    assert smiles_to_iupac("O=C1CNC(=O)N1") == "imidazolidine-2,4-dione"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,4-related two-heteroatom 6-membered saturated ring, one of the
        # six retained/systematic-name shapes (morpholine/piperazine/
        # thiomorpholine/1,4-dioxane/1,4-oxathiane/1,4-dithiane) --
        # PubChem-verified. The higher-priority
        # heteroatom (P-22.2.1 order O > S > N) is always locant 1; for
        # piperazine's two identical nitrogens, 1,4-dioxane's two
        # identical oxygens, and 1,4-dithiane's two identical sulfurs,
        # both are tried as the locant-1 candidate.
        ("O=C1COCCN1", "morpholin-3-one"),
        ("O=C1CNCCN1", "piperazin-2-one"),
        ("O=C1CSCCN1", "thiomorpholin-3-one"),
        ("O=C1COCCO1", "1,4-dioxan-2-one"),
        ("O=C1COCCS1", "1,4-oxathian-3-one"),
        ("O=C1CSCCS1", "1,4-dithian-2-one"),
        # Symmetric diketone: same multiplying-prefix + full-locant-set
        # pattern as every other hetero-ring ketone case.
        ("O=C1CNC(=O)CN1", "piperazine-2,5-dione"),
        ("O=C1COC(=O)CO1", "1,4-dioxane-2,5-dione"),
        ("O=C1CSC(=O)CS1", "1,4-dithiane-2,5-dione"),
    ],
)
def test_two_hetero_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_hetero_ring_ketone_substituted_heteroatom_raises():
    # An N-methyl ring heteroatom is out of scope, same as the
    # single-heteroatom path.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1COCCN1C")


def test_two_hetero_ring_ketone_wrong_element_pair_raises():
    # An O+Se pair has no retained name in this module's scope (Se/Te
    # pairs are out of scope) --
    # falls through to the existing, unrelated rejection for whatever
    # other module (if any) matches this shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1COCC[Se]1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,4-related two-heteroatom 7-membered saturated ring, one of the
        # six retained/systematic-name shapes (1,4-diazepane/1,4-oxazepane/
        # 1,4-thiazepane/1,4-dioxepane/1,4-oxathiepane/1,4-dithiepane) --
        # PubChem-verified. Unlike the 6-membered ring, the two arcs
        # between the heteroatoms differ in length (2 carbons vs. 3), so
        # only one numbering direction per candidate start actually gives
        # locant 4 to the other heteroatom.
        ("C1CNCCNC1=O", "1,4-diazepan-5-one"),
        ("C1CNCC(=O)NC1", "1,4-diazepan-2-one"),
        ("C1COCCNC1=O", "1,4-oxazepan-5-one"),
        ("C1CSCCNC1=O", "1,4-thiazepan-5-one"),
        ("C1COCC(=O)OC1", "1,4-dioxepan-2-one"),
        ("C1COC(=O)CSC1", "1,4-oxathiepan-2-one"),
    ],
)
def test_seven_membered_1_4_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_membered_1_4_ring_ketone_wrong_element_pair_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CNCC[Se]C1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,3-related two-heteroatom 7-membered saturated ring, one of the
        # six retained/systematic-name shapes (1,3-diazepane/1,3-oxazepane/
        # 1,3-thiazepane/1,3-dioxepane/1,3-oxathiepane/1,3-dithiepane) --
        # PubChem-verified. The sole ketone always sits on the bridging
        # carbon between the two heteroatoms (locant 2), regardless of
        # numbering direction, same as the 5-membered 1,3-ring shape.
        ("C1CCNC(=O)NC1", "1,3-diazepan-2-one"),
        ("C1CCOC(=O)NC1", "1,3-oxazepan-2-one"),
        ("C1CCSC(=O)NC1", "1,3-thiazepan-2-one"),
        ("C1CCOC(=O)OC1", "1,3-dioxepan-2-one"),
        ("C1CCSC(=O)OC1", "1,3-oxathiepan-2-one"),
        ("C1CCSC(=O)SC1", "1,3-dithiepan-2-one"),
    ],
)
def test_seven_membered_1_3_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_membered_1_3_ring_ketone_second_ketone_on_far_arc_raises():
    # A second ketone anywhere on the 4-carbon far arc is out of scope --
    # no PubChem-registered example was found to confirm its locant.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1NC(=O)CCCN1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,2-related (directly bonded) two-heteroatom 7-membered
        # saturated ring, one of the six retained/systematic-name shapes
        # (1,2-diazepane/1,2-oxazepane/1,2-thiazepane/1,2-dioxepane/
        # 1,2-oxathiepane/1,2-dithiepane) -- PubChem-verified by
        # connectivity (its own computed name drops the '1,2-' locant for
        # all six pairs, contradicted by the same P-22.2.2.1.7/stem-clash
        # evidence already used for the unsubstituted parent names, so
        # '1,2-' is kept here too). The sole ketone always sits on the
        # ring carbon immediately after the second heteroatom (locant 3).
        ("N1NC(=O)CCCC1", "1,2-diazepan-3-one"),
        ("O1NC(=O)CCCC1", "1,2-oxazepan-3-one"),
        ("S1NC(=O)CCCC1", "1,2-thiazepan-3-one"),
        ("O1OC(=O)CCCC1", "1,2-dioxepan-3-one"),
        ("O1SC(=O)CCCC1", "1,2-oxathiepan-3-one"),
        ("S1SC(=O)CCCC1", "1,2-dithiepan-3-one"),
    ],
)
def test_seven_membered_1_2_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_seven_membered_1_2_ring_ketone_second_ketone_on_far_arc_raises():
    # A second ketone anywhere on the 4-carbon far arc is out of scope --
    # no PubChem-registered example was found to confirm its locant.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N1NC(=O)CCC(=O)C1")


def test_seven_membered_1_2_ring_ketone_wrong_element_pair_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC[Se]N1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A ketone carbonyl sitting directly between the two heteroatoms of
        # a five-membered 1,3-related saturated ring -- always locant 2,
        # regardless of numbering direction, since it's the only ring atom
        # between them on the short arc. PubChem-verified: imidazolidin-2-one
        # (CID 8453), 1,3-oxazolidin-2-one (CID 73949), 1,3-thiazolidin-2-one
        # (CID 97431), 1,3-dioxolan-2-one (CID 7303), 1,3-oxathiolan-2-one
        # (CID 72822), 1,3-dithiolan-2-one (CID 123140).
        ("O=C1NCCN1", "imidazolidin-2-one"),
        ("O=C1OCCN1", "1,3-oxazolidin-2-one"),
        ("O=C1SCCN1", "1,3-thiazolidin-2-one"),
        ("O=C1OCCO1", "1,3-dioxolan-2-one"),
        ("O=C1SCCO1", "1,3-oxathiolan-2-one"),
        ("O=C1SCCS1", "1,3-dithiolan-2-one"),
    ],
)
def test_five_membered_1_3_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_3_ring_ketone_heteroatoms_adjacent_not_this_path():
    # The two heteroatoms are directly bonded (a 1,2-relationship, not
    # 1,3-) -- out of scope for this module's 1,3-path, but named
    # correctly via the separate five-membered 1,2-path instead (see
    # test_five_membered_1_2_ring_ketone_names).
    assert smiles_to_iupac("O=C1CCNN1") == "pyrazolidin-3-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A second ketone on the long arc of the same five-membered
        # 1,3-ring shape (the hydantoin family) -- when the two
        # heteroatoms differ, P-22.2.1 element seniority (O > S > N) fixes
        # locant 1 regardless of which locant that gives the second
        # ketone (no minimization); when they're identical, either may be
        # locant 1, so the lower ketone locant set wins. PubChem-verified:
        # imidazolidine-2,4-dione (CID 10006), 1,3-oxazolidine-2,5-dione
        # (CID 75136), 1,3-thiazolidine-2,5-dione (CID 542718),
        # 1,3-dioxolane-2,4-dione (CID 12793796),
        # 1,3-oxathiolane-2,4-dione (CID 67415624),
        # 1,3-dithiolane-2,4-dione (CID 637793).
        ("O=C1NC(=O)CN1", "imidazolidine-2,4-dione"),
        ("O=C1OC(=O)CN1", "1,3-oxazolidine-2,5-dione"),
        ("O=C1SC(=O)CN1", "1,3-thiazolidine-2,5-dione"),
        ("O=C1OC(=O)CO1", "1,3-dioxolane-2,4-dione"),
        ("O=C1SC(=O)CO1", "1,3-oxathiolane-2,4-dione"),
        ("O=C1SC(=O)CS1", "1,3-dithiolane-2,4-dione"),
    ],
)
def test_five_membered_1_3_ring_dione_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unbranched, unsubstituted alkyl substituent on a ring
        # nitrogen of the five-membered 1,3-two-heteroatom shape -- its
        # own locant is whatever the already-fixed ketone-locant numbering
        # gives it (P-44.4.1.8: suffix locants are decided first and win
        # outright), not separately minimized. PubChem-verified:
        # 3-methylimidazolidine-2,4-dione (CID 138851),
        # 1-methylimidazolidine-2,4-dione (CID 69217),
        # 1,3-dimethylimidazolidine-2,4-dione (CID 123410),
        # 1-methylimidazolidin-2-one (CID 567600, same pattern on the
        # single-ketone shape).
        ("O=C1N(C)C(=O)CN1", "3-methylimidazolidine-2,4-dione"),
        ("O=C1NC(=O)CN1C", "1-methylimidazolidine-2,4-dione"),
        ("O=C1N(C)C(=O)CN1C", "1,3-dimethylimidazolidine-2,4-dione"),
        ("O=C1N(C)CCN1", "1-methylimidazolidin-2-one"),
    ],
)
def test_five_membered_1_3_ring_ketone_n_substituent_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_3_ring_ketone_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1N(C(C)C)C(=O)CN1")


def test_five_membered_1_3_ring_ketone_substituted_n_substituent_raises():
    # A hydroxyl on the N-substituent itself -- only a plain, unbranched,
    # unsubstituted alkyl N-substituent is in scope (same restriction as
    # `_amide.py`/`_hydrazide.py`).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1N(CCO)C(=O)CN1")


def test_five_membered_1_3_ring_ketone_single_ketone_ring_carbon_substituent():
    # Different heteroatom elements fix the numbering direction outright
    # (O outranks N for locant 1), so the substituent locant just falls
    # out of that -- no tie-break needed. PubChem CID 136837.
    assert smiles_to_iupac("O=C1OC(C)CN1") == "5-methyl-1,3-oxazolidin-2-one"


def test_five_membered_1_3_ring_ketone_single_ketone_ring_carbon_substituent_tie_break():
    # Identical heteroatoms (both N): the lone ketone sits on the bridging
    # carbon (locant 2) regardless of direction, so the two candidate
    # numberings tie on the suffix locant and fall through to P-14.5.2's
    # next criterion -- lowest locant to the substituent set, giving '4-'
    # rather than '5-'. PubChem CID 97832.
    assert smiles_to_iupac("O=C1NC(C)CN1") == "4-methylimidazolidin-2-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # An alkyl substituent on the sole remaining plain ring carbon of
        # the dione (hydantoin) shape -- unlike the single-ketone case
        # above, the ketone locants alone already fully fix the
        # numbering, so the substituent's own locant is simply whatever
        # that fixed numbering gives it (same P-44.4.1.8 principle as the
        # N-substituent case). PubChem-verified:
        # 5-methylimidazolidine-2,4-dione (CID 69216),
        # 4-methyl-1,3-oxazolidine-2,5-dione (CID 70938, its locant '4'
        # matches the same carbon's locant in the unsubstituted
        # numbering), 5-methyl-1,3-dioxolane-2,4-dione (CID 22227598).
        ("O=C1NC(=O)C(C)N1", "5-methylimidazolidine-2,4-dione"),
        ("O=C1OC(=O)C(C)N1", "4-methyl-1,3-oxazolidine-2,5-dione"),
        ("O=C1OC(=O)C(C)O1", "5-methyl-1,3-dioxolane-2,4-dione"),
    ],
)
def test_five_membered_1_3_ring_dione_carbon_substituent_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_3_ring_ketone_wrong_element_pair_raises():
    # An O+Se pair has no retained name in this module's scope (Se/Te
    # pairs are out of scope) -- falls through to the existing, unrelated
    # rejection for whatever other module (if any) matches this shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1OCC[Se]1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A ketone carbonyl on the five-membered 1,2-two-heteroatom ring
        # shape (the two heteroatoms directly bonded, locants 1/2 fixed;
        # the ketone lands on one of the three remaining ring carbons,
        # locants 3/4/5). PubChem-verified for the N-containing pairs:
        # pyrazolidin-3-one (CID 151497), 1,2-oxazolidin-3-one (CID
        # 192737), 1,2-thiazolidin-3-one (CID 21878697). The O/S-only
        # pairs follow the same '1,2-' locant-citation rule established
        # in test_five_membered_1_2_two_heteroatom_ring_names (Blue Book
        # Table 2.3/P-22.2.2.1.3, not PubChem, which drops the locant for
        # these three specifically): 1,2-dioxolan-3-one,
        # 1,2-oxathiolan-3-one, 1,2-dithiolan-3-one.
        ("O=C1CCNN1", "pyrazolidin-3-one"),
        ("O=C1CCON1", "1,2-oxazolidin-3-one"),
        ("O=C1CCSN1", "1,2-thiazolidin-3-one"),
        ("O=C1CCOO1", "1,2-dioxolan-3-one"),
        ("O=C1CCOS1", "1,2-oxathiolan-3-one"),
        ("O=C1CCSS1", "1,2-dithiolan-3-one"),
    ],
)
def test_five_membered_1_2_ring_ketone_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_2_ring_ketone_wrong_element_pair_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC[Se]N1")


def test_five_membered_1_2_ring_ketone_single_ketone_ring_carbon_substituent():
    # Unlike the 1,3-ring case, the heteroatom-then-ketone locant rules
    # already fully determine the numbering direction even with a single
    # ketone (the two candidate directions give the ketone locant 3 vs 5,
    # never a tie), so no separate substituent tie-break is needed here.
    # PubChem CID 312666.
    assert smiles_to_iupac("CC1CC(=O)NN1") == "5-methylpyrazolidin-3-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # An alkyl substituent on the sole remaining plain ring carbon of
        # the 1,2-dione shape -- mirrors the 1,3-ring dione case above.
        # PubChem-verified: 4-methylpyrazolidine-3,5-dione (CID 12391681).
        # The O/S-containing pairs follow the same '1,2-' locant-citation
        # rule as their unsubstituted parents (not independently
        # PubChem-verified for the substituted form, same evidentiary
        # gap as the unsubstituted O/S-only names above).
        ("O=C1C(C)C(=O)NN1", "4-methylpyrazolidine-3,5-dione"),
        ("O=C1C(C)C(=O)ON1", "4-methyl-1,2-oxazolidine-3,5-dione"),
        ("O=C1C(C)C(=O)SN1", "4-methyl-1,2-thiazolidine-3,5-dione"),
    ],
)
def test_five_membered_1_2_ring_dione_carbon_substituent_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_alcohol_hetero_mix_names_hydroxy_prefix():
    # 'one' outranks 'ol' in Table 3.3, so a coexisting -OH is cited as the
    # 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("OCC(=O)C") == "1-hydroxypropan-2-one"


def test_ketone_alcohol_mix_on_longer_chain():
    assert smiles_to_iupac("CC(=O)CCO") == "4-hydroxybutan-2-one"


def test_cyclic_ketone_alcohol_mix_names_hydroxy_prefix():
    assert smiles_to_iupac("OC1CCC(=O)CC1") == "4-hydroxycyclohexan-1-one"


def test_ketone_enol_mix_raises():
    # A hydroxyl on a C=C carbon (an enol) is a tautomer of a more senior
    # carbonyl form and out of scope, same as `_alcohol.py`'s own enol check.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(=O)C")


def test_unsaturated_ring_ketone_with_substituent_raises():
    # A substituent alongside both a ring double bond and a ketone needs
    # more careful numbering-priority verification than this first pass
    # covers.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC=C1C")


def test_unsaturated_ring_ketone_triple_bond_raises():
    # A ring triple bond (cycloalkyne) alongside a ketone is out of scope
    # for this first pass -- only a ring double bond is supported.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCC#C1")


def test_ring_substituent_chain_ketone():
    # A ketone entirely on a chain hanging off a plain saturated ring (the
    # ring itself bears no ketone) -- the ring is cited as a "cyclo..."
    # substituent prefix on the chain, mirroring `_alcohol.py`'s
    # `_name_ring_substituent_chain_alcohol`. PubChem PUG REST-verified
    # "1-cyclohexylethanone" (CID 13207).
    assert smiles_to_iupac("CC(=O)C1CCCCC1") == "1-cyclohexylethan-1-one"


def test_ring_substituent_chain_ketone_longer_chain():
    # PubChem PUG REST: "1-cyclohexylpropan-1-one" (CID 70748) -- an
    # exact match here since the '1' locant is a genuine choice (three or
    # more chain carbons), unlike the two-carbon case above.
    assert smiles_to_iupac("CCC(=O)C1CCCCC1") == "1-cyclohexylpropan-1-one"


def test_ring_substituent_chain_ketone_ring_with_substituent_raises():
    # A ring atom other than the chain attachment carrying its own
    # substituent is out of scope for this first pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)C1CCC(C)CC1")


def test_ring_substituent_chain_ketone_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)C1CCCC=C1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter on the chain (P-92),
        # same pattern as `_carboxylic_acid.py`/`_aldehyde.py` (CIP
        # computed entirely by RDKit's `rdCIPLabeler`). PubChem CID
        # 92284545.
        ("CC[C@@H](Cl)C(C)=O", "(3R)-3-chloropentan-2-one"),
        ("CC[C@H](Cl)C(C)=O", "(3S)-3-chloropentan-2-one"),
        # Two specified stereocenters, ascending-locant group (P-91.3).
        # PubChem CID 12688076, name matches exactly.
        ("C[C@H](Cl)[C@H](Cl)C(C)=O", "(3R,4S)-3,4-dichloropentan-2-one"),
    ],
)
def test_acyclic_ketone_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_ketone_stereocenter():
    # A stereocenter on the ring itself (P-92), same pattern as
    # `_alcohol.py`'s `_name_cyclic_alcohol`. PubChem CID 641138.
    assert smiles_to_iupac("O=C1CCCC[C@H]1Cl") == "(2R)-2-chlorocyclohexan-1-one"


def test_cyclic_ketone_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_alcohol.py`'s
    # `_name_cyclic_alcohol` -- the branch sits at C4 (para to the
    # carbonyl), a ring position symmetric enough that the ring atom
    # itself needs no locant/stereo label of its own. PubChem's own
    # computed IUPACName for this exact SMILES agrees.
    assert smiles_to_iupac("O=C1CCC(CC1)[C@@H](C)CC") == "4-[(2S)-butan-2-yl]cyclohexan-1-one"


def test_ketone_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(Cl)C(C)=O") == "3-chloropentan-2-one"


def test_ketone_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](Cl)C(Cl)C(C)=O")


def test_phenyl_substituent_ketone_branched_chain():
    # A branch at the ring-adjacent carbon is absorbed into the parent
    # chain (`longest_branched_chain_through`, P-44.3.2), same principle
    # as `_carboxylic_acid.py`'s ibuprofen case, but here the ketone
    # carbon can sit anywhere along the resulting chain rather than
    # always at C1. PubChem PUG REST-verified "3-phenylbutan-2-one".
    assert smiles_to_iupac("c1ccccc1C(C)C(=O)C") == "3-phenylbutan-2-one"
