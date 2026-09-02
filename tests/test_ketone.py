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
    # An ether coexisting with a ketone (as opposed to a plain ether on its
    # own, now handled by the separate ether module) is still out of scope
    # for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOCC(=O)C")


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


def test_hetero_ring_ketone_substituted_heteroatom_raises():
    # An N-methyl ring heteroatom is out of scope for this module's narrow
    # first pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCN(C)CC1")


def test_hetero_ring_ketone_ring_substituent_raises():
    # A plain alkyl substituent elsewhere on the ring is out of scope for
    # this module's narrow first pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC(C)NCC1")


def test_hetero_ring_ketone_two_heteroatoms_raises():
    # A ring-fused urea (e.g. hydantoin, two N heteroatoms in a 5-membered
    # ring) doesn't fit the 6-membered 1,4-two-heteroatom shape either, so
    # it falls through to the existing, unrelated urea-module rejection.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CNC(=O)N1")


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


def test_five_membered_1_3_ring_ketone_heteroatoms_adjacent_raises():
    # The two heteroatoms are directly bonded (a 1,2-relationship, not
    # 1,3-) -- out of scope, falls through to whatever other module (if
    # any) matches this shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCNN1")


def test_five_membered_1_3_ring_ketone_extra_carbonyl_raises():
    # Hydantoin (both ring carbons flanking the heteroatoms carbonylated,
    # not just the one directly between them) is out of scope for this
    # module's narrow single-ketone-at-C2 path.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1NC(=O)CN1")


def test_five_membered_1_3_ring_ketone_wrong_element_pair_raises():
    # An O+Se pair has no retained name in this module's scope (Se/Te
    # pairs are out of scope) -- falls through to the existing, unrelated
    # rejection for whatever other module (if any) matches this shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1OCC[Se]1")


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


def test_ketone_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(Cl)C(C)=O") == "3-chloropentan-2-one"


def test_ketone_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](Cl)C(Cl)C(C)=O")
