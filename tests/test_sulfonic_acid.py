import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanesulfonic_acid():
    # '-sulfonic acid' suffix construction mirrors '-thiol'/'-ol' locant
    # rules (see test_thiol.py/test_alcohol.py). 'methanesulfonic acid' is
    # also a well-known real compound (a common industrial acid catalyst),
    # independently verifiable.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_ethanesulfonic_acid():
    assert smiles_to_iupac("CCS(=O)(=O)O") == "ethanesulfonic acid"


def test_propane_1_sulfonic_acid():
    assert smiles_to_iupac("CCCS(=O)(=O)O") == "propane-1-sulfonic acid"


def test_propane_2_sulfonic_acid():
    assert smiles_to_iupac("CC(S(=O)(=O)O)C") == "propane-2-sulfonic acid"


def test_chlorobutanesulfonic_acid():
    assert smiles_to_iupac("ClCCCCS(=O)(=O)O") == "4-chlorobutane-1-sulfonic acid"


def test_pent_4_ene_1_sulfonic_acid():
    assert smiles_to_iupac("C=CCCCS(=O)(=O)O") == "pent-4-ene-1-sulfonic acid"


def test_disulfonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CCS(=O)(=O)O")


def test_cyclohexanesulfonic_acid():
    # PubChem CID 428836.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCCC1") == "cyclohexanesulfonic acid"


def test_unsaturated_ring_sulfonic_acid():
    # Monocyclic ring, single -SO3H, single ring double bond (P-31.1.3):
    # the sulfonic acid always gets locant 1 (suffix priority), the ring
    # double bond's locant is minimized by choosing direction.
    # Cross-checked against PubChem (CID 15311747/20471403).
    assert smiles_to_iupac("OS(=O)(=O)C1CCCC=C1") == "cyclohex-2-ene-1-sulfonic acid"
    assert smiles_to_iupac("OS(=O)(=O)C1CC=CCC1") == "cyclohex-3-ene-1-sulfonic acid"


def test_unsaturated_ring_sulfonic_acid_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC=C1C")


def test_unsaturated_ring_sulfonic_acid_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CCCC#C1")


def test_2_methylcyclohexane_1_sulfonic_acid():
    # PubChem CID 121004858.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCCC1C") == "2-methylcyclohexane-1-sulfonic acid"


def test_cyclopentanesulfonic_acid():
    # PubChem CID 15707015.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCC1") == "cyclopentanesulfonic acid"


def test_2_chlorocyclohexane_1_sulfonic_acid():
    # PubChem CID 129994011.
    assert smiles_to_iupac("OS(=O)(=O)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfonic acid"


def test_polycyclic_sulfonic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)C1CC2CCC1CC2")


def test_ring_substituent_chain_sulfonic_acid():
    # A sulfonic acid entirely on a chain hanging off a plain saturated
    # ring (the ring itself bears no sulfonic acid) -- the ring is cited
    # as a "cyclo..." substituent prefix on the chain, mirroring
    # `_name_phenyl_chain_sulfonic_acid`/`_ketone.py`'s `_name_ring_
    # substituent_chain_ketone`. PubChem PUG REST-verified
    # "cyclohexylmethanesulfonic acid" (CID 18406397).
    assert smiles_to_iupac("OS(=O)(=O)CC1CCCCC1") == "cyclohexylmethanesulfonic acid"


def test_ring_substituent_chain_sulfonic_acid_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CC1CCC(C)CC1")


def test_ring_substituent_chain_sulfonic_acid_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CC1CCCC=C1")


def test_sulfonic_acid_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CCO")


def test_phenyl_chain_sulfonic_acid():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#276's
    # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/acyl-
    # halide chains): the ring is cited as a "phenyl" substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCCS(=O)(=O)O") == "3-phenylpropane-1-sulfonic acid"


def test_phenyl_chain_sulfonic_acid_internal_locant():
    # The -SO3H locant is a genuine choice on the chain (unlike a
    # terminus-only suffix like -CHO), same as the base acyclic module.
    assert smiles_to_iupac("c1ccccc1CC(C)S(=O)(=O)O") == "1-phenylpropane-2-sulfonic acid"


def test_benzenesulfonic_acid():
    # -SO3H directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1S(=O)(=O)O") == "benzenesulfonic acid"  # CID 7371


def test_substituted_benzenesulfonic_acid():
    # A substituent on a different ring atom than the -SO3H: the
    # mancude-ring numbering is free to start at the -SO3H carbon, so its
    # own locant is never cited, unlike the cycloalkane case.
    assert smiles_to_iupac("Cc1ccccc1S(=O)(=O)O") == "2-methylbenzenesulfonic acid"  # CID 6925
    assert smiles_to_iupac("Cc1ccc(cc1)S(=O)(=O)O") == "4-methylbenzenesulfonic acid"  # CID 6101


def test_phenyl_substituted_benzene_ring_sulfonic_acid_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-2.md) -- PubChem PUG
    # REST-verified "2-(2-methylphenyl)ethanesulfonic acid".
    assert smiles_to_iupac("Cc1ccccc1CCS(=O)(=O)O") == "2-(2-methylphenyl)ethane-1-sulfonic acid"


def test_phenyl_chain_sulfonic_acid_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propane-1-sulfonic
    # acid" for this structure.
    assert smiles_to_iupac("Clc1ccc(CCCS(=O)(=O)O)cc1") == "3-(4-chlorophenyl)propane-1-sulfonic acid"


def test_phenyl_chain_sulfonic_acid_ring_methyl():
    # PubChem PUG REST-verified "3-(4-methylphenyl)propane-1-sulfonic acid".
    assert smiles_to_iupac("Cc1ccc(CCCS(=O)(=O)O)cc1") == "3-(4-methylphenyl)propane-1-sulfonic acid"


def test_phenyl_chain_sulfonic_acid_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CS(=O)(=O)O") == "(4-ethylphenyl)methanesulfonic acid"


def test_phenyl_chain_sulfonic_acid_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCS(=O)(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Structure verified against PubChem: CID 71650202.
        ("OS(=O)(=O)Cc1cccnc1", "(pyridin-3-yl)methanesulfonic acid"),
        # Structure verified against PubChem: CID 23361618.
        ("OS(=O)(=O)Cc1cccs1", "(thiophen-2-yl)methanesulfonic acid"),
        # Structure verified against PubChem: CID 18323441.
        ("OS(=O)(=O)Cc1ccco1", "(furan-2-yl)methanesulfonic acid"),
        # Structure verified against PubChem: CID 91804704.
        ("OS(=O)(=O)Cc1cc[nH]c1", "(1H-pyrrol-3-yl)methanesulfonic acid"),
    ],
)
def test_heteroaromatic_chain_sulfonic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_heteroaromatic_direct_attachment_sulfonic_acid_raises():
    # Unlike benzene, a heteroaromatic ring's numbering must fix the
    # heteroatom at locant 1 and search for the -SO3H's own lowest locant
    # relative to it -- ring-locant-search machinery out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)c1cccnc1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_carboxylic_acid.py`/`_aldehyde.py`/`_ketone.py` (CIP
        # computed entirely by RDKit's `rdCIPLabeler`). PubChem CID
        # 46398806.
        ("CC[C@@H](C)S(=O)(=O)O", "(2R)-butane-2-sulfonic acid"),
        ("CC[C@H](C)S(=O)(=O)O", "(2S)-butane-2-sulfonic acid"),
    ],
)
def test_acyclic_sulfonic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_sulfonic_acid_stereocenter_with_coexisting_substituent():
    # A stereocenter that also bears a halogen substituent: the suffix's
    # own locant is still cited on the 2-carbon chain despite the
    # substituent sharing its position (this project's own established
    # convention, unaffected by this PR -- see the identical pre-existing
    # 'CC(Cl)S(=O)(=O)O' -> '1-chloroethane-1-sulfonic acid'). PubChem's
    # own auto-generated name omits that locant ('(1R)-1-chloroethanesulfonic
    # acid', CID 124389840) -- a non-PIN quirk already documented elsewhere
    # in this project -- so only the structure is cross-checked there.
    assert smiles_to_iupac("C[C@@H](Cl)S(=O)(=O)O") == "(1R)-1-chloroethane-1-sulfonic acid"


def test_cyclic_sulfonic_acid_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_ketone.py`'s `_name_cyclic_ketone`. PubChem CID 101028852 (name
    # carries a redundant 'trans-' relative descriptor this project drops
    # once full R/S is given, same policy as the existing ring-alcohol
    # stereocenter task).
    assert (
        smiles_to_iupac("O=S(=O)(O)[C@H]1CCCC[C@@H]1Cl")
        == "(1S,2S)-2-chlorocyclohexane-1-sulfonic acid"
    )


def test_cyclic_sulfonic_acid_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_alcohol.py`'s
    # `_name_cyclic_alcohol`/`_ketone.py`/`_thiol.py` -- the branch and the
    # -SO3H share the same ring carbon (C1), same shape as the alcohol/
    # thiol precedent. PubChem has no cached record for this exact
    # structure (CID 0, sparse data gap), but the unstereo parent
    # ('1-ethylcyclohexane-1-sulfonic acid') matches PubChem exactly, and
    # the descriptor-construction code itself is identical to the
    # already-verified alcohol/ketone/thiol cases.
    assert (
        smiles_to_iupac("OS(=O)(=O)C1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-sulfonic acid"
    )


def test_sulfonic_acid_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(Cl)S(=O)(=O)O") == "1-chloropropane-1-sulfonic acid"


def test_sulfonic_acid_partially_specified_stereocenters_raises():
    # Only one of the ring's two genuine stereocenters is marked -- the
    # bug this PR fixes used to silently drop the marker and emit an
    # incomplete/wrong name; it must now raise instead.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(O)[C@H]1CCCCC1Cl")
