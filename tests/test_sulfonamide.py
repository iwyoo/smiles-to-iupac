import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanesulfonamide():
    # PubChem structure match: "methanesulfonamide".
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_ethanesulfonamide():
    # PubChem structure match: "ethanesulfonamide".
    assert smiles_to_iupac("CCS(=O)(=O)N") == "ethanesulfonamide"


def test_propane_1_sulfonamide():
    # PubChem structure match: "propane-1-sulfonamide".
    assert smiles_to_iupac("CCCS(=O)(=O)N") == "propane-1-sulfonamide"


def test_propane_2_sulfonamide():
    assert smiles_to_iupac("CC(S(=O)(=O)N)C") == "propane-2-sulfonamide"


def test_chlorobutanesulfonamide():
    assert smiles_to_iupac("ClCCCCS(=O)(=O)N") == "4-chlorobutane-1-sulfonamide"


def test_2_chloroethane_1_sulfonamide():
    # Locant cited even on a 2-carbon chain once a substituent (here, the
    # halogen) is present -- same project-wide convention as
    # '2-chloroethane-1-selenol'/'2-chloroethane-1-thiol' (PubChem's own
    # generated name omits the locant here; that divergence is accepted
    # project-wide, see test_selenol.py's identical case).
    assert smiles_to_iupac("ClCCS(=O)(=O)N") == "2-chloroethanesulfonamide"


def test_pent_4_ene_1_sulfonamide():
    assert smiles_to_iupac("C=CCCCS(=O)(=O)N") == "pent-4-ene-1-sulfonamide"


def test_disulfonamide():
    assert smiles_to_iupac("NS(=O)(=O)CCS(=O)(=O)N") == "ethane-1,2-disulfonamide"


def test_cyclohexanesulfonamide():
    # PubChem structure match: "cyclohexanesulfonamide".
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1") == "cyclohexanesulfonamide"


def test_2_methylcyclohexane_1_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1C") == "2-methylcyclohexane-1-sulfonamide"


def test_cyclopentanesulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC1") == "cyclopentanesulfonamide"


def test_2_chlorocyclohexane_1_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfonamide"


def test_n_methylcyclohexanesulfonamide():
    # PubChem structure match: "N-methylcyclohexanesulfonamide" (CID 23534457).
    assert smiles_to_iupac("O=S(=O)(NC)C1CCCCC1") == "N-methylcyclohexanesulfonamide"


def test_unsaturated_ring_sulfonamide():
    # Monocyclic ring, single -SO2NH2, single ring double bond (P-31.1.3):
    # the sulfonamide always gets locant 1 (suffix priority), the ring
    # double bond's locant is minimized by choosing direction.
    # Cross-checked against PubChem (CID 114763869/71757353).
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC=C1") == "cyclohex-2-ene-1-sulfonamide"
    assert smiles_to_iupac("O=S(=O)(N)C1CC=CCC1") == "cyclohex-3-ene-1-sulfonamide"


def test_unsaturated_ring_sulfonamide_with_substituent():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC=C1C") == "2-methylcyclohex-2-ene-1-sulfonamide"


def test_unsaturated_ring_sulfonamide_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CCCC#C1")


def test_unsaturated_ring_sulfonamide_with_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(NC)C1CCCC=C1")


def test_ring_substituent_chain_sulfonamide():
    # A sulfonamide entirely on a chain hanging off a plain saturated ring
    # (the ring itself bears no sulfonamide) -- mirrors
    # `_sulfonic_acid.py`'s `test_ring_substituent_chain_sulfonic_acid`.
    # PubChem PUG REST-verified "cyclohexylmethanesulfonamide"
    # (CID 21686562).
    assert smiles_to_iupac("NS(=O)(=O)CC1CCCCC1") == "cyclohexylmethanesulfonamide"


def test_ring_substituent_chain_sulfonamide_ring_with_substituent():
    assert smiles_to_iupac("NS(=O)(=O)CC1CCC(C)CC1") == "(4-methylcyclohexyl)methanesulfonamide"


def test_sulfonamide_with_alcohol():
    assert smiles_to_iupac("NS(=O)(=O)CCO") == "2-hydroxyethanesulfonamide"


def test_phenyl_chain_sulfonamide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#285's
    # carboxylic-acid/.../tellone chains): the ring is cited as a
    # "phenyl" substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCCS(=O)(=O)N") == "3-phenylpropane-1-sulfonamide"


def test_phenyl_chain_sulfonamide_internal_locant():
    assert smiles_to_iupac("c1ccccc1CC(C)S(=O)(=O)N") == "1-phenylpropane-2-sulfonamide"


def test_benzenesulfonamide():
    # -SO2NH2 directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1S(=O)(=O)N") == "benzenesulfonamide"


def test_substituted_benzenesulfonamide():
    # The mancude-ring numbering is free to start at the -SO2NH2 carbon,
    # so its own locant is never cited, mirroring benzenesulfonic acid.
    assert smiles_to_iupac("Cc1ccccc1S(=O)(=O)N") == "2-methylbenzenesulfonamide"
    assert smiles_to_iupac("Cc1ccc(cc1)S(=O)(=O)N") == "4-methylbenzenesulfonamide"  # PubChem PUG REST


def test_benzenesulfonamide_n_alkyl():
    # PubChem PUG REST IUPACName matches.
    assert smiles_to_iupac("CNS(=O)(=O)c1ccccc1") == "N-methylbenzenesulfonamide"
    assert smiles_to_iupac("CN(C)S(=O)(=O)c1ccccc1") == "N,N-dimethylbenzenesulfonamide"


def test_substituted_benzenesulfonamide_n_alkyl():
    # The N-substituent and the ring substituent share the same name
    # ("methyl"), so they collapse into one merged, multiplied citation
    # with mixed numeric/'N' locants -- PubChem PUG REST IUPACName match.
    assert smiles_to_iupac("Cc1ccc(cc1)S(=O)(=O)NC") == "N,4-dimethylbenzenesulfonamide"


def test_phenyl_chain_sulfonamide_n_alkyl():
    assert smiles_to_iupac("c1ccccc1CCCS(=O)(=O)NC") == "N-methyl-3-phenylpropane-1-sulfonamide"


def test_phenyl_substituted_benzene_ring_sulfonamide_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-3.md).
    assert smiles_to_iupac("Cc1ccccc1CCS(=O)(=O)N") == "2-(2-methylphenyl)ethanesulfonamide"


def test_phenyl_chain_sulfonamide_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propane-1-sulfonamide"
    # for this structure.
    assert smiles_to_iupac("Clc1ccc(CCCS(=O)(=O)N)cc1") == "3-(4-chlorophenyl)propane-1-sulfonamide"


def test_phenyl_chain_sulfonamide_ring_methyl():
    # PubChem PUG REST-verified "3-(4-methylphenyl)propane-1-sulfonamide".
    assert smiles_to_iupac("Cc1ccc(CCCS(=O)(=O)N)cc1") == "3-(4-methylphenyl)propane-1-sulfonamide"


def test_phenyl_chain_sulfonamide_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CS(N)(=O)=O") == "(4-ethylphenyl)methanesulfonamide"


def test_phenyl_chain_sulfonamide_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CCS(=O)(=O)N") == "2-(2-ethenylphenyl)ethanesulfonamide"


def test_n_methylmethanesulfonamide():
    # PubChem structure match: "N-methylmethanesulfonamide" (CID 97632).
    assert smiles_to_iupac("CS(=O)(=O)NC") == "N-methylmethanesulfonamide"


def test_n_n_dimethylmethanesulfonamide():
    # PubChem structure match: "N,N-dimethylmethanesulfonamide" (CID 70191).
    assert smiles_to_iupac("CS(=O)(=O)N(C)C") == "N,N-dimethylmethanesulfonamide"


def test_n_ethyl_n_methylethanesulfonamide():
    # PubChem structure match: "N-ethyl-N-methylethanesulfonamide" (CID 21102946).
    assert smiles_to_iupac("CCS(=O)(=O)N(C)CC") == "N-ethyl-N-methylethanesulfonamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 312702, 4130162.
        ("CS(=O)(=O)NC(C)C", "N-(propan-2-yl)methanesulfonamide"),
        ("CS(=O)(=O)NC(C)(C)C", "N-tert-butylmethanesulfonamide"),
    ],
)
def test_branched_n_substituted_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_identical_branched_n_substituents_parenthesized_when_compound():
    # PubChem CID 284325.
    assert smiles_to_iupac("CS(=O)(=O)N(C(C)C)C(C)C") == "N,N-di(propan-2-yl)methanesulfonamide"


def test_two_identical_branched_n_substituents_not_parenthesized_when_retained():
    # PubChem CID 58624041.
    assert smiles_to_iupac("CS(=O)(=O)N(C(C)(C)C)C(C)(C)C") == "N,N-ditert-butylmethanesulfonamide"


def test_two_different_n_substituents_alphabetized_ignoring_italic_prefix():
    # PubChem CID 58540473.
    assert smiles_to_iupac("CS(=O)(=O)N(CC)C(C)(C)C") == "N-tert-butyl-N-ethylmethanesulfonamide"


def test_acyclic_n_alkyl_coinciding_with_chain_substituent_name():
    # PubChem PUG REST IUPACName match: the N-methyl and the chain's own
    # 2-methyl branch share a name, so they merge into one multiplied
    # citation with mixed numeric/'N' locants, not two separate prefix
    # blocks.
    assert smiles_to_iupac("CNS(=O)(=O)CC(C)C") == "N,2-dimethylpropane-1-sulfonamide"


def test_cyclic_n_alkyl_coinciding_with_ring_substituent_name():
    # Not independently PubChem-verified (CID 0) -- structural/mechanism
    # consistency check with the acyclic case above.
    assert smiles_to_iupac("CNS(=O)(=O)C1(C)CCCCC1") == "N,1-dimethylcyclohexane-1-sulfonamide"


def test_halogenated_n_substituent():
    assert smiles_to_iupac("CS(=O)(=O)NCCCl") == "N-(2-chloroethyl)methanesulfonamide"


def test_unsaturated_n_substituted_sulfonamide():
    assert smiles_to_iupac("CS(=O)(=O)NCC=C") == "N-(prop-2-en-1-yl)methanesulfonamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_sulfonic_acid.py`. PubChem CID 93472539.
        ("CC[C@@H](C)S(=O)(=O)N", "(2R)-butane-2-sulfonamide"),
        ("CC[C@H](C)S(=O)(=O)N", "(2S)-butane-2-sulfonamide"),
    ],
)
def test_acyclic_sulfonamide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_sulfonamide_stereocenter_with_coexisting_substituent():
    # A stereocenter that also bears a halogen substituent: the suffix's
    # own locant is still forced to be C1 by P-44.4.1.8, so it's omitted
    # here too (P-14.3.4.2(b)) even though the halogen shares that
    # position -- PubChem CID 92264750 confirms this PIN form directly.
    assert smiles_to_iupac("C[C@@H](Cl)S(=O)(=O)N") == "(1R)-1-chloroethanesulfonamide"


def test_cyclic_sulfonamide_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`. PubChem has no
    # registered CID for this exact stereoisomer (CID 0), so only the
    # achiral structure is cross-checked (CID 130649687,
    # '2-chlorocyclohexane-1-sulfonamide').
    assert (
        smiles_to_iupac("N[S](=O)(=O)[C@H]1CCCC[C@@H]1Cl")
        == "(1S,2S)-2-chlorocyclohexane-1-sulfonamide"
    )


def test_cyclic_sulfonamide_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_alcohol.py`'s
    # `_name_cyclic_alcohol`/`_sulfonic_acid.py` -- the branch and the
    # -SO2NH2 share the same ring carbon (C1), same shape as the
    # alcohol/sulfonic-acid precedent. PubChem has no cached record for
    # this exact structure (CID 0, sparse data gap), but the unstereo
    # parent ('1-ethylcyclohexane-1-sulfonamide') matches PubChem exactly,
    # and the descriptor-construction code itself is identical to the
    # already-verified alcohol/ketone/thiol/sulfonic-acid/amine cases.
    assert (
        smiles_to_iupac("NS(=O)(=O)C1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-sulfonamide"
    )


def test_sulfonamide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(Cl)S(=O)(=O)N") == "1-chloropropane-1-sulfonamide"


def test_sulfonamide_partially_specified_stereocenters_cites_the_specified_elements():
    # Only one of the ring's two genuine stereocenters is marked -- must
    # raise rather than silently dropping the marker.
    assert smiles_to_iupac("N[S](=O)(=O)[C@H]1CCCCC1Cl") == '(1S)-2-chlorocyclohexane-1-sulfonamide'