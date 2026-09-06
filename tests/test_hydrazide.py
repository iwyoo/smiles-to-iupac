import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName) and matching the
        # Blue Book's own worked examples (P6a.pdf, P-66.3.1.1/.1.2.3):
        # 'butanehydrazide (PIN)', 'pentanehydrazide (PIN)'. The
        # hydrazide carbon is always C1, and its own locant is never
        # cited (P-14.3.3), same as `_amide.py`.
        ("CCCC(=O)NN", "butanehydrazide"),
        ("CCCCC(=O)NN", "pentanehydrazide"),
        ("CCC(=O)NN", "propanehydrazide"),
        # -hydrazide + halogen substituent prefix.
        ("ClCCC(=O)NN", "3-chloropropanehydrazide"),
        # -hydrazide combined with existing unsaturation support.
        ("CC=CC(=O)NN", "but-2-enehydrazide"),
        # 'hydrazide' outranks 'ol' in Table 3.3, so a coexisting
        # standalone -OH is cited as the 'hydroxy' prefix, same as
        # `_amide.py`.
        ("OCCC(=O)NN", "3-hydroxypropanehydrazide"),
    ],
)
def test_hydrazide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formohydrazide():
    # P-66.3.1.2.1: 'formohydrazide' (a retained name), not the
    # systematic 'methanehydrazide', is the actual PIN here. PubChem CID
    # 12229.
    assert smiles_to_iupac("C(=O)NN") == "formohydrazide"


def test_acetohydrazide():
    # 'acetohydrazide' is the PIN, not 'ethanehydrazide'. PubChem CID
    # 14039.
    assert smiles_to_iupac("CC(=O)NN") == "acetohydrazide"


def test_2_chloroacetohydrazide():
    # A substituent on acetohydrazide's terminal carbon is cited as an
    # ordinary prefix. PubChem CID 101883.
    assert smiles_to_iupac("ClCC(=O)NN") == "2-chloroacetohydrazide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # N-/N'-substituted hydrazide (P-66.3.1.1): the carbonyl-adjacent
        # nitrogen is "N", the terminal one "N'", mirroring `_amide.py`'s
        # "N-"/`_urea.py`'s "N-"/"N'-" convention. All cross-checked
        # against PubChem PUG REST (compound/smiles/<smiles>/property/
        # IUPACName).
        ("CC(=O)N(C)N", "N-methylacetohydrazide"),  # CID 19051
        ("CC(=O)NNC", "N'-methylacetohydrazide"),  # CID 122488
        ("CC(=O)N(C)NC", "N,N'-dimethylacetohydrazide"),  # CID 12596294
        ("CC(=O)NN(C)C", "N',N'-dimethylacetohydrazide"),  # CID 80385
        ("CC(=O)N(C)N(C)C", "N,N',N'-trimethylacetohydrazide"),  # CID 12362621
        ("CC(=O)N(CC)NC", "N-ethyl-N'-methylacetohydrazide"),  # CID 89037658
        ("CC(=O)N(C)NCC", "N'-ethyl-N-methylacetohydrazide"),  # CID 119096326
        # The systematic (non-retained-name) chain lengths behave the same
        # way as the retained 'aceto-' case above.
        ("CCCC(=O)N(C)N", "N-methylbutanehydrazide"),  # CID 53649946
        ("CCCC(=O)NNC", "N'-methylbutanehydrazide"),  # CID 88556637
    ],
)
def test_n_substituted_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_stereocenter_hydrazide():
    # Mechanical combination of two independently-verified axes: N-prefix
    # placement (see above) and stereodescriptor-outermost ordering
    # (`test_hydrazide_stereocenter` below) -- not itself registered on
    # PubChem, so treated as a reviewed/eyeballed result rather than an
    # independently-verified one (see the base '(2R)-2-methylbutanehydrazide'
    # case this extends, PubChem CID 30066157).
    assert smiles_to_iupac("CC[C@@H](C)C(=O)N(C)N") == "(2R)-N-methyl-2-methylbutanehydrazide"


def test_n_substituted_formohydrazide_raises():
    # Unlike the dinuclear (aceto-) and longer chain cases, the mononuclear
    # retained name itself doesn't survive N-substitution -- PubChem
    # switches to plain 'formamide' + amino-substituent naming instead
    # (`O=CN(C)N` -> "N-amino-N-methylformamide", CID 12962), a different
    # naming paradigm this module doesn't attempt.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CN(C)N")


def test_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(N)C(C)C")


def test_diacylhydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCC(=O)NNC(=O)C")


def test_amide_not_misnamed_as_hydrazide():
    # A plain primary amide (`_amide.py`) has only one nitrogen and must
    # not be routed here.
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): a hydrazide
        # carbon is always chain-terminal (fixed C1), same pattern as
        # `_amide.py` (CIP computed entirely by RDKit's `rdCIPLabeler`).
        # PubChem CID 30066157.
        ("CC[C@@H](C)C(=O)NN", "(2R)-2-methylbutanehydrazide"),
        ("CC[C@H](C)C(=O)NN", "(2S)-2-methylbutanehydrazide"),
    ],
)
def test_hydrazide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hydrazide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(=O)NN") == "2-methylbutanehydrazide"


def test_phenyl_chain_hydrazide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#288's
    # carboxylic-acid/.../hydroperoxide chains): the ring is cited as a
    # "phenyl" substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCC(=O)NN") == "3-phenylpropanehydrazide"


def test_phenyl_chain_hydrazide_retained_name():
    # Dinuclear case: the retained 'acetohydrazide' name is still the PIN
    # even with the phenyl substituent, mirroring the module's own
    # '2-chloroacetohydrazide' pattern (module docstring).
    assert smiles_to_iupac("c1ccccc1CC(=O)NN") == "2-phenylacetohydrazide"


def test_phenyl_directly_attached_hydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)NN")


def test_phenyl_chain_hydrazide_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)N(C)N")


def test_n_prime_phenylhydrazide():
    # PubChem PUG REST IUPACName match: a plain phenyl ring hanging
    # directly off the hydrazide's terminal (N') nitrogen -- distinct
    # from the chain-substituent case above (ring far from any
    # nitrogen).
    assert smiles_to_iupac("CC(=O)NNc1ccccc1") == "N'-phenylacetohydrazide"
    assert smiles_to_iupac("CCC(=O)NNc1ccccc1") == "N'-phenylpropanehydrazide"


def test_n_phenylhydrazide():
    # PubChem PUG REST IUPACName match: the phenyl ring on the
    # carbonyl-adjacent (N) nitrogen instead.
    assert smiles_to_iupac("CC(=O)N(c1ccccc1)N") == "N-phenylacetohydrazide"


def test_n_and_n_prime_substituted_phenylhydrazide():
    # PubChem PUG REST IUPACName match: an alkyl on N and a phenyl on N'
    # coexist.
    assert smiles_to_iupac("CC(=O)N(C)Nc1ccccc1") == "N-methyl-N'-phenylacetohydrazide"


def test_n_phenylhydrazide_with_second_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(c1ccccc1)Nc1ccccc1")


def test_phenyl_substituted_benzene_ring_hydrazide_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-5.md).
    assert smiles_to_iupac("Cc1ccccc1CC(=O)NN") == "2-(2-methylphenyl)acetohydrazide"


def test_phenyl_chain_hydrazide_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propanehydrazide" for
    # this structure.
    assert smiles_to_iupac("Clc1ccc(CCC(=O)NN)cc1") == "3-(4-chlorophenyl)propanehydrazide"


def test_phenyl_chain_hydrazide_ring_methyl():
    assert smiles_to_iupac("Cc1ccc(CCC(=O)NN)cc1") == "3-(4-methylphenyl)propanehydrazide"


def test_phenyl_chain_hydrazide_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CCC(=O)NN") == "3-(4-ethylphenyl)propanehydrazide"


def test_phenyl_chain_hydrazide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)NN")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Symmetric diacylhydrazide (P-66.3.3.3), PubChem PUG REST
        # cross-checked: 'N'-acetylacetohydrazide' (CID 72884, both sides
        # the dinuclear retained name) and 'N'-propanoylpropanehydrazide'
        # (CID 73715, both sides a systematic three-carbon chain).
        ("CC(=O)NNC(=O)C", "N'-acetylacetohydrazide"),
        ("CCC(=O)NNC(=O)CC", "N'-propanoylpropanehydrazide"),
    ],
)
def test_symmetric_diacyl_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diacyl_hydrazide_formyl():
    # Mononuclear retained name on both sides.
    assert smiles_to_iupac("O=CNNC=O") == "N'-formylformohydrazide"


def test_asymmetric_diacyl_hydrazide_raises():
    # The two acyl groups name differently -- deciding which is senior
    # needs Table 3.3/4.4 acid-seniority handling, out of scope here
    # (mirrors `_imide.py`'s identical symmetric-only restriction).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NNC(=O)CC")


def test_diacyl_hydrazide_aromatic_acyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NNC(=O)c1ccccc1")


def test_diacyl_hydrazide_branched_acyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)C(=O)NNC(=O)C(C)C")


def test_halogen_substituted_n_alkyl_raises():
    # Same bug class as `_amidine.py`'s identical `_collect_n_alkyl`:
    # a halogen hanging off an N-substituent is invisible to the
    # carbon-only chain walk and passed through silently -- found via
    # real-data testing: 'CCC(=O)NNCC(F)(F)F' was misnamed
    # "N'-ethylpropanehydrazide", the -CH2CF3 substituent's three
    # fluorines vanishing entirely.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC(=O)NNCC(F)(F)F")
