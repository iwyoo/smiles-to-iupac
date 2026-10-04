import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). Note the 'al'
        # locant is never cited (see module docstring): the -CHO carbon is
        # always the chain terminus, so 'propan-1-al' is never written.
        ("CC=O", "ethanal"),
        ("CCC=O", "propanal"),
        ("CC(C)C=O", "2-methylpropanal"),
        # Dialdehyde: multiplying prefix, still no locants at all (both ends
        # of the chain, unambiguous), cross-checked against PubChem.
        ("O=CCCCC=O", "pentanedial"),
        # -al combined with existing unsaturation support: the 'ene' locant
        # is still cited (counted from the -CHO end, suffix priority,
        # P-44.4.1.8), cross-checked against PubChem.
        ("C=CCCC=O", "pent-4-enal"),
        # -al combined with a halogen substituent prefix, cross-checked
        # against PubChem.
        ("ClCCC=O", "3-chloropropanal"),
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/
        # #271/#272's carboxylic-acid/ketone/alcohol/ester chains): the
        # ring is cited as a "phenyl" substituent prefix. Cross-checked
        # against PubChem CID 7707 ("3-phenylpropanal").
        ("c1ccccc1CCC=O", "3-phenylpropanal"),
        # Two-carbon chain: PubChem's own name for this SMILES
        # ('2-phenylacetaldehyde', CID 998) uses the retained
        # 'acetaldehyde' stem, but this module always uses the systematic
        # 'ethanal' stem (see 'ethanal' above), so this is an accepted,
        # reviewed result rather than a PubChem-confirmed one -- same
        # policy as the carboxylic-acid module's '2-phenylethanoic acid'
        # (PR #269).
        ("c1ccccc1CC=O", "2-phenylethanal"),
    ],
)
def test_aldehyde_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formaldehyde_raises():
    # A carbonyl carbon with zero carbon neighbors (formaldehyde) is out of
    # scope for this module (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=O")


def test_ketene_raises():
    # A ketene's carbonyl carbon (C=C=O) is itself doubly bonded to its
    # carbon neighbor (a cumulated double bond, sp-hybridized) -- it is not
    # an aldehyde's -CHO and must not be silently misnamed as one.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCCCCCCCCCCC=C=O")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=C=O")


def test_carboxylic_acid_not_misread_as_aldehyde():
    # A carboxylic acid is routed to the dedicated carboxylic-acid module
    # (see test_carboxylic_acid.py) instead of falling through here.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_benzaldehyde():
    # -CHO attached directly to a benzene ring carbon (P-66.6.1.1.3). The
    # retained name 'benzaldehyde' stands for the whole ring+CHO system
    # (like 'phenol'/'aniline'), so the aldehyde's own ring locant is
    # never cited, only other substituents'. PubChem PUG REST:
    # "benzaldehyde"/"4-methylbenzaldehyde"/"2-chlorobenzaldehyde" -- all
    # exact matches.
    assert smiles_to_iupac("O=Cc1ccccc1") == "benzaldehyde"
    assert smiles_to_iupac("O=Cc1ccc(C)cc1") == "4-methylbenzaldehyde"
    assert smiles_to_iupac("O=Cc1ccccc1Cl") == "2-chlorobenzaldehyde"


def test_phenyl_substituted_benzene_ring_aldehyde_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout.md) -- PubChem PUG
    # REST-verified "2-(2-methylphenyl)acetaldehyde".
    assert smiles_to_iupac("Cc1ccccc1CC=O") == "2-(2-methylphenyl)ethanal"


def test_phenyl_chain_aldehyde_ring_halogen():
    # PubChem PUG REST computes "2-(4-chlorophenyl)acetaldehyde" for this
    # structure -- this project's own systematic-stem convention (see
    # test_carboxylic_acid.py's identical halogenated-ring cases) prefers
    # 'ethanal' once substituted.
    assert smiles_to_iupac("O=CCc1ccc(Cl)cc1") == "2-(4-chlorophenyl)ethanal"


def test_phenyl_chain_aldehyde_ring_dihalogen():
    assert smiles_to_iupac("Clc1cc(Cl)ccc1CC=O") == "2-(2,4-dichlorophenyl)ethanal"


def test_phenyl_chain_aldehyde_ring_methyl():
    # PubChem PUG REST-verified "2-(4-methylphenyl)acetaldehyde".
    assert smiles_to_iupac("Cc1ccc(CC=O)cc1") == "2-(4-methylphenyl)ethanal"


def test_phenyl_chain_aldehyde_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CC=O") == "2-(4-ethylphenyl)ethanal"


def test_phenyl_chain_aldehyde_with_hydroxyl():
    assert smiles_to_iupac("OCc1ccccc1CC=O") == "2-[2-(hydroxymethyl)phenyl]ethanal"


def test_phenyl_chain_aldehyde_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC=O") == "2-(2-ethenylphenyl)ethanal"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single aromatic (benzo or heteroaromatic monocycle) ring
        # substituent on a chain that also carries the aldehyde -- the
        # shared heteroaromatic-substituent naming engine (#620) now
        # reaches this shape once `_aldehyde.py`'s own single-ring gate
        # accepts a heteroaromatic ring alongside plain benzene (M2 step
        # 3, mirroring `_nitrile.py`'s #631). PubChem-confirmed:
        # `O=CCc1cccnc1` -> "2-pyridin-3-ylacetaldehyde" (CID 11073404;
        # this project's convention parenthesizes the compound
        # substituent, matching `_carboxylic_acid.py`'s
        # "2-(pyridin-3-yl)ethanoic acid").
        ("O=CCc1cccnc1", "2-(pyridin-3-yl)ethanal"),
        # PubChem-confirmed: `O=CCc1cccs1` -> "2-thiophen-2-ylacetaldehyde"
        # (CID 6430722).
        ("O=CCc1cccs1", "2-(thiophen-2-yl)ethanal"),
        # PubChem-confirmed: `O=CCc1cc[nH]c1` ->
        # "2-(1H-pyrrol-3-yl)acetaldehyde" (CID 21314538).
        ("O=CCc1cc[nH]c1", "2-(1H-pyrrol-3-yl)ethanal"),
        # Two separate simple monocycles joined by one direct bond, one a
        # plain benzo/heteroaromatic ring with no substituent of its own
        # -- the aldehyde-bearing ring's own dispatch (`_name_ring_
        # aldehyde`) is reused unchanged, with the aromatic ring cited as
        # a substituent via `name_branch`, the same generalization
        # #622/#624/#628/#631 made for `_ketone.py`/`_alcohol.py`/
        # `_thiol.py`/`_nitrile.py`. PubChem-confirmed:
        # `O=CC1CCCCC1c1ccccc1` -> "2-phenylcyclohexane-1-carbaldehyde"
        # (CID 15554880).
        ("O=CC1CCCCC1c1ccccc1", "2-phenylcyclohexane-1-carbaldehyde"),
        # PubChem-confirmed: `O=CC1CCCCCC1c1ccccc1` ->
        # "2-phenylcycloheptane-1-carbaldehyde" (CID 87968077).
        ("O=CC1CCCCCC1c1ccccc1", "2-phenylcycloheptane-1-carbaldehyde"),
        # The aromatic ring need not sit adjacent to the aldehyde -- ring
        # numbering still picks the lower-locant direction (P-14.5.2).
        ("O=CC1CCCC(c2ccccc2)C1", "3-phenylcyclohexane-1-carbaldehyde"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_aldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_aromatic_substituent_aldehyde_substituted_ring():
    assert smiles_to_iupac("O=CC1CCCC(c2ccc(C)cc2)C1") == "3-(4-methylphenyl)cyclohexane-1-carbaldehyde"


def test_two_ring_aromatic_substituent_aldehyde_three_rings():
    assert smiles_to_iupac("O=CC1CCCC(c2ccccc2)C1c3ccccc3") == "2,3-diphenylcyclohexane-1-carbaldehyde"


def test_two_ring_aromatic_substituent_aldehyde_chain_aldehyde():
    assert smiles_to_iupac("O=CCC1CCCCC1c1ccccc1") == "2-(2-phenylcyclohexyl)ethanal"


def test_heteroaromatic_ring_directly_attached_aldehyde():
    assert smiles_to_iupac("O=Cc1cccnc1") == "pyridine-3-carbaldehyde"


def test_alcohol_aldehyde_mix_names_hydroxy_prefix():
    # 'al' outranks 'ol' in Table 3.3, so a coexisting -OH is cited as the
    # 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("OCC=O") == "2-hydroxyethanal"


def test_aldehyde_alcohol_mix_on_longer_chain():
    assert smiles_to_iupac("OCCCC=O") == "4-hydroxybutanal"


def test_aldehyde_enol_mix():
    assert smiles_to_iupac("OC=CC=O") == "3-hydroxyprop-2-enal"


def test_ring_aldehyde():
    # -CHO attached directly to a saturated monocyclic ring carbon
    # (P-66.6.1.1.3, the 'carbaldehyde' suffix) -- e.g.
    # 'cyclohexanecarbaldehyde'. PubChem PUG REST:
    # "cyclohexanecarbaldehyde"/"4-methylcyclohexane-1-carbaldehyde" --
    # both exact matches.
    assert smiles_to_iupac("O=CC1CCCCC1") == "cyclohexanecarbaldehyde"
    assert smiles_to_iupac("O=CC1CCC(C)CC1") == "4-methylcyclohexane-1-carbaldehyde"
    # A substituent sharing the same ring atom as the -CHO (a quaternary
    # ring carbon), verified via PubChem PUG REST.
    assert smiles_to_iupac("O=CC1(C)CCCCC1") == "1-methylcyclohexane-1-carbaldehyde"


def test_ring_aldehyde_multiple_groups():
    assert smiles_to_iupac("O=CC1CCC(C=O)CC1") == "cyclohexane-1,4-dicarbaldehyde"


def test_aldehyde_carbon_off_the_longest_chain():
    # P-44.1.1: the principal chain must carry the -CHO even when a longer
    # chain exists elsewhere.
    assert smiles_to_iupac("CCCCC(C=O)CCCC") == "2-butylhexanal"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_carboxylic_acid.py` (CIP computed entirely by RDKit's
        # `rdCIPLabeler`, not reimplemented here). PubChem CID 76956407.
        ("C[C@@H](Cl)C=O", "(2R)-2-chloropropanal"),
        ("C[C@H](Cl)C=O", "(2S)-2-chloropropanal"),
        # Two specified stereocenters, ascending-locant group (P-91.3).
        # PubChem CID 92160220, name matches exactly.
        ("C[C@H](Cl)[C@H](Cl)C=O", "(2R,3S)-2,3-dichlorobutanal"),
    ],
)
def test_aldehyde_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aldehyde_stereocenter_with_hydroxy_coexistence():
    # Two stereocenters coexisting with a standalone hydroxyl (already
    # supported by this module as the 'hydroxy' prefix). PubChem doesn't
    # have this exact stereoisomer registered (CID 0 for every @/@@
    # combination tried), so only the non-stereo parent structure is
    # verified (PubChem CID 24973902, "2-chloro-3-hydroxybutanal");
    # the R/S computation itself is RDKit's `rdCIPLabeler`, already
    # verified elsewhere (`_alcohol.py`/`_carboxylic_acid.py`).
    assert smiles_to_iupac("C[C@H](O)[C@H](Cl)C=O") == "(2S,3S)-2-chloro-3-hydroxybutanal"


def test_aldehyde_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CC(Cl)C=O") == "2-chloropropanal"


def test_aldehyde_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("C[C@H](Cl)C(Cl)C=O") == '(3S)-2,3-dichlorobutanal'


def test_acyclic_aldehyde_specified_ez_double_bond():
    # Standalone specified C=C E/Z stereo, no stereocenter -- mirrors
    # `_ketone.py`'s identical PR #710 case. PubChem CID 447466 confirms
    # the structure (SMILES "C/C=C/C=O"), though PubChem's own
    # auto-generated name is "(E)-but-2-enal" (no locant); this project's
    # existing convention is to always cite the locant explicitly.
    assert smiles_to_iupac("O=C/C=C/C") == "(2E)-but-2-enal"


def test_acyclic_aldehyde_stereocenter_with_ez_double_bond_coexistence():
    # P-91.3: a specified tetrahedral stereocenter and a specified C=C
    # double-bond E/Z element on the same acyclic aldehyde chain are cited
    # together in one ascending-locant group, mirroring `_alcohol.py`'s/
    # `_ketone.py`'s identical case. PubChem CIDs 101710438/131246396
    # confirm both structures (SMILES "C/C=C/[C@@H](C)C=O" /
    # "C/C=C/[C@H](C)C=O"), though PubChem's own names group by type
    # ("(E,2R)-...") rather than by ascending locant -- the expected
    # locant-ascending form here is the primary source's own P-91.3 rule.
    assert smiles_to_iupac("O=C[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enal"
    assert smiles_to_iupac("O=C[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enal"


def test_cyclic_aldehyde_ring_stereocenter():
    # A specified stereocenter on the ring itself, alongside the ring's
    # sole -CHO substituent (P-92), same pattern as
    # `_carboxylic_acid.py`'s `_name_ring_carboxylic_acid`. PubChem has
    # no cached record for this exact stereoisomer, but the CIP labels
    # are RDKit's `rdCIPLabeler`, already verified elsewhere, and the
    # unstereo parent ('2-chlorocyclohexane-1-carbaldehyde') matches this
    # project's own existing output.
    assert (
        smiles_to_iupac("O=C[C@H]1CCCC[C@@H]1Cl")
        == "(1R,2S)-2-chlorocyclohexane-1-carbaldehyde"
    )


def test_aldehyde_ring_branch_stereocenter():
    assert smiles_to_iupac("O=C[C@H]1CC[C@H](C[C@@H](C)CC)C1") == "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carbaldehyde"


def test_aldehyde_ring_stereocenter_unspecified_unaffected():
    assert smiles_to_iupac("O=CC1CCCCC1Cl") == "2-chlorocyclohexane-1-carbaldehyde"


def test_aldehyde_ring_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("O=C[C@H]1CCCCC1Cl") == '(1R)-2-chlorocyclohexane-1-carbaldehyde'