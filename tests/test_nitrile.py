import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). Note the 'nitrile'
        # locant is never cited (see module docstring): the -C#N carbon is
        # always the chain terminus, so 'propane-1-nitrile' is never written.
        ("CC#N", "ethanenitrile"),
        ("CCC#N", "propanenitrile"),
        ("CC(C)C#N", "2-methylpropanenitrile"),
        # -nitrile combined with existing unsaturation support: the 'ene'
        # locant is still cited (counted from the -C#N end, suffix priority,
        # P-44.4.1.8), cross-checked against PubChem. Note the doubled 'e'
        # (no elision before a consonant-initial suffix, see module
        # docstring).
        ("C=CCC#N", "but-3-enenitrile"),
        # -nitrile combined with a halogen substituent prefix, cross-checked
        # against PubChem.
        ("ClCCC#N", "3-chloropropanenitrile"),
    ],
)
def test_nitrile_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formonitrile_raises():
    # A nitrile carbon with zero carbon neighbors (bare HC#N) is out of
    # scope for this module (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C#N")


def test_dinitrile():
    # PubChem "butanedinitrile" (succinonitrile) -- both nitrile carbons
    # are chain termini, so the existing count-parameterized naming
    # helpers (`_suffix_body`/`_candidate_key`) already generalize; only
    # the validation gate needed loosening (see module docstring).
    assert smiles_to_iupac("N#CCCC#N") == "butanedinitrile"


def test_dinitrile_longer_chain():
    assert smiles_to_iupac("N#CCCCC#N") == "pentanedinitrile"


def test_dinitrile_branched_substituent():
    assert smiles_to_iupac("N#CCC(C)(C#N)C") == "2,2-dimethylbutanedinitrile"


def test_dinitrile_halogen_substituent():
    assert smiles_to_iupac("ClC(C#N)CC#N") == "2-chlorobutanedinitrile"


def test_trinitrile():
    assert smiles_to_iupac("N#CC(CC#N)CC#N") == "propane-1,2,3-tricarbonitrile"


def test_dinitrile_alongside_ring():
    assert smiles_to_iupac("N#CC1CCC(C#N)CC1") == "cyclohexane-1,4-dicarbonitrile"


def test_benzonitrile():
    # -C#N attached directly to a benzene ring carbon (P-66.5.1.1.3). The
    # retained name 'benzonitrile' stands for the whole ring+CN system
    # (like 'benzaldehyde'), so the nitrile's own ring locant is never
    # cited, only other substituents'. PubChem PUG REST:
    # "benzonitrile"/"4-methylbenzonitrile"/"2-chlorobenzonitrile" -- all
    # exact matches.
    assert smiles_to_iupac("N#Cc1ccccc1") == "benzonitrile"
    assert smiles_to_iupac("N#Cc1ccc(C)cc1") == "4-methylbenzonitrile"
    assert smiles_to_iupac("N#Cc1ccccc1Cl") == "2-chlorobenzonitrile"


def test_ring_nitrile():
    # -C#N attached directly to a saturated monocyclic ring carbon
    # (P-66.5.1.1.3, the 'carbonitrile' suffix) -- e.g.
    # 'cyclohexanecarbonitrile'. PubChem PUG REST:
    # "cyclohexanecarbonitrile"/"4-methylcyclohexane-1-carbonitrile" --
    # both exact matches.
    assert smiles_to_iupac("N#CC1CCCCC1") == "cyclohexanecarbonitrile"
    assert smiles_to_iupac("N#CC1CCC(C)CC1") == "4-methylcyclohexane-1-carbonitrile"
    # A substituent sharing the same ring atom as the -C#N (a quaternary
    # ring carbon), verified via PubChem PUG REST.
    assert smiles_to_iupac("N#CC1(C)CCCCC1") == "1-methylcyclohexane-1-carbonitrile"


def test_nitrile_carbon_off_the_longest_chain():
    assert smiles_to_iupac("CCCCC(C#N)CCCC") == "2-butylhexanenitrile"


def test_amine_still_routes_to_amine_module():
    # A primary amine with no nitrile present must still reach
    # `name_amine` unchanged (see core.py's nitrogen-branch routing).
    assert smiles_to_iupac("CCN") == "ethanamine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/
        # #271/#272/#273/#274's carboxylic-acid/ketone/alcohol/ester/
        # aldehyde/amide chains): the ring is cited as a "phenyl"
        # substituent prefix. Cross-checked against PubChem CID 22795
        # ("3-phenylpropanenitrile").
        ("c1ccccc1CCC#N", "3-phenylpropanenitrile"),
        # Two-carbon chain: PubChem's own name for this SMILES uses the
        # retained "phenylacetonitrile" stem, but this module always uses
        # the systematic 'ethanenitrile' stem (see test_nitrile_names
        # above), so this is an accepted, reviewed result rather than a
        # PubChem-confirmed one -- same policy as the aldehyde module's
        # '2-phenylethanal' (PR #273).
        ("c1ccccc1CC#N", "2-phenylethanenitrile"),
    ],
)
def test_phenyl_chain_nitrile_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_nitrile():
    # Same benzonitrile-type case as test_benzonitrile, reached via a
    # different ring-atom ordering in the SMILES.
    assert smiles_to_iupac("c1ccccc1C#N") == "benzonitrile"


def test_phenyl_substituted_benzene_ring_nitrile_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout.md) -- PubChem PUG
    # REST-verified "2-(2-methylphenyl)acetonitrile".
    assert smiles_to_iupac("Cc1ccccc1CC#N") == "2-(2-methylphenyl)ethanenitrile"


def test_phenyl_chain_nitrile_ring_halogen():
    # PubChem PUG REST computes "2-(4-chlorophenyl)acetonitrile" for this
    # structure -- this project's own systematic-stem convention (see
    # test_carboxylic_acid.py's identical halogenated-ring cases) prefers
    # 'ethanenitrile' once substituted.
    assert smiles_to_iupac("N#CCc1ccc(Cl)cc1") == "2-(4-chlorophenyl)ethanenitrile"


def test_phenyl_chain_nitrile_ring_dihalogen():
    assert smiles_to_iupac("Clc1cc(Cl)ccc1CC#N") == "2-(2,4-dichlorophenyl)ethanenitrile"


def test_phenyl_chain_nitrile_ring_methyl():
    # PubChem PUG REST-verified "2-(4-methylphenyl)acetonitrile".
    assert smiles_to_iupac("Cc1ccc(CC#N)cc1") == "2-(4-methylphenyl)ethanenitrile"


def test_phenyl_chain_nitrile_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CC#N") == "2-(4-ethylphenyl)ethanenitrile"


def test_phenyl_chain_nitrile_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC#N") == "2-(2-ethenylphenyl)ethanenitrile"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single aromatic (benzo or heteroaromatic monocycle) ring
        # substituent on a chain that also carries the nitrile -- the
        # shared heteroaromatic-substituent naming engine (#620) now
        # reaches this shape once `_nitrile.py`'s own single-ring gate
        # accepts a heteroaromatic ring alongside plain benzene (M2 step
        # 2, mirroring `_thiol.py`'s #628). PubChem-confirmed:
        # `N#CCc1cccnc1` -> "2-pyridin-3-ylacetonitrile" (CID 80923; this
        # project's convention parenthesizes the compound substituent,
        # matching `_carboxylic_acid.py`'s "2-(pyridin-3-yl)ethanoic
        # acid").
        ("N#CCc1cccnc1", "2-(pyridin-3-yl)ethanenitrile"),
        # PubChem-confirmed: `N#CCc1cccs1` -> "2-thiophen-2-ylacetonitrile"
        # (CID 72880).
        ("N#CCc1cccs1", "2-(thiophen-2-yl)ethanenitrile"),
        # PubChem-confirmed: `N#CCc1cc[nH]c1` ->
        # "2-(1H-pyrrol-3-yl)acetonitrile" (CID 23138077).
        ("N#CCc1cc[nH]c1", "2-(1H-pyrrol-3-yl)ethanenitrile"),
        # Two separate simple monocycles joined by one direct bond, one a
        # plain benzo/heteroaromatic ring with no substituent of its own
        # -- the nitrile-bearing ring's own dispatch (`_name_ring_
        # nitrile`) is reused unchanged, with the aromatic ring cited as a
        # substituent via `name_branch`, the same generalization #622/
        # #624/#628 made for `_ketone.py`/`_alcohol.py`/`_thiol.py`.
        # PubChem-confirmed: `N#CC1CCCCC1c1ccccc1` ->
        # "2-phenylcyclohexane-1-carbonitrile" (CID 13645775).
        ("N#CC1CCCCC1c1ccccc1", "2-phenylcyclohexane-1-carbonitrile"),
        # PubChem-confirmed: `N#CC1CCCCC1c1cccnc1` ->
        # "2-pyridin-3-ylcyclohexane-1-carbonitrile" (CID 82241004).
        ("N#CC1CCCCC1c1cccnc1", "2-(pyridin-3-yl)cyclohexane-1-carbonitrile"),
        ("N#CC1CCCCC1c1cccs1", "2-(thiophen-2-yl)cyclohexane-1-carbonitrile"),
        ("N#CC1CCCCC1c1ccc[nH]1", "2-(1H-pyrrol-2-yl)cyclohexane-1-carbonitrile"),
        # PubChem-confirmed: `N#CC1CCCCCC1c1ccccc1` ->
        # "2-phenylcycloheptane-1-carbonitrile" (CID 82241213).
        ("N#CC1CCCCCC1c1ccccc1", "2-phenylcycloheptane-1-carbonitrile"),
        # The aromatic ring need not sit adjacent to the nitrile -- ring
        # numbering still picks the lower-locant direction (P-14.5.2).
        ("N#CC1CCCC(c2ccccc2)C1", "3-phenylcyclohexane-1-carbonitrile"),
    ],
)
def test_two_ring_and_heteroaromatic_chain_nitrile(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_ring_aromatic_substituent_nitrile_substituted_ring():
    assert smiles_to_iupac("N#CC1CCCC(c2ccc(C)cc2)C1") == "3-(4-methylphenyl)cyclohexane-1-carbonitrile"


def test_two_ring_aromatic_substituent_nitrile_three_rings():
    assert smiles_to_iupac("N#CC1CCCC(c2ccccc2)C1c3ccccc3") == "2,3-diphenylcyclohexane-1-carbonitrile"


def test_two_ring_aromatic_substituent_nitrile_chain_nitrile():
    assert smiles_to_iupac("N#CCC1CCCCC1c1ccccc1") == "2-(2-phenylcyclohexyl)ethanenitrile"


def test_heteroaromatic_ring_directly_attached_nitrile():
    assert smiles_to_iupac("N#Cc1cccnc1") == "pyridine-3-carbonitrile"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92) -- the nitrile
        # carbon itself (sp, triple-bonded to nitrogen) is never a
        # potential stereocenter, confirmed via RDKit `FindPotentialStereo`.
        # PubChem CID 5479121.
        ("CC[C@@H](C)C#N", "(2R)-2-methylbutanenitrile"),
        ("CC[C@H](C)C#N", "(2S)-2-methylbutanenitrile"),
    ],
)
def test_nitrile_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nitrile_stereocenter_with_coexisting_halogen():
    # PubChem CID 88154506.
    assert smiles_to_iupac("CC[C@@H](Cl)C#N") == "(2R)-2-chlorobutanenitrile"


def test_nitrile_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)C#N") == "2-methylbutanenitrile"


def test_acyclic_nitrile_specified_ez_double_bond():
    # Standalone specified C=C E/Z stereo, no stereocenter -- mirrors
    # `_ketone.py`/`_aldehyde.py`/`_carboxylic_acid.py`'s identical
    # PR #710/#713/#715 case. PubChem CID 637921 confirms the structure
    # (SMILES "C/C=C/C#N"), though PubChem's own auto-generated name is
    # "(E)-but-2-enenitrile" (no locant); this project's existing
    # convention is to always cite the locant explicitly.
    assert smiles_to_iupac("N#C/C=C/C") == "(2E)-but-2-enenitrile"


def test_acyclic_nitrile_stereocenter_with_ez_double_bond_coexistence():
    # P-91.3: a specified tetrahedral stereocenter and a specified C=C
    # double-bond E/Z element on the same acyclic nitrile chain are cited
    # together in one ascending-locant group, mirroring `_alcohol.py`'s/
    # `_ketone.py`'s/`_aldehyde.py`'s/`_carboxylic_acid.py`'s identical
    # case. PubChem CIDs 92976641/92976642 confirm both structures
    # (SMILES "C/C=C/[C@@H](C)C#N" / "C/C=C/[C@H](C)C#N"), though
    # PubChem's own names group by type ("(E,2R)-...") rather than by
    # ascending locant -- the expected locant-ascending form here is the
    # primary source's own P-91.3 rule.
    assert smiles_to_iupac("N#C[C@H](C)/C=C/C") == "(2R,3E)-2-methylpent-3-enenitrile"
    assert smiles_to_iupac("N#C[C@@H](C)/C=C/C") == "(2S,3E)-2-methylpent-3-enenitrile"


def test_cyclic_nitrile_ring_stereocenter():
    # A specified stereocenter on the ring itself, alongside the ring's
    # sole -C#N substituent (P-92), same pattern as
    # `_carboxylic_acid.py`'s `_name_ring_carboxylic_acid`.
    assert (
        smiles_to_iupac("N#C[C@H]1CCCC[C@@H]1Cl")
        == "(1R,2S)-2-chlorocyclohexane-1-carbonitrile"
    )


def test_nitrile_ring_branch_stereocenter():
    assert smiles_to_iupac("N#C[C@H]1CC[C@H](C[C@@H](C)CC)C1") == "(1S,3R)-3-[(2S)-2-methylbutyl]cyclopentane-1-carbonitrile"


def test_nitrile_ring_stereocenter_unspecified_unaffected():
    assert smiles_to_iupac("N#CC1CCCCC1Cl") == "2-chlorocyclohexane-1-carbonitrile"


def test_nitrile_ring_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("N#C[C@H]1CCCCC1Cl") == '(1R)-2-chlorocyclohexane-1-carbonitrile'