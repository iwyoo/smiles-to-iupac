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


def test_dinitrile_raises():
    # More than one nitrile group is future work (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#CCCC#N")


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


def test_bicyclic_carbon_skeleton_with_stray_nitrile_raises():
    # A -C#N group whose carbon is not on any single longest chain of the
    # molecule (here, a branch off a longer chain) is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC(C#N)CCCC")


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


def test_phenyl_substituted_benzene_ring_nitrile_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC#N")


def test_phenyl_chain_nitrile_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC#N")


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
