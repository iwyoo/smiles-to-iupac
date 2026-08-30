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


def test_aryl_nitrile_raises():
    # Benzonitrile: an aromatic ring elsewhere in the molecule is out of
    # scope for this module (separate, in-progress aromatic-ring module's
    # territory).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#Cc1ccccc1")


def test_ring_nitrile_raises():
    # -C#N on a ring is the 'carbonitrile' suffix (P-66.5.1.2), a different
    # naming pattern this module deliberately excludes (see module
    # docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N#CC1CCCCC1")


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
