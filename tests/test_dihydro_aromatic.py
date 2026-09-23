import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,2-dihydronaphthalene: CAS 447-53-0, NIST WebBook and
        # Sigma-Aldrich both list "1,2-Dihydronaphthalene" as the name.
        ("C1CC=Cc2ccccc12", "1,2-dihydronaphthalene"),
        # 1,4-dihydronaphthalene: CAS 612-17-9, a distinct real compound
        # from the 1,2- isomer (NIST WebBook separate entry).
        ("C1C=CCc2ccccc12", "1,4-dihydronaphthalene"),
        # Same 1,4-dihydronaphthalene molecule as above, written starting
        # from the intact ring instead of the reduced ring -- the
        # lowest-locants tie-break (P-31.1.4.2) must still normalize this
        # to "1,4-", not "2,3-"/"5,8-"/"6,7-" (whichever the raw SMILES
        # atom order would otherwise suggest).
        ("c1ccc2c(c1)CC=CC2", "1,4-dihydronaphthalene"),
    ],
)
def test_smiles_to_iupac_dihydronaphthalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_naphthalene_itself_is_unaffected():
    # exactly zero hydro pairs still goes through _aromatic.py's own path,
    # not this module (find_dihydronaphthalene_core requires the reduced
    # ring to have exactly two sp3 atoms).
    assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


def test_tetrahydronaphthalene_raises():
    # more than one hydro pair (tetralin) is out of scope for this module
    # (see docstring) -- must not be mistaken for the single-pair case.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCc2ccccc12")


def test_substituted_dihydronaphthalene_raises():
    # any substituent is out of scope -- every ring atom must have exactly
    # its "bare" degree.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC=Cc2ccccc12")


def test_halogen_substituted_dihydronaphthalene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC1CC=Cc2ccccc12")


def test_decahydronaphthalene():
    # decahydronaphthalene (decalin): PubChem CID 7044 gives
    # "1,2,3,4,4a,5,6,7,8,8a-decahydronaphthalene" -- this project follows
    # the Blue Book's own plain worked example instead (P-31.2.3.3.2:
    # total hydrogenation locants are redundant and omitted, P-14.3.4.5).
    assert smiles_to_iupac("C1CCC2CCCCC2C1") == "decahydronaphthalene"


def test_decahydronaphthalene_not_von_baeyer():
    # Same molecule as above, written starting from a different atom --
    # must not be mistaken for a von Baeyer bicyclic (see #814).
    assert smiles_to_iupac("C1CCCC2CCCCC12") == "decahydronaphthalene"


def test_bridged_bicyclic_is_not_decahydronaphthalene():
    # bicyclo[2.2.2]octane has the same atom/bond count shape narrowed
    # down by num_rings==2, but a nonzero third bridge -- must still be
    # named via _bicyclic.py, not mistaken for naphthalene's skeleton.
    assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Octahydro (1 double bond remaining, the fusion pair 7-8): both
        # PubChem-confirmed (CID 565677's own structure, cited directly
        # in #828's own scoping).
        ("C1CC2CCCCC2C=C1", "1,2,3,4,4a,5,6,8a-octahydronaphthalene"),
        # Hexahydro (2 double bonds remaining): PubChem CID 561877.
        ("C1=CC2CCCCC2C=C1", "1,2,3,4,4a,8a-hexahydronaphthalene"),
        # Hexahydro, a different double-bond arrangement (5-6 and 7-8
        # kept, 1-2/3-4/4a-8a saturated): PubChem CID 22035261.
        ("C1CC2C=CCCC2C=C1", "1,2,4a,5,6,8a-hexahydronaphthalene"),
        # Hexahydro, yet another arrangement (1-2 and 5-6 kept):
        # PubChem CID 14120284.
        ("C1CC2CCC=CC2C=C1", "1,2,4a,7,8,8a-hexahydronaphthalene"),
        # Tetrahydro (3 double bonds remaining): PubChem CID 13056319.
        ("C1=CC2C=CCCC2C=C1", "1,2,4a,8a-tetrahydronaphthalene"),
        # Dihydro via the fusion-bond-only pair (4a,8a both saturated,
        # all four "peripheral" double bonds kept) -- a real, distinct
        # compound from `test_smiles_to_iupac_dihydronaphthalene`'s
        # 1,2-/1,4-dihydronaphthalene above (those keep the fusion bond
        # double and saturate a peripheral pair instead; this module's
        # own aromatic-ring requirement means neither ring here is
        # RDKit-aromatic-flagged, so `find_dihydronaphthalene_core`
        # correctly leaves this one to this newer, more general
        # detector). PubChem CID 15555390.
        ("C1=CC2C=CC=CC2C=C1", "4a,8a-dihydronaphthalene"),
    ],
)
def test_partially_unsaturated_naphthalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_arbitrary_double_bond_placement_not_matched():
    # A double-bond arrangement that isn't reachable by removing bonds
    # from naphthalene's own reference Kekule structure (#828's own scope
    # note) -- falls through to the pre-existing (unrelated to this
    # module) von Baeyer path rather than being misnamed as a hydro
    # derivative.
    assert smiles_to_iupac("C1CC2CC=CCC2C=C1") == "bicyclo[4.4.0]deca-2,8-diene"


def test_partially_unsaturated_naphthalene_substituent_falls_through():
    # A substituent means `find_partially_unsaturated_naphthalene_core`'s
    # own 10-atom bare-skeleton check returns None (mirrors the existing
    # dihydro/decahydro detectors' identical restriction), so this still
    # falls through to the pre-existing (unrelated to this task) von
    # Baeyer path rather than raising.
    assert smiles_to_iupac("CC1CC2CCCCC2C=C1") == "4-methylbicyclo[4.4.0]dec-2-ene"
