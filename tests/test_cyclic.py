import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC1", "cyclopropane"),
        ("C1CCC1", "cyclobutane"),
        ("C1CCCCC1", "cyclohexane"),
        ("CC1CCCCC1", "methylcyclohexane"),
        ("CC1CCCCC1C", "1,2-dimethylcyclohexane"),
        ("CC1CC(C)CCC1", "1,3-dimethylcyclohexane"),
    ],
)
def test_smiles_to_iupac_cyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polycyclic_raises():
    # pentagonal prism: a hexacyclic ring system (cyclomatic number 6, ten
    # branch atoms of degree 3, two pentagons joined by five bridging bonds):
    # out of scope for bicyclic through pentacyclic support alike (see
    # _bicyclic.py, _polycyclic.py). Built from scratch via RDKit's RWMol
    # (two independent 5-cycles plus one bond between each corresponding
    # pair of atoms), not copied from a database.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C1C2C2C3C4C12")


def test_ring_compound_substituent():
    # sec-butyl-like branch on the ring (P-29.4); the only substituent on an
    # otherwise unsubstituted ring, so its locant is omitted (P-14.3.3).
    assert smiles_to_iupac("CC(CC)C1CCCCC1") == "(1-methylpropyl)cyclohexane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # cis-/trans-1,2-dimethylcyclohexane (P-93.5.1.3 general
        # nomenclature, tasks/ring-cis-trans-naming.md, 2026-08-25):
        # cross-checked against PubChem's own isomeric SMILES for CID
        # 16628 (cis) and CID 23313 (trans).
        ("C[C@@H]1CCCC[C@@H]1C", "cis-1,2-dimethylcyclohexane"),
        ("C[C@@H]1CCCC[C@H]1C", "trans-1,2-dimethylcyclohexane"),
    ],
)
def test_smiles_to_iupac_ring_cis_trans(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_ring_stereo_unaffected():
    # no `@`/`@@` at all: not a new rejection case -- named exactly as
    # before (no cis-/trans- prefix), matching this project's prior
    # behavior for every other ring module.
    assert smiles_to_iupac("CC1CCCCC1C") == "1,2-dimethylcyclohexane"


def test_non_stereogenic_ring_position_unaffected():
    # a single substituted ring carbon on an otherwise plain ring has a
    # constitutional mirror plane through it (the two ring paths away from
    # it are equivalent), so it's never a real stereocenter regardless of
    # any `@`/`@@` marker -- not a new rejection case either.
    assert smiles_to_iupac("C[C@H]1CCCCC1") == "methylcyclohexane"


def test_asymmetric_1_2_disubstituted_stereo_raises():
    # methyl and ethyl (different substituents) on adjacent ring carbons,
    # both specified: P-93.5.1.3's full reference-substituent selection
    # rules would be needed, out of this module's minimal scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@@H]1CCCC[C@H]1CC")


def test_1_3_disubstituted_stereo_raises():
    # same symmetric substituents, but not adjacent (1,3- instead of
    # 1,2-) -- the geometric-plane-check logic isn't verified for this
    # position pattern, out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@@H]1C[C@H](C)CCC1")


def test_single_ring_stereocenter_raises():
    # a real stereocenter (four different groups) specified alone, with no
    # second one to compare against -- not the cis/trans shape this module
    # handles.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H]1CCCCC1CC")
