import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed (CID 10345, 107363, 139709, 276875, 137688,
        # 12550422 respectively).
        ("O=C1CC2CCC1C2", "bicyclo[2.2.1]heptan-2-one"),
        ("O=C1CCC2CCC1C2", "bicyclo[3.2.1]octan-2-one"),
        ("O=C1CC2CCC(C1)C2", "bicyclo[3.2.1]octan-3-one"),
        ("O=C1CC2CCCC(C1)C2", "bicyclo[3.3.1]nonan-3-one"),
        ("O=C1CC2CCC1CC2", "bicyclo[2.2.2]octan-2-one"),
        ("O=C1CC2CCC2C1", "bicyclo[3.2.0]heptan-3-one"),
        # 1-adamantanone -- this project's adamantane already uses the
        # systematic 'tricyclo[3.3.1.1^3,7]decane' name rather than the
        # retained 'adamantane' one (see `test_von_baeyer_alcohol.py`'s
        # identical 1-/2-adamantanol precedent), so this module's
        # ring_count>=3 path follows the same systematic convention
        # rather than PubChem's own retained-name 'adamantan-2-one'.
        ("O=C1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.1^3,7]decan-2-one"),
    ],
)
def test_von_baeyer_ketone_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_ketones_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2CCC1C(=O)C2")


def test_ketone_hydroxyl_combination_on_ring_system():
    # A coexisting hydroxyl on a von Baeyer bicyclic/polycyclic ketone is
    # supported (#1029) -- cited as a "hydroxy" prefix, same PubChem
    # structure as tests/test_ketone_hydroxyl_vonbaeyer.py's own
    # CID 85551307 case, written from a different starting atom here.
    assert smiles_to_iupac("O=C1CC2CCC1C2O") == "7-hydroxybicyclo[2.2.1]heptan-2-one"


def test_unsaturated_von_baeyer_ketone_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CC2C=CC1C2")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-confirmed real structures, each PubChem's own IUPACName
        # property matches exactly (#1078, M5 step 1): the bicyclic suffix
        # engine never learned to cite a specified stereocenter at all
        # before this, unconditionally rejecting every one of these.
        ("C[C@@]12CC[C@@H](C1(C)C)CC2=O", "(1R,4R)-1,7,7-trimethylbicyclo[2.2.1]heptan-2-one"),  # (+)-camphor, CID 159055
        ("C[C@]12CC[C@H](C1(C)C)CC2=O", "(1S,4S)-1,7,7-trimethylbicyclo[2.2.1]heptan-2-one"),  # (-)-camphor, CID 444294
        ("C[C@]12CC[C@H](C1)C(C2=O)(C)C", "(1S,4R)-1,3,3-trimethylbicyclo[2.2.1]heptan-2-one"),  # fenchone, CID 1201521
        ("C[C@@H]1[C@H]2C[C@]2(CC1=O)C(C)C", "(1S,4R,5R)-4-methyl-1-(propan-2-yl)bicyclo[3.1.0]hexan-3-one"),  # thujone, CID 261491
        ("C[C@@H]1[C@@H]2C[C@@H](C2(C)C)CC1=O", "(1S,2R,5R)-2,6,6-trimethylbicyclo[3.1.1]heptan-3-one"),  # pinocamphone, CID 6427105
    ],
)
def test_von_baeyer_ketone_specified_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_von_baeyer_ketone_stereocenter_with_coexisting_hydroxyl_still_raises():
    # A specified stereocenter alongside a *coexisting hydroxyl* is a
    # separate, still-unsupported combination (M5 step 2) -- this step
    # only wires stereo through the no-coexisting-hydroxyl path.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1C[C@H]2CC[C@H]1[C@@H]2O")


def test_von_baeyer_ketone_stereocenter_with_ring_unsaturation_still_raises():
    # Ring unsaturation alongside a stereocenter is separately unsupported
    # (M6) -- confirm the unsaturation guard still fires first.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1=CC(=O)[C@@H]2C[C@H]1C2(C)C")  # verbenone, CID 65724
