import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # cross-checked against PubChem PUG REST IUPACName:
        # CID 6331 (B), CID 5326078 (CB), CID 15879358 (CBC),
        # CID 68979 (CB(C)C), CID 7357 (CCB(CC)CC), CID 143352201 (CCBC),
        # CID 517955 (CCB(C)C), CID 543198 (CCB(CC)C).
        ("B", "borane"),
        ("CB", "methylborane"),
        ("CBC", "dimethylborane"),
        ("CB(C)C", "trimethylborane"),
        ("CCB(CC)CC", "triethylborane"),
        ("CCBC", "ethyl(methyl)borane"),
        ("CCB(C)C", "ethyl(dimethyl)borane"),
        ("CCB(CC)C", "diethyl(methyl)borane"),
        # halogen bonded directly to boron: cross-checked against PubChem
        # PUG REST IUPACName, CID 140714 (ClB), CID 137221 (ClB(C)C).
        ("ClB", "chloroborane"),
        ("ClB(C)C", "chloro(dimethyl)borane"),
    ],
)
def test_smiles_to_iupac_simple_borane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituted_alkyl_chain_raises():
    # a halogen embedded partway along a carbon chain (rather than bonded
    # directly to boron) is out of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClCCB")


def test_borane_chain_raises():
    # diborane (B-B bond): two boron atoms, a separate "chain" nomenclature
    # problem this module doesn't handle.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("BB")


def test_branched_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)B")


def test_unsaturated_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CB")


def test_aromatic_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1B")
