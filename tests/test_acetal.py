import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName) and matching the
        # Blue Book's own worked example (P6a.pdf, P-66.6.5.1:
        # '1,1-diethoxypropane (PIN)').
        ("CCC(OCC)OCC", "1,1-diethoxypropane"),
        ("COC(OC)C", "1,1-dimethoxyethane"),
        # A ketal (neither R nor R' is hydrogen) -- named the same way.
        ("COC(C)(C)OC", "2,2-dimethoxypropane"),
        # A symmetric acetal carbon that is itself the parent chain's
        # midpoint (a real positional/multiplying-prefix check).
        ("CCC(OCC)(OCC)CC", "3,3-diethoxypentane"),
        # Two different alkoxy substituents, cited individually in
        # alphabetical order rather than combined with a multiplying
        # prefix.
        ("COC(OCC)C", "1-ethoxy-1-methoxyethane"),
    ],
)
def test_acetal_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_acetal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(OC)(OC)CC1")


def test_branched_alkoxy_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC(OC(C)C)OCC")


def test_unsaturated_acetal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(OCC)OCC")


def test_hemiacetal_raises():
    # RR'C(OH)(O-R'') (P-66.6.5.2) has only one alkoxy oxygen plus a
    # hydroxyl -- a different, unverified shape, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(OCC)CC")


def test_ether_not_misnamed_as_acetal():
    # A plain single ether (`_ether.py`) has only one oxygen and must not
    # be routed here.
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_ketone_not_misnamed_as_acetal():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"
