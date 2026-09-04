import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 123195: unsubstituted parent.
        ("N=N", "diazene"),
        # PubChem CID 123421: single substituent.
        ("CN=N", "methyldiazene"),
        # PubChem CID 10421: symmetric disubstitution.
        ("CN=NC", "dimethyldiazene"),
        # PubChem CID 526060: asymmetric disubstitution -- alphabetically
        # first ('ethyl') unparenthesized, the other ('methyl')
        # parenthesized (P-16.5.1.3.1, mirrors _phosphane.py's identical
        # rule, already verified there for 'ethyl(methyl)phosphane').
        ("CN=NCC", "ethyl(methyl)diazene"),
        # PubChem CID 13183.
        ("CCN=NCC", "diethyldiazene"),
    ],
)
def test_diazene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 22166172, 300540.
        ("CC(C)N=N", "propan-2-yldiazene"),
        ("CC(C)N=NC", "methyl(propan-2-yl)diazene"),
    ],
)
def test_branched_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(CC1)N=NC")


def test_aromatic_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1N=NC")


def test_chloroethyldiazene_name():
    # PubChem auto-generated name matches exactly.
    assert smiles_to_iupac("ClCCN=N") == "2-chloroethyldiazene"


def test_bis_chloroethyldiazene_name():
    # Two identical compound (halogen-bearing) substituents combine with
    # the compound 'bis' multiplying prefix (P-14.2.2), not the plain
    # 'di' used for simple substituents. PubChem auto-generated name
    # matches exactly.
    assert smiles_to_iupac("ClCCN=NCCCl") == "bis(2-chloroethyl)diazene"


def test_halogen_on_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClN=N")


def test_mixed_compound_and_simple_substituent():
    # Two different substituents where one is a compound (halogen-
    # bearing) name (not independently PubChem-registered for this exact
    # structure, CID 0, but an accepted, reviewed result following the
    # same mechanism confirmed elsewhere). alpha_sort_key strips the
    # leading '2-' locant from '2-chloroethyl', leaving 'chloroethyl' <
    # 'methyl', so the chloroethyl group is cited first -- but per
    # P-16.5.1.3.1 ("the first cited substituent never has enclosing
    # marks unless it is a compound substituent group or includes a
    # locant"), a compound first substituent still needs its own
    # parentheses.
    assert smiles_to_iupac("ClCCN=NC") == "(2-chloroethyl)(methyl)diazene"
