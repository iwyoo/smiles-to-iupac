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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem PUG REST confirmed: `c1ccccc1N=N` -> "phenyldiazene"
        # (CID 141902), `c1ccccc1N=Nc1ccccc1` -> "diphenyldiazene" (CID
        # 2272), `c1ccccc1N=NC` -> "methyl(phenyl)diazene" (CID 6451992).
        ("c1ccccc1N=N", "phenyldiazene"),
        ("c1ccccc1N=Nc1ccccc1", "diphenyldiazene"),
        ("c1ccccc1N=NC", "methyl(phenyl)diazene"),
        # A plain saturated ring substituent (e.g. cyclohexyl) comes for
        # free from the same `name_branch` fix, since that helper already
        # recognizes both ring shapes -- PubChem PUG REST confirmed:
        # `C1CCCCC1N=N` -> "cyclohexyldiazene" (CID 19772949),
        # `C1CCC(CC1)N=NC` -> "cyclohexyl(methyl)diazene" (CID 21031023).
        ("C1CCCCC1N=N", "cyclohexyldiazene"),
        ("C1CCC(CC1)N=NC", "cyclohexyl(methyl)diazene"),
    ],
)
def test_ring_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_phenyl():
    assert smiles_to_iupac("Cc1ccccc1N=N") == "(2-methylphenyl)diazene"


def test_chloroethyldiazene_name():
    # PubChem auto-generated name matches exactly.
    assert smiles_to_iupac("ClCCN=N") == "(2-chloroethyl)diazene"


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
