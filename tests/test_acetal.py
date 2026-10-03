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


def test_branched_alkoxy_named():
    assert smiles_to_iupac("CCC(OC(C)C)OCC") == "1-ethoxy-1-(propan-2-yloxy)propane"


def test_unsaturated_acetal_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(OCC)OCC")


def test_hemiacetal_named_as_alkoxy_alcohol():
    # RR'C(OH)(O-R'') (P-66.6.5.2) isn't a distinct suffix construction --
    # it's just an alcohol (`_alcohol.py`) with a plain alkoxy ether
    # substituent prefix, same as '2-methoxyethanol'. PubChem CID 93269.
    assert smiles_to_iupac("OC(OCC)CC") == "1-ethoxypropan-1-ol"


def test_ether_not_misnamed_as_acetal():
    # A plain single ether (`_ether.py`) has only one oxygen and must not
    # be routed here.
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_ketone_not_misnamed_as_acetal():
    assert smiles_to_iupac("CC(=O)C") == "propan-2-one"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter on the chain, away
        # from the acetal carbon (P-92) -- CIP computed entirely by
        # RDKit's `rdCIPLabeler`, same pattern as `_carboxylic_acid.py`.
        # PubChem CID 90324169.
        ("CC[C@@H](C)C(OC)OC", "(2R)-1,1-dimethoxy-2-methylbutane"),
        ("CC[C@H](C)C(OC)OC", "(2S)-1,1-dimethoxy-2-methylbutane"),
        # The stereocenter can be the acetal carbon itself, when its two
        # alkoxy groups differ (otherwise it's not a genuine stereocenter
        # at all -- confirmed via RDKit's `Chem.FindPotentialStereo`).
        # PubChem CID 97550416.
        ("CC[C@H](OC)OCC", "(1R)-1-ethoxy-1-methoxypropane"),
    ],
)
def test_acetal_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acetal_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(OC)OC") == "1,1-dimethoxy-2-methylbutane"
