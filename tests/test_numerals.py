import pytest

from smiles_to_iupac._numerals import alkane_name, alkyl_name, multiplying_prefix, numerical_term


@pytest.mark.parametrize(
    "n,expected",
    [
        (1, "mono"),
        (2, "di"),
        (7, "hepta"),
        (10, "deca"),
        (11, "undeca"),
        (12, "dodeca"),
        (14, "tetradeca"),
        (20, "icosa"),
        (21, "henicosa"),
        (22, "docosa"),
        (23, "tricosa"),
        (24, "tetracosa"),
        (41, "hentetraconta"),
        (52, "dopentaconta"),
        (70, "heptaconta"),
        (100, "hecta"),
        (111, "undecahecta"),
        (200, "dicta"),
        (363, "trihexacontatricta"),
        (486, "hexaoctacontatetracta"),
    ],
)
def test_numerical_term(n, expected):
    assert numerical_term(n) == expected


@pytest.mark.parametrize(
    "n,expected",
    [
        (1, "methane"),
        (2, "ethane"),
        (3, "propane"),
        (4, "butane"),
        (5, "pentane"),
        (7, "heptane"),
        (11, "undecane"),
        (20, "icosane"),
        (23, "tricosane"),
        (70, "heptacontane"),
    ],
)
def test_alkane_name(n, expected):
    assert alkane_name(n) == expected


@pytest.mark.parametrize(
    "n,expected",
    [
        (1, "methyl"),
        (2, "ethyl"),
        (3, "propyl"),
        (4, "butyl"),
        (5, "pentyl"),
        (10, "decyl"),
    ],
)
def test_alkyl_name(n, expected):
    assert alkyl_name(n) == expected


@pytest.mark.parametrize(
    "n,compound,expected",
    [
        (2, False, "di"),
        (3, False, "tri"),
        (2, True, "bis"),
        (3, True, "tris"),
        (4, True, "tetrakis"),
        (5, True, "pentakis"),
        (6, True, "hexakis"),
        (10, True, "decakis"),
    ],
)
def test_multiplying_prefix(n, compound, expected):
    assert multiplying_prefix(n, compound=compound) == expected
