import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure

P = "P(c1ccccc1)(c1ccccc1)c1ccccc1"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[NH3][Pt@SP1](Cl)(Cl)[NH3]", "(SP-4-2)-diamminedichloridoplatinum"),
        ("[NH3][Pt@SP1](Cl)([NH3])Cl", "(SP-4-1)-diamminedichloridoplatinum"),
        ("[NH3][Co@OH1](Cl)([NH3])(Cl)([NH3])Cl", "(OC-6-21)-triamminetrichloridocobalt"),
        ("[NH3][Co@OH4](Cl)([NH3])(Cl)([NH3])Cl", "(OC-6-22)-triamminetrichloridocobalt"),
        (
            f"[O+]#[C][Fe@TB17]([C]#[O+])([C]#[O+])({P}){P}",
            "(TBPY-5-11)-tricarbonylbis(triphenylphosphane)iron",
        ),
        ("[O+]#[C][Fe@TB1]([C]#[O+])(I)(" + P + ")Cl", "(TBPY-5-24-C)-dicarbonylchloridoiodido(triphenylphosphane)iron"),
        ("I[Fe@](Cl)(Br)F", "(T-4-R)-bromidochloridofluoridoiodidoiron"),
        ("I[Fe@@](Cl)(Br)F", "(T-4-S)-bromidochloridofluoridoiodidoiron"),
    ],
)
def test_configuration_index_and_chirality(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_geometry_has_no_descriptor():
    assert smiles_to_iupac("[NH3][Pt]([NH3])(Cl)Cl") == "diamminedichloridoplatinum"


def _co_smiles(cls):
    return f"[Co@OH{cls}]%10%11%12%13%14%15.Br%10.Br%11.[NH2]%12CC[NH2]%13.[NH3]%14.[NH3]%15"


@pytest.mark.parametrize(
    "cls,expected",
    [
        (1, "(OC-6-32-A)-diamminedibromido(ethane-1,2-diamine-κ2N,N')cobalt"),
        (2, "(OC-6-32-C)-diamminedibromido(ethane-1,2-diamine-κ2N,N')cobalt"),
        (12, "(OC-6-22)-diamminedibromido(ethane-1,2-diamine-κ2N,N')cobalt"),
        (25, "(OC-6-13)-diamminedibromido(ethane-1,2-diamine-κ2N,N')cobalt"),
    ],
)
def test_octahedral_chelate_configuration_index(cls, expected):
    assert smiles_to_iupac(_co_smiles(cls)) == expected


def _bis_tridentate(cls):
    a = "[NH2]%10CC[NH]%11CC[NH2]%12"
    b = "[NH2]%13CC[NH]%14CC[NH2]%15"
    return f"[Co@OH{cls}]%10%11%12%13%14%15.{a}.{b}"


@pytest.mark.parametrize("cls,prefix", [(4, "(OC-6-1′2′)-"), (19, "(OC-6-1′2)-"), (1, "(OC-6-2′2-C)-"), (2, "(OC-6-2′2-A)-")])
def test_bis_tridentate_uses_the_priming_convention(cls, prefix):
    assert smiles_to_iupac(_bis_tridentate(cls)).startswith(prefix)
