import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). The amide carbon is
        # always C1, and its own locant is never cited (P-14.3.3).
        ("C(N)=O", "methanamide"),
        ("CC(N)=O", "ethanamide"),
        ("CCC(N)=O", "propanamide"),
        # A real positional choice for a substituent: the amide carbon fixes
        # C1 regardless.
        ("CC(C)CC(N)=O", "3-methylbutanamide"),
        # -amide combined with existing unsaturation support, cross-checked
        # against PubChem (crotonamide's PIN).
        ("CC=CC(N)=O", "but-2-enamide"),
        # -amide + halogen substituent prefix, cross-checked against
        # PubChem.
        ("NC(=O)CCCl", "3-chloropropanamide"),
    ],
)
def test_amide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_amide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(C)=O")


def test_lactam_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCN1")


def test_diamide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)CC(N)=O")


def test_ester_not_misnamed_as_amide():
    # -COO- (ester, `_ester.py`) is not amide-shaped: its carbonyl carbon's
    # other oxygen neighbor is carbon-bonded, not a nitrogen, so this must
    # not be routed here and misnamed.
    assert smiles_to_iupac("CC(=O)OC") == "methyl ethanoate"


def test_carboxylic_acid_not_misnamed_as_amide():
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_aryl_amide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)c1ccccc1")


def test_alcohol_mix_names_hydroxy_prefix():
    # 'amide' outranks 'ol' in Table 3.3, so a coexisting standalone -OH is
    # cited as the 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("NC(=O)CCO") == "3-hydroxypropanamide"


def test_amide_enol_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(N)=O")
