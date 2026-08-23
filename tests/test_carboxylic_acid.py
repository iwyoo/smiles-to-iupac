import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). The -COOH carbon is
        # always C1, and its own locant is never cited (P-14.3.3).
        ("C(=O)O", "methanoic acid"),
        ("CC(=O)O", "ethanoic acid"),
        ("CCC(=O)O", "propanoic acid"),
        # A real positional choice for a substituent, cross-checked against
        # PubChem: the -COOH carbon fixes C1 regardless.
        ("CC(C)CC(=O)O", "3-methylbutanoic acid"),
        # Dicarboxylic acid: both chain ends are -COOH carbons, multiplying
        # prefix + no locants cited, cross-checked against PubChem (adipic
        # acid's PIN).
        ("OC(=O)CCCCC(=O)O", "hexanedioic acid"),
        # -oic acid combined with existing unsaturation support, cross-checked
        # against PubChem (crotonic acid's PIN).
        ("CC=CC(=O)O", "but-2-enoic acid"),
        # Unsaturated dicarboxylic acid, cross-checked against PubChem
        # (fumaric/maleic acid's PIN).
        ("OC(=O)C=CC(=O)O", "but-2-enedioic acid"),
        # -oic acid + halogen substituent prefix, cross-checked against
        # PubChem.
        ("OC(=O)CCCl", "3-chloropropanoic acid"),
        # Diene acid: two multiplied 'ene' locants ahead of the acid suffix,
        # cross-checked against PubChem (sorbic acid's PIN).
        ("CC=CC=CC(=O)O", "hexa-2,4-dienoic acid"),
    ],
)
def test_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ether_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOCC")


def test_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OC")


def test_ring_attached_carboxylic_acid_raises():
    # P-65.1.1.2: a -COOH on a ring uses the separate 'carboxylic acid'
    # suffix construction, out of scope for this acyclic-only module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCCCC1")


def test_aryl_carboxylic_acid_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)c1ccccc1")


def test_amine_coexisting_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(=O)O")


def test_alcohol_mix_names_hydroxy_prefix():
    # '-oic acid' outranks 'ol' in Table 3.3, so a coexisting standalone -OH
    # is cited as the 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("OC(=O)CCO") == "3-hydroxypropanoic acid"


def test_carboxylic_acid_enol_mix_raises():
    # A hydroxyl on a C=C carbon (an enol) is a tautomer of a more senior
    # carbonyl form and out of scope, same as `_alcohol.py`'s own enol check.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(=O)O")
