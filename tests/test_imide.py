import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A symmetric imide is "N-acyl amide" substitutive nomenclature
        # (P-66.6.3, see _imide.py docstring): the acyl part reuses
        # _carboxylic_acid.py's already-verified '-oyl' stem and the base
        # part reuses _amide.py's already-verified '-amide' stem, joined by
        # a fixed 'N-' prefix -- no new locant/alphabetization judgment.
        ("O=CNC=O", "N-methanoylmethanamide"),
        ("CC(=O)NC(=O)C", "N-ethanoylethanamide"),
        ("CCC(=O)NC(=O)CC", "N-propanoylpropanamide"),
        ("CCCC(=O)NC(=O)CCC", "N-butanoylbutanamide"),
    ],
)
def test_symmetric_imide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_symmetric_imide():
    assert smiles_to_iupac("CC(C)C(=O)NC(=O)C(C)C") == "N-(2-methylpropanoyl)-2-methylpropanamide"


def test_halogen_substituent():
    assert smiles_to_iupac("ClCC(=O)NC(=O)CCl") == "N-(2-chloroethanoyl)-2-chloroethanamide"


def test_unsymmetric_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC(=O)CC")


def test_n_substituted_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(C)C(=O)C")


def test_unsaturated_imide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(=O)NC(=O)C=C")


def test_cyclic_symmetric_imide_is_named_via_ketone_suffix():
    # An unsubstituted, symmetric cyclic imide (succinimide) doesn't use
    # this module's acyclic "N-acyl amide" construction at all -- P-66.6.3
    # only covers the acyclic case. It's routed to `_ketone.py`'s
    # hetero-ring ketone path instead, which already produces
    # the correct PIN as a plain ring dione. PubChem-verified: CID 11439.
    assert smiles_to_iupac("O=C1CCC(=O)N1") == "pyrrolidine-2,5-dione"


def test_n_substituted_cyclic_imide_named_via_ketone_module():
    # Routed to `_ketone.py`'s hetero-ring ketone path (same as the plain
    # cyclic imide), which supports a single plain alkyl substituent on
    # the ring nitrogen. PubChem-verified: N-methylsuccinimide (CID 11621).
    assert smiles_to_iupac("CN1C(=O)CCC1=O") == "1-methylpyrrolidine-2,5-dione"
