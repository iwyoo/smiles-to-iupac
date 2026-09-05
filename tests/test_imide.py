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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Both acyl chains end in a plain benzene ring (module docstring's
        # "aromatic ring substituent" extension). PubChem PUG REST confirms
        # the structure/naming pattern (CID 220148, CID 53874044), though
        # its raw string uses retained names ('acetamide') and a different
        # prefix order than this project's own established convention
        # (`_amide.py` already gives 'c1ccccc1CC(=O)N' as
        # '2-phenylethanamide', not '2-phenylacetamide') -- expected values
        # below follow that existing systematic-naming/N-acyl-first
        # convention instead of PubChem's raw string.
        ("c1ccccc1CC(=O)NC(=O)Cc1ccccc1", "N-(2-phenylethanoyl)-2-phenylethanamide"),
        ("c1ccccc1CCC(=O)NC(=O)CCc1ccccc1", "N-(3-phenylpropanoyl)-3-phenylpropanamide"),
    ],
)
def test_symmetric_imide_phenyl_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_imide_phenyl_chain_unsymmetric_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)NC(=O)CCc1ccccc1")


def test_imide_phenyl_ring_directly_on_acyl_carbon_raises():
    # Benzoyl-style ring attachment (chain length 1) uses a separate
    # construction, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)NC(=O)c1ccccc1")


def test_imide_one_sided_phenyl_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)NC(=O)CC")


def test_imide_substituted_benzene_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)NC(=O)Cc1ccccc1C")
