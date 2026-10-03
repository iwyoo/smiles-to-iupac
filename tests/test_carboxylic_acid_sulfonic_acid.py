import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the 'sulfo'
        # prefix). PubChem CID 18466245 gives the retained name
        # "sulfoformic acid"; systematic-stem substitution as below.
        ("OC(=O)S(=O)(=O)O", "sulfomethanoic acid"),
        # PubChem CID 31257 gives the retained-name "2-sulfoacetic acid";
        # this codebase always uses the systematic '...oic acid' stem
        # instead of a retained name once any substituent is present
        # (see `_carboxylic_acid.py`'s own established convention, e.g.
        # 'ClCC(=O)O' -> '2-chloroethanoic acid' not '2-chloroacetic
        # acid'), so 'ethanoic' replaces 'acetic' here too.
        ("OC(=O)CS(=O)(=O)O", "2-sulfoethanoic acid"),
        # PubChem CID 409694.
        ("OC(=O)CCS(=O)(=O)O", "3-sulfopropanoic acid"),
        # PubChem CID 21206169.
        ("OC(=O)CCCS(=O)(=O)O", "4-sulfobutanoic acid"),
        # Halogen coexisting alongside the pair. PubChem CID 4308193 gives
        # "2-chloro-2-sulfoacetic acid"; same systematic-stem substitution
        # as above.
        ("OC(=O)C(Cl)S(=O)(=O)O", "2-chloro-2-sulfoethanoic acid"),
        # Multiple sulfonic acids (multiplying prefix). PubChem CID
        # 89390889.
        ("OC(=O)CC(S(=O)(=O)O)CS(=O)(=O)O", "3,4-disulfobutanoic acid"),
    ],
)
def test_carboxylic_acid_sulfonic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_carboxylic_acid_still_routes_normally():
    # No sulfonic acid present -- must still reach `_carboxylic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_plain_sulfonic_acid_still_routes_normally():
    # No carboxylic acid present -- must still reach `_sulfonic_acid.py`'s
    # own path.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_multiple_carboxylic_acids():
    assert smiles_to_iupac("OC(=O)CC(S(=O)(=O)O)C(=O)O") == "2-sulfobutanedioic acid"


def test_other_heteroatom():
    assert smiles_to_iupac("OC(=O)C(O)CS(=O)(=O)O") == "2-hydroxy-3-sulfopropanoic acid"


def test_unsaturated_chain():
    assert smiles_to_iupac("OC(=O)C=CCS(=O)(=O)O") == "4-sulfobut-2-enoic acid"


def test_ring():
    assert smiles_to_iupac("OC(=O)C1CCC(S(=O)(=O)O)CC1") == "4-sulfocyclohexane-1-carboxylic acid"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the carboxylic acid/
        # sulfonic acid chain, mirroring `_carboxylic_acid.py`'s/
        # `_sulfonic_acid_thiol.py`'s own benzene-ring-substituent path,
        # cross-checked against PubChem PUG REST.
        ("c1ccccc1C(S(=O)(=O)O)CC(=O)O", "3-phenyl-3-sulfopropanoic acid"),  # CID 20265801
        ("c1ccccc1CC(S(=O)(=O)O)C(=O)O", "3-phenyl-2-sulfopropanoic acid"),  # CID 129643411
        ("c1ccccc1CCC(S(=O)(=O)O)C(=O)O", "4-phenyl-2-sulfobutanoic acid"),  # CID 159879802
    ],
)
def test_phenyl_chain_carboxylic_acid_sulfonic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_ring_with_second_substituent():
    assert smiles_to_iupac("OC(=O)c1ccccc1S(=O)(=O)O") == "2-sulfobenzoic acid"


def test_phenyl_chain_carboxylic_acid_sulfonic_acid_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CC(S(=O)(=O)O)C(=O)O") == "3-(2-ethenylphenyl)-2-sulfopropanoic acid"
