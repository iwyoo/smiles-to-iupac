import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._seniority import SUFFIX_CLASS_RANK, senior_class


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Mononuclear parent (P-14.3.4.2(a), no locant on the -SO3H).
        ("SCS(=O)(=O)O", "sulfanylmethanesulfonic acid"),
        # Two-carbon chain: not the P-14.3.4.2(b) no-locant case, since a
        # thiol substituent is present in addition to the -SO3H suffix.
        ("SCCS(=O)(=O)O", "2-sulfanylethane-1-sulfonic acid"),
        ("SCCCS(=O)(=O)O", "3-sulfanylpropane-1-sulfonic acid"),
        # Halogen coexisting alongside the sulfonic acid/thiol pair.
        ("ClC(S)CS(=O)(=O)O", "2-chloro-2-sulfanylethane-1-sulfonic acid"),
        # Multiple thiols (multiplying prefix).
        ("SCC(S)CS(=O)(=O)O", "2,3-disulfanylpropane-1-sulfonic acid"),
    ],
)
def test_sulfonic_acid_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_sulfonic_acid_still_routes_normally():
    # No thiol present -- must still reach `_sulfonic_acid.py`'s own path.
    assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"


def test_multiple_sulfonic_acids_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OS(=O)(=O)CCS(=O)(=O)O")


def test_other_heteroatom_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCCS(=O)(=O)O")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC(S)S(=O)(=O)O")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(S)(CC1)S(=O)(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the sulfonic acid/thiol
        # chain, mirroring `_sulfonic_acid.py`'s/`_thiol.py`'s own
        # benzene-ring-substituent path, cross-checked against PubChem
        # PUG REST.
        ("c1ccccc1C(S)CCS(=O)(=O)O", "3-phenyl-3-sulfanylpropane-1-sulfonic acid"),  # CID 57312026
        ("c1ccccc1CC(S)CS(=O)(=O)O", "3-phenyl-2-sulfanylpropane-1-sulfonic acid"),  # CID 57473419
    ],
)
def test_phenyl_chain_sulfonic_acid_thiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_ring_with_second_substituent_raises():
    # A benzene ring with two exocyclic attachments (a thiol directly on
    # the ring plus a separate sulfonic-acid-bearing chain) is out of
    # scope for this single-chain-substituent module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Sc1ccccc1CS(=O)(=O)O")


def test_phenyl_chain_sulfonic_acid_thiol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(S)S(=O)(=O)O")


def test_senior_class_ranks_sulfonic_acid_over_alcohol():
    assert senior_class("sulfonic_acid", "alcohol") == "sulfonic_acid"
    assert senior_class("alcohol", "sulfonic_acid") == "sulfonic_acid"
    assert senior_class("carboxylic_acid", "sulfonic_acid") == "carboxylic_acid"


def test_senior_class_rejects_unranked_names():
    with pytest.raises(KeyError):
        senior_class("sulfonic_acid", "not_a_real_class")


def test_suffix_class_rank_has_no_duplicate_ranks_within_distinct_classes():
    # Chalcogen analogues intentionally share a rank (see module
    # docstring); this only guards against an accidental typo assigning
    # two genuinely different classes the same number.
    assert len(SUFFIX_CLASS_RANK) == len(set(SUFFIX_CLASS_RANK.values()))
