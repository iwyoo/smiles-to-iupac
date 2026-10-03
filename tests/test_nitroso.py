import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 70075: mononuclear parent (P-14.3.4.2(a), no locant).
        ("CN=O", "nitrosomethane"),
        # PubChem CID 79124: two-carbon chain, locant omitted (same
        # P-14.3.4.2(b) rule `_nitro.py`'s own docstring confirms for
        # 'nitroethane', without the PubChem-vs-PIN mismatch nitro has).
        ("CCN=O", "nitrosoethane"),
        # PubChem CID 21528903: nitroso coexisting with a halogen
        # substituent.
        ("ClCCN=O", "1-chloro-2-nitrosoethane"),
        # PubChem CID 13116308: two nitroso groups (multiplying prefix).
        ("O=NCCN=O", "1,2-dinitrosoethane"),
    ],
)
def test_nitroso(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_nitro_not_misnamed_as_nitroso():
    assert smiles_to_iupac("C[N+](=O)[O-]") == "nitromethane"


def test_ring():
    assert smiles_to_iupac("C1CCCCC1N=O") == "nitrosocyclohexane"


def test_unsaturated_chain():
    assert smiles_to_iupac("C=CN=O") == "nitrosoethene"


def test_phenyl_nitroso_direct_bond():
    # P-44.1.2.2 rule (1): 'nitroso' has no suffix form, so the ring is
    # always senior to a chain of the same class -- confirmed by PubChem
    # CID 11473.
    assert smiles_to_iupac("c1ccccc1N=O") == "nitrosobenzene"


def test_phenyl_nitroso_chain():
    # PubChem CID 12267972/21470433 give "nitrosomethylbenzene"/
    # "2-nitrosoethylbenzene" (no parentheses); this codebase follows
    # `_azide.py`'s identical, PIN-verified rule instead (its own test
    # confirms "(azidomethyl)benzene" for the same 1-carbon shape) --
    # same PubChem inconsistency already documented for `_nitro.py`
    # (PR #321).
    assert smiles_to_iupac("c1ccccc1CN=O") == "(nitrosomethyl)benzene"
    assert smiles_to_iupac("c1ccccc1CCN=O") == "(2-nitrosoethyl)benzene"


def test_phenyl_nitroso_substituted_ring():
    assert smiles_to_iupac("Cc1ccccc1CN=O") == "1-methyl-2-(nitrosomethyl)benzene"


def test_phenyl_nitroso_unsaturation():
    assert smiles_to_iupac("C=Cc1ccccc1CN=O") == "1-ethenyl-2-(nitrosomethyl)benzene"
