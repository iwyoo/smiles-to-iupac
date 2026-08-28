import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 9550: mononuclear parent (P-14.3.4.2(a), no locant).
        ("C=[N+]=[N-]", "diazomethane"),
        # PubChem CID 70695: two-carbon chain, locant omitted
        # (P-14.3.4.2(b)).
        ("CC=[N+]=[N-]", "diazoethane"),
        # PubChem CID 12224242: three-carbon chain, terminal position,
        # locant required.
        ("CCC=[N+]=[N-]", "1-diazopropane"),
        # PubChem CID 137678: internal carbon (two skeletal-carbon
        # neighbors plus the diazo group) -- confirms the substituted
        # carbon behaves as an ordinary chain position.
        ("CC(=[N+]=[N-])C", "2-diazopropane"),
        # PubChem CID 53674471: diazo coexisting with a halogen
        # substituent, terminal position.
        ("[N-]=[N+]=CCCl", "1-chloro-2-diazoethane"),
        # PubChem CID 150045656: diazo coexisting with a halogen
        # substituent, internal position.
        ("ClCC(=[N+]=[N-])C", "1-chloro-2-diazopropane"),
    ],
)
def test_diazo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1(=[N+]=[N-])CCCCC1")


def test_unsaturated_chain_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CC=[N+]=[N-]")


def test_two_diazo_groups_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[N-]=[N+]=CC=[N+]=[N-]")
