import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Symmetric ethers: P-14.3.4.2(b)'s omitted-locant convention for a
        # homogeneous two-carbon chain with exactly one substituent applies
        # here too, cross-checked against PubChem.
        ("COC", "methoxymethane"),
        ("CCOCC", "ethoxyethane"),
        # Asymmetric ethers: the longer chain is the parent (P-44.3), the
        # shorter side becomes the 'oxy' prefix (P-63.2.1). Cross-checked
        # against PubChem.
        ("COCC", "methoxyethane"),
        ("COCCC", "1-methoxypropane"),
        ("CCOCCC", "1-ethoxypropane"),
        # Longer unbranched alkoxy: no contracted retained name past
        # 'butoxy' (P-63.2.2.1.1), so the plain alkyl name plus 'oxy' is
        # used unchanged.
        ("COCCCCC", "1-methoxypentane"),
        # A branched parent chain is unaffected by this module's R'-only
        # restriction; only the (unbranched) shorter side is the prefix.
        ("CC(C)OCC", "2-ethoxypropane"),
    ],
)
def test_ether(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(C)OC(C)C",  # both sides branched and tied in length
        "CC(C)OC(C)(C)C",  # branched, shorter side is also branched
    ],
)
def test_ether_branched_prefix_out_of_scope(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
