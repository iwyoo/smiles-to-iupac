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
        # A branched R' (the shorter, prefix side): P-63.2.2.1.1 encloses
        # only R' in parentheses, with 'oxy' outside them. Structure
        # verified against PubChem: CCCCOC(C)C -> CID 137240
        # ("1-propan-2-yloxybutane" -- PubChem's own PIN-style name, this
        # project keeps its usual CAS-style substituent name instead).
        ("CCCCOC(C)C", "1-(1-methylethyl)oxybutane"),
        # Longer R' side also branched: PubChem CID 28488 confirms the
        # structure ("2-methyl-2-propan-2-yloxypropane").
        ("CC(C)OC(C)(C)C", "2-(1-methylethyl)oxy-2-methylpropane"),
    ],
)
def test_ether(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ether_both_sides_branched_and_tied_out_of_scope():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)OC(C)C")
