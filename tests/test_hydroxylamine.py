import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 787: unsubstituted hydroxylamine (P-68.3.1.1.1).
        ("NO", "hydroxylamine"),
        # PubChem CID 4113: O-substitution only (P-68.3.1.1.1.2's own
        # worked example, "O-methylhydroxylamine (PIN)").
        ("CON", "O-methylhydroxylamine"),
        # PubChem CID 69357.
        ("CCON", "O-ethylhydroxylamine"),
        # PubChem CID 33240: a longer unbranched O-substituent chain.
        ("CCCON", "O-propylhydroxylamine"),
    ],
)
def test_hydroxylamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_not_supported():
    # CH3-NH-OH's actual PIN is 'N-hydroxymethanamine' (an amine-parent
    # name, P-68.3.1.1.1.1) -- a different mechanism this module doesn't
    # implement (see the module docstring); it must not be misnamed as
    # 'N-methylhydroxylamine'.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNO")


def test_n_disubstituted_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(C)O")


def test_n_o_disubstituted_not_supported():
    # CH3-NH-O-CH3's actual PIN is 'N-methoxymethanamine' (also an
    # amine-parent name, P-68.3.1.1.1.3), not 'N,O-dimethylhydroxylamine'.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNOC")


def test_branched_o_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)ON")


def test_ring_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ONC1CCCCC1")


def test_unsaturated_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CON")


def test_phenyl_o_substituent():
    # A plain, unsubstituted benzene ring bonded directly to the
    # hydroxylamine oxygen (P-2/P-3 aromatic-ring-substituent extension) --
    # named as a "phenyl" O-substituent, not as a chain-parent case (unlike
    # `_thiol.py` and friends, `_hydroxylamine.py`'s R is a functional-class
    # substituent name, so no chain is involved here). PubChem PUG REST:
    # "O-phenylhydroxylamine".
    assert smiles_to_iupac("NOc1ccccc1") == "O-phenylhydroxylamine"


def test_benzyl_o_substituent_still_not_supported():
    # A chain that merely ends in a ring further out (O-benzyl) is a
    # compound O-substituent -- already out of scope for an unrelated
    # reason (mirrors `_ether.py`'s enclosure-mark limitation), unaffected
    # by the new direct-ring exemption.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NOCc1ccccc1")


def test_substituted_phenyl_o_substituent_not_supported():
    # A substituted benzene ring (o-tolyl) doesn't match
    # `is_plain_benzene_ring`, so it's out of scope, same as every other
    # module's plain-benzene-only support.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NOc1ccccc1C")
