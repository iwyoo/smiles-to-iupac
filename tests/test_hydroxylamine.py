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
