import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_l_histidine():
    # PubChem CID 6274.
    assert smiles_to_iupac("C1=C(NC=N1)C[C@@H](C(=O)O)N") == "L-histidine"


def test_d_histidine():
    # PubChem CID 71083.
    assert smiles_to_iupac("C1=C(NC=N1)C[C@H](C(=O)O)N") == "D-histidine"


def test_histidine_unspecified_stereocenter_no_ld_prefix():
    assert smiles_to_iupac("C1=C(NC=N1)CC(C(=O)O)N") == "histidine"


def test_ring_substituted_histidine_is_named_as_an_amino_acid_with_a_ring_prefix():
    assert (
        smiles_to_iupac("Cn1cnc(CC(N)C(=O)O)c1")
        == "2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid"
    )


def test_phenylalanine_still_resolves():
    # A plain benzene side chain, not imidazole -- unaffected by this
    # module, resolved by `_amino_acid.py`'s own table now (see
    # test_amino_acid.py) rather than `_carboxylic_acid_amine.py`'s
    # systematic phenyl-chain path.
    assert smiles_to_iupac("c1ccccc1CC(N)C(=O)O") == "phenylalanine"


def test_tryptophan_still_resolves():
    # A bicyclic (indole) side chain -- now resolved by `_amino_acid.py`'s
    # table (see test_amino_acid.py), unaffected by this module.
    assert smiles_to_iupac("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O") == "L-tryptophan"


def test_plain_imidazole_still_resolves():
    # An imidazole ring with no amino-acid backbone at all must still
    # route to the generic heteroaromatic-ring naming unchanged.
    assert smiles_to_iupac("c1c[nH]cn1") == "1H-imidazole"
