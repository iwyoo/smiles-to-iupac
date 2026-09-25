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


def test_ring_substituted_histidine_still_unsupported():
    # A substituent on the imidazole ring itself (would need the special
    # pi/tau numbering) is out of scope for this step, deferred to a
    # follow-up.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cn1cnc(CC(N)C(=O)O)c1")


def test_phenylalanine_still_resolves():
    # A plain benzene side chain, not imidazole -- unaffected by this
    # module, still routed to `_carboxylic_acid_amine.py`'s existing
    # phenyl-chain path.
    assert smiles_to_iupac("c1ccccc1CC(N)C(=O)O") == "2-amino-3-phenylpropanoic acid"


def test_tryptophan_still_unsupported():
    # A bicyclic (indole) side chain is out of scope here, deferred to a
    # separate later step.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O")


def test_plain_imidazole_still_resolves():
    # An imidazole ring with no amino-acid backbone at all must still
    # route to the generic heteroaromatic-ring naming unchanged.
    assert smiles_to_iupac("c1c[nH]cn1") == "1H-imidazole"
