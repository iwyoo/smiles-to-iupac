import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_benzo_f_quinoline():
    # PubChem CID 6796, structure-verified.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC=N3") == "benzo[f]quinoline"


def test_benzo_g_quinoline():
    # PubChem CID 520238, structure-verified.
    assert smiles_to_iupac("C1=CC=C2C=C3C(=CC2=C1)C=CC=N3") == "benzo[g]quinoline"


def test_benzo_h_quinoline():
    # PubChem CID 9191, structure-verified.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2N=CC=C3") == "benzo[h]quinoline"


def test_benzo_f_isoquinoline():
    # PubChem CID 123043, structure-verified.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CN=C3") == "benzo[f]isoquinoline"


def test_benzo_g_isoquinoline():
    # PubChem CID 601692, structure-verified -- also the Blue Book's own
    # P-25.3.1.3 worked example.
    assert smiles_to_iupac("C1=CC=C2C=C3C=NC=CC3=CC2=C1") == "benzo[g]isoquinoline"


def test_benzo_h_isoquinoline():
    # PubChem CID 160447, structure-verified.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=NC=C3") == "benzo[h]isoquinoline"


def test_bare_quinoline_and_isoquinoline_still_work():
    assert smiles_to_iupac("c1ccc2ncccc2c1") == "quinoline"
    assert smiles_to_iupac("c1ccc2cnccc2c1") == "isoquinoline"


def test_substituted_variant_names_with_locant():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2N=CC(C)=C3") == "3-methylbenzo[h]quinoline"
