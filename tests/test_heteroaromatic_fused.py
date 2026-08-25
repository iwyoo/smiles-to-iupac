import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_quinoline():
    # C9H7N, cross-checked against PubChem CID 7047's canonical SMILES.
    assert smiles_to_iupac("c1ccc2ncccc2c1") == "quinoline"


def test_indole():
    # C8H7N, cross-checked against PubChem CID 798's canonical SMILES.
    # The PIN cites the indicated hydrogen explicitly (1H-indole), since
    # the non-aromatic 3H-indole tautomer is a distinct real compound.
    assert smiles_to_iupac("c1ccc2[nH]ccc2c1") == "1H-indole"


def test_benzofuran():
    # 1-benzofuran, C8H6O: cross-checked against PubChem CID 9223's
    # canonical SMILES. The PIN uses a leading numeral (not the bracket
    # letter "benzo[b]furan") to disambiguate from the isomeric
    # isobenzofuran -- see module docstring.
    assert smiles_to_iupac("c1ccc2occc2c1") == "1-benzofuran"


def test_isobenzofuran():
    # 2-benzofuran (isobenzofuran), C8H6O: the oxygen sits at the "meso"
    # position between the two fusion carbons rather than adjacent to one,
    # a distinct real compound from 1-benzofuran.
    assert smiles_to_iupac("c1ccc2cocc2c1") == "2-benzofuran"


def test_benzothiophene():
    # 1-benzothiophene, C8H6S: the sulfur analogue of 1-benzofuran, same
    # numeral-disambiguation convention.
    assert smiles_to_iupac("c1ccc2sccc2c1") == "1-benzothiophene"


def test_isobenzothiophene():
    assert smiles_to_iupac("c1ccc2cscc2c1") == "2-benzothiophene"


def test_substituted_benzofuran_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2occc2c1")


def test_substituted_quinoline_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2ncccc2c1")


def test_substituted_indole_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2[nH]ccc2c1")


def test_3h_indole_tautomer_raises():
    # 3H-indole (indolenine): a distinct, non-aromatic-at-C3 tautomer,
    # not the same molecule as 1H-indole -- must not accidentally match.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1(C=Nc2ccccc12)")


def test_other_heteroaromatic_fused_shape_raises():
    # purine: a different retained-name heteroaromatic fused system, not
    # in this module's scope yet -- must not accidentally match quinoline
    # or indole.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ncc2[nH]cnc2n1")
