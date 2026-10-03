import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_terthiophene():
    # Real PubChem structure, CID 65067 (2,2':5',2''-terthiophene).
    assert (
        smiles_to_iupac("c1ccc(-c2ccc(-c3cccs3)s2)s1")
        == "12,22:25,32-terthiophene"
    )


def test_terfuran():
    # Real PubChem structure, CID 10921534 (2,2':5',2''-terfuran).
    assert (
        smiles_to_iupac("C1=COC(=C1)C2=CC=C(O2)C3=CC=CO3")
        == "12,22:25,32-terfuran"
    )


def test_terselenophene():
    assert (
        smiles_to_iupac("c1cc[se]c1-c1ccc([se]1)-c1ccc[se]1")
        == "12,22:25,32-terselenophene"
    )


def test_tertellurophene():
    assert (
        smiles_to_iupac("c1cc[te]c1-c1ccc([te]1)-c1ccc[te]1")
        == "12,22:25,32-tertellurophene"
    )


def test_terpyrimidine():
    assert (
        smiles_to_iupac("c1ccnc(-c2ccnc(-c3ccncn3)n2)n1")
        == "12,24:22,34-terpyrimidine"
    )


def test_terpyridazine():
    assert (
        smiles_to_iupac("c1ccc(-c2ccc(-c3cccnn3)nn2)nn1")
        == "13,23:26,33-terpyridazine"
    )


def test_terpyrazine():
    assert (
        smiles_to_iupac("c1cnc(-c2cnc(-c3cnccn3)cn2)cn1")
        == "12,22:25,32-terpyrazine"
    )


def test_halogen_substituent_on_thiophene_chain():
    assert (
        smiles_to_iupac("c1ccc(-c2ccc(-c3cc(Cl)cs3)s2)s1")
        == "14-chloro-12,22:25,32-terthiophene"
    )


def test_two_thiophene_rings_routes_to_ring_assembly_module():
    # N=2 is `_ring_assembly.py`'s own job (P-28.2.1, not this module's
    # P-28.3 composite-locant scheme) -- confirms the N=3-6 module here no
    # longer claims (and mis-fails on) the N=2 shape.
    assert smiles_to_iupac("c1ccc(-c2cccs2)s1") == "2,2'-bithiophene"


def test_mixed_thiophene_and_furan():
    assert smiles_to_iupac("c1ccc(-c2ccc(-c3ccco3)s2)s1") == "2-([2,2'-bithiophen]-5-yl)furan"
