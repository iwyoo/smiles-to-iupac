import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_terpyrrole_no_indicated_hydrogen():
    # Every ring's own junction sits at its own N -- mirrors
    # `1,1'-bipyrrole`'s own N=2 case (#949): no indicated hydrogen needed
    # anywhere.
    assert (
        smiles_to_iupac("c1ccn(-n2ccc(-n3cccc3)c2)c1")
        == "11,21:23,31-terpyrrole"
    )


def test_terpyrrole_middle_ring_indicated_hydrogen():
    # End rings both attach through their own N (no indicated hydrogen
    # needed); the middle ring's two junctions are both through carbon, so
    # its own N keeps its H and is cited at the front.
    assert (
        smiles_to_iupac("c1ccn(-c2cc(-n3cccc3)c[nH]2)c1")
        == "21H-11,22:24,31-terpyrrole"
    )


def test_quaterpyrrole_two_middle_rings_indicated_hydrogen():
    # Both middle rings attach through carbon only; both end rings attach
    # through their own N.
    assert (
        smiles_to_iupac("c1ccn(-c2cc(-c3cc(-n4cccc4)c[nH]3)c[nH]2)c1")
        == "21H,31H-11,22:24,32:34,41-quaterpyrrole"
    )


def test_two_pyrrole_rings_routes_to_ring_assembly_module():
    # N=2 stays `_ring_assembly.py`'s own job (P-28.2.1's primed-locant
    # scheme, not this module's P-28.3 composite-locant one).
    assert smiles_to_iupac("c1ccn(-n2cccc2)c1") == "1,1'-bipyrrole"


def test_mixed_pyrrole_and_furan_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccn(-c2ccc(-c3ccco3)[nH]2)c1")
