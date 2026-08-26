import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # pentaprismane (flagship case): C10H10, ten branch atoms all of
        # degree 3, two pentagons joined by five bridging bonds -- built
        # from scratch via RDKit's RWMol (as `_polycyclic.py`'s own
        # `find_polycyclic_core`/`iter_polycyclic_candidates` were already
        # generic in `ring_count`; the only actual gap was `core.py`'s
        # dispatch loop stopping at 5). Cross-checked externally against
        # PubChem CID 138295's computed IUPAC name via the PUG REST API.
        ("C12C3C4C1C1C2C2C3C4C12", "hexacyclo[4.4.0.0^2,5.0^3,9.0^4,8.0^7,10]decane"),
        # 1-methylpentaprismane: all ten vertices of unsubstituted
        # pentaprismane are equivalent, so the lowest-locant tie-break
        # (P-14.4/P-45.2) puts the substituted one at position 1, exactly
        # as it does for cubane's equivalent vertices (see
        # tests/test_pentacyclic.py). PubChem-verified (CID 86129046).
        ("CC12C3C4C5C6C4C1C6C2C53", "1-methylhexacyclo[4.4.0.0^2,5.0^3,9.0^4,8.0^7,10]decane"),
        # Same, with a halogen substituent -- proves `halogens` threads
        # through find_polycyclic_core/name_polycycloalkane the same way it
        # does for the tricyclic/tetracyclic/pentacyclic cases this engine
        # also covers. Not independently found on PubChem, but exercises
        # the same halogen-substituent machinery already PubChem-verified
        # for the pentacyclic case above (single varying axis: ring count).
        ("ClC12C3C4C5C6C4C1C6C2C53", "1-chlorohexacyclo[4.4.0.0^2,5.0^3,9.0^4,8.0^7,10]decane"),
    ],
)
def test_smiles_to_iupac_hexacyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_pentacyclic_is_not_hexacyclic():
    # cubane: genuinely pentacyclic (cyclomatic number 5, eight branch
    # atoms), still resolved as such unaffected by hexacyclic support.
    assert smiles_to_iupac("C12C3C4C1C1C2C3C41") == "pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]octane"


def test_hexacyclic_propellane_like_degree_four_branch_atom_raises():
    # Cubane with one extra direct bond added between two face-diagonal
    # branch atoms: still eight skeletal atoms, but now cyclomatic number 6
    # (hexacyclic) with two branch atoms pushed to degree 4 by the added
    # bond -- a propellane-like topology, out of scope regardless of ring
    # count (mirrors _tricyclic.py's own degree-4 exclusion for the
    # non-propellane-shaped case; this exact SMILES already existed as a
    # regression test in tests/test_pentacyclic.py before hexacyclic
    # support was added -- moved here since it's specifically a hexacyclic
    # (cyclomatic number 6) case).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C15C2C31C45")


def test_heptacyclic_raises():
    # A genuinely heptacyclic system (cyclomatic number 7) is still out of
    # this task's scope (`core.py`'s plain-hydrocarbon dispatch loop stops
    # at ring_count=6) -- built from scratch via RDKit's RWMol: pentaprismane
    # (ten branch atoms, cyclomatic 6) plus one extra bond between two
    # non-adjacent branch atoms, verified to give cyclomatic number 7 with
    # two branch atoms at degree 4 (so it's also propellane-like, out of
    # scope on that basis too, but included here to document that
    # ring_count=7 itself is simply never attempted).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C15C4C14C3C2C54")
