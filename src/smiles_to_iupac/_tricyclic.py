"""Von Baeyer naming of propellane-type saturated tricyclic hydrocarbons:
two branch atoms, both of degree 4, directly bonded to each other *in
addition to* three bridges between them. The direct bond is itself an
independent secondary bridge (P-23.2.5.1), with length 0 (`0^x,y`, the same
zero-length-bridge notation `_polycyclic.py` uses, e.g.
tricyclo[4.4.0.0^3,8]decane) and its two attachment points are the main
bridgeheads themselves.

The general (non-propellane) tricyclic case -- four branch atoms of degree 3
joined by six bridges, per P-23.2.5.1/.2 -- is handled by
`_polycyclic.find_polycyclic_core`/`name_polycycloalkane` with
`ring_count=3`; see that module's docstring for the shared rules (P-23.2.1,
P-23.2.4, P-23.2.6.2.1/.2.4/.2.5, P-14.4/P-45.2) both this module and that
one rely on.

Flagship validation cases:
- [1.1.1]propellane, tricyclo[1.1.1.0^1,3]pentane (C5H6 -- cross-checked
  against Wikipedia/ChemSpider/ACS "Molecule of the Week").
- [2.2.2]propellane, tricyclo[2.2.2.0^1,4]octane (C8H12 -- cross-checked
  against Wikipedia/Wikidata).

Degree-4-or-higher branch atoms in any topology other than this exact
two-bridgehead propellane shape, and tetracyclic-or-higher propellane-like
systems, remain out of scope and raise UnsupportedStructure -- see
`find_propellane_core`'s degree/branch-atom-count checks.
"""

from itertools import permutations

from ._cyclic import _substituents_for_ring
from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    non_single_bonds,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._polycyclic import _candidate_key, _strip_leaves, _walk_to_branch


def find_propellane_core(mol):
    """Return (bh1, bh2, bridges) if `mol`'s carbon skeleton, after stripping
    acyclic branches, reduces to exactly two branch atoms of degree 4 (and no
    atom of higher degree), directly bonded to each other, joined by exactly
    three further bridges through degree-2 atoms -- else None.

    `bridges` is a list of three `(bh1, bh2, path)` tuples, `path` ordered
    nearest-`bh1`-first. The cyclomatic-number-3 check makes this genuinely
    tricyclic (not tetracyclic-or-higher) for any bridge lengths: with two
    degree-4 branch atoms and three bridges plus the direct bond, edges -
    vertices + 1 is always 3 regardless of bridge length, so no extra guard
    against tetracyclic-or-higher is needed beyond the checks already here."""
    graph = adjacency(mol)
    core = _strip_leaves(graph)
    if not core:
        return None
    vertices = len(core)
    edge_count = sum(len(neighbors) for neighbors in core.values()) // 2
    if edge_count - vertices + 1 != 3:
        return None
    if any(len(neighbors) not in (2, 4) for neighbors in core.values()):
        return None
    branch_atoms = {atom for atom, neighbors in core.items() if len(neighbors) == 4}
    if len(branch_atoms) != 2:
        return None
    bh1, bh2 = branch_atoms
    if bh2 not in core[bh1]:
        return None

    bridges = []
    for first in core[bh1] - {bh2}:
        result = _walk_to_branch(core, branch_atoms, bh1, first)
        if result is None:
            return None
        v, path = result
        if v != bh2:
            return None
        bridges.append((bh1, bh2, path))
    if len(bridges) != 3:
        return None

    return bh1, bh2, bridges


def name_propellane(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated propellane-type ring systems are not supported yet "
            "(see P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    bh1, bh2, bridges = core
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    best_name = None
    for start, other in ((bh1, bh2), (bh2, bh1)):
        oriented_paths = [
            path if start == bh1 else list(reversed(path)) for _, _, path in bridges
        ]
        for perm in permutations(oriented_paths):
            if not (len(perm[0]) >= len(perm[1]) >= len(perm[2])):
                continue
            main_ring_first, main_ring_second, main_bridge = perm
            main_order = (
                [start] + list(main_ring_first) + [other]
                + list(reversed(main_ring_second)) + list(main_bridge)
            )
            position = {atom: idx + 1 for idx, atom in enumerate(main_order)}
            lo, hi = sorted((position[start], position[other]))
            substituents = _substituents_for_ring(graph, main_order, halogens)
            a, b, c = len(main_ring_first), len(main_ring_second), len(main_bridge)
            total_atoms = a + b + c + 2
            # The direct bond between the two bridgeheads is the independent
            # secondary bridge (P-23.2.5.1), length 0, attached at the main
            # bridgeheads themselves -- see module docstring.
            parent = f"tricyclo[{a}.{b}.{c}.0^{lo},{hi}]{alkane_name(total_atoms)}"
            # P-23.2.6.2.4/.2.5 (lowest secondary-bridge locants) then
            # P-14.4/P-45.2 (lowest substituent locants).
            key = ((lo, hi),) + _candidate_key(parent, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    if best_name is None:
        raise UnsupportedStructure(
            "this propellane topology is not supported yet (see "
            "_tricyclic.py's module docstring for the scope this module "
            "covers)"
        )
    return best_name
