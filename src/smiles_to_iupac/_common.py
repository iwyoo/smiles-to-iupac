"""Shared graph utilities and validation for saturated hydrocarbon parent
hydrides, used by both the acyclic (`_acyclic.py`) and monocyclic (`_cyclic.py`)
naming modules.

- P-35.2.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'fluoro',
  'chloro', 'bromo', and 'iodo' are the preselected substituent prefixes for
  -F, -Cl, -Br, and -I respectively. These halogen atoms are always
  monovalent, non-skeletal substituents (P-44.3's "skeletal atoms" always
  means carbon in this repository) — never part of a parent hydride's
  counted chain/ring — so `carbon_adjacency` below lets chain/ring-skeleton
  search ignore them while substituent detection still finds them.
"""

from rdkit import Chem

from ._numerals import numerical_term


class UnsupportedStructure(NotImplementedError):
    pass


HALOGEN_PREFIXES = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_ALLOWED_ATOMIC_NUMS = {6, *HALOGEN_PREFIXES}

ENE_BOND_ORDER = 2.0
YNE_BOND_ORDER = 3.0


def validate_atoms_and_bonds(mol):
    """Structure-independent checks shared by every parent hydride kind: a
    single-fragment, all-carbon skeleton optionally bearing monovalent
    halogen substituents (P-35.2.1). Bond order (all single, or exactly one
    double/triple bond) and ring shape (none, one simple ring, or more) are
    checked separately by each naming module, since what's allowed there
    differs (see `_acyclic.py`, `_cyclic.py`, `_unsaturated.py`)."""
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than halogen substituents (F, Cl, Br, I; "
                "see P-35.2.1) are not supported yet (see P-21.2.3, skeletal "
                "replacement nomenclature)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
        elif atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1); "
                "polyvalent or bridging halogen structures are not supported"
            )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute (see P-44.3)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )


def non_single_bonds(mol):
    """List of (begin_atom_idx, end_atom_idx, bond_order) for every bond whose
    order isn't 1.0 (single). Used to classify a molecule's degree of
    unsaturation for dispatch (see `core.py`)."""
    return [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0
    ]


def adjacency(mol):
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[a].append(b)
        graph[b].append(a)
    return graph


def carbon_adjacency(mol):
    """Like `adjacency`, but restricted to carbon atoms and the bonds directly
    between them. Chain/ring-skeleton search (P-44.3.2) must use this instead
    of `adjacency` so a terminal halogen substituent is never mistaken for a
    chain-extending skeletal atom; substituent-detection code should keep
    using the full `adjacency(mol)` so it can still find that halogen."""
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def halogen_substituents(mol):
    """{atom_idx -> substituent prefix name} for every halogen atom in `mol`
    (P-35.2.1). Passed down into `name_branch` so it can name a halogen leaf
    and exclude halogens from a compound substituent's own internal chain
    search, the same way `carbon_adjacency` does for a parent hydride."""
    return {
        atom.GetIdx(): HALOGEN_PREFIXES[atom.GetAtomicNum()]
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() in HALOGEN_PREFIXES
    }


def bfs(graph, start):
    dist = {start: 0}
    parent = {start: None}
    queue = [start]
    while queue:
        next_queue = []
        for node in queue:
            for neighbor in graph[node]:
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    next_queue.append(neighbor)
        queue = next_queue
    return dist, parent


def path_between(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return list(reversed(path))


def linear_branch(graph, root, coming_from):
    """Walk a branch outward; return its atom count, or None if it forks
    (a "compound" substituent, P-29.4, not yet supported)."""
    length = 1
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) == 0:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]
        length += 1


def lowest_locant_set(locants):
    return tuple(sorted(locants))


def multiplied_word(count, base):
    """Multiplying-prefix word for `count` occurrences of a suffix like 'ol'/
    'one'/'al' (P-14.2.1): omitted for zero, bare for one, else prefixed with
    the basic numerical term ('di', 'tri', ...)."""
    if count == 0:
        return ""
    if count == 1:
        return base
    return numerical_term(count) + base


def group_substituents(substituents):
    """{position -> [(name, is_compound), ...]} -> {name -> {"locants": [...],
    "compound": bool}}, merging same-named substituents at different
    positions so they can be cited once with a multiplying prefix (P-16)."""
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def longest_chains(graph):
    """Every longest path (by atom count) through an undirected acyclic
    `graph` (P-44.3, the principal chain must be a longest chain candidate)."""
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(path_between(parents[u], u, v))
    return chains


def bond_locant(chain, bond_atoms):
    """1-based position along `chain` of the bond between `bond_atoms` (a
    2-tuple of atom indices), or None if that bond doesn't lie on `chain`."""
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def bond_locants(chain, bonds):
    """(ene_locants, yne_locants) for every (a, b, order) bond in `bonds` that
    lies on `chain`, or None if any bond doesn't (a shorter/wrong chain
    candidate)."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return ene, yne
