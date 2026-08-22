"""Naming of acyclic saturated hydrocarbons (alkanes) whose branches are themselves
unbranched, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-44.3 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the principal
  chain is the one with the greater number of skeletal atoms.
- P-45.2 (same chapter): remaining ties are broken, in order, by (1) the maximum
  number of substituent prefixes, (2) the lowest locant set for those prefixes, and
  (3) the lowest locants in the prefixes' order of citation.
- P-14.3.5 / P-14.4 / P-14.5 (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf):
  lowest-locant-set comparison, numbering, and alphanumerical order of prefixes.
- P-29.3.2.1 (Chapter P-2): unbranched substituent groups (methyl, ethyl, propyl, ...).

Rings, unsaturation, heteroatoms, and branched ("compound") substituents are out of
scope for this module and raise NotImplementedError.
"""

from rdkit import Chem

from ._numerals import alkane_name, alkyl_name, numerical_term


class UnsupportedStructure(NotImplementedError):
    pass


def _validate_carbon_skeleton_tree(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms are not supported yet (see P-21.2.3, skeletal "
                "replacement nomenclature)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure(
                "unsaturation is not supported yet (see P-31.1, alkenes and alkynes)"
            )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported yet (see P-22/P-23, cyclic and polyalicyclic parent hydrides)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )


def _adjacency(mol):
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[a].append(b)
        graph[b].append(a)
    return graph


def _bfs(graph, start):
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


def _path_between(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return list(reversed(path))


def _longest_chains(graph):
    """All maximum-length simple paths in the tree (P-44.3.2: greater number of
    skeletal atoms)."""
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = _bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(_path_between(parents[u], u, v))
    return chains


def _linear_branch(graph, root, coming_from):
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


def _substituents_for_chain(graph, chain):
    """Return {position (1-based) -> [substituent lengths]} for a candidate chain,
    or None if a compound substituent makes this chain unusable."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        lengths = []
        for root in branch_roots:
            length = _linear_branch(graph, root, atom)
            if length is None:
                return None
            lengths.append(length)
        substituents[position] = lengths
    return substituents


def _lowest_locant_set(locants):
    return tuple(sorted(locants))


def _name_from_substituents(chain_length, substituents_by_name):
    """substituents_by_name: {alkyl_name: [locants]} -> full alkane name."""
    prefixes = []
    for name in sorted(substituents_by_name):
        locants = sorted(substituents_by_name[name])
        multiplier = numerical_term(len(locants)) if len(locants) > 1 else ""
        prefixes.append(f"{','.join(str(loc) for loc in locants)}-{multiplier}{name}")
    return "-".join(prefixes) + alkane_name(chain_length) if prefixes else alkane_name(chain_length)


def _candidate_key(chain_length, substituents):
    """Sort key implementing P-45.2.1-P-45.2.3, most-preferred first."""
    substituents_by_name = {}
    for position, lengths in substituents.items():
        for length in lengths:
            substituents_by_name.setdefault(alkyl_name(length), []).append(position)

    total_count = sum(len(v) for v in substituents_by_name.values())
    locant_set = _lowest_locant_set(
        loc for locants in substituents_by_name.values() for loc in locants
    )
    citation_locants = tuple(
        loc
        for name in sorted(substituents_by_name)
        for loc in sorted(substituents_by_name[name])
    )
    name = _name_from_substituents(chain_length, substituents_by_name)
    # Higher substituent count and lower locants are preferred, so negate the count
    # to sort every field in ascending "most preferred first" order.
    return (-total_count, locant_set, citation_locants, name), name


def name_acyclic_alkane(mol) -> str:
    _validate_carbon_skeleton_tree(mol)

    if mol.GetNumAtoms() == 1:
        return alkane_name(1)

    graph = _adjacency(mol)
    chains = _longest_chains(graph)
    chain_length = len(chains[0])

    best_key = None
    best_name = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            substituents = _substituents_for_chain(graph, candidate)
            if substituents is None:
                continue
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    if best_name is None:
        raise UnsupportedStructure(
            "branched (compound) substituent groups are not supported yet (see P-29.4)"
        )
    return best_name
