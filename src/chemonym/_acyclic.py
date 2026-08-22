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

Branched ("compound") substituents are out of scope for this module and raise
NotImplementedError.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    linear_branch,
    lowest_locant_set,
    path_between,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, alkyl_name, numerical_term


def _longest_chains(graph):
    """All maximum-length simple paths in the tree (P-44.3.2: greater number of
    skeletal atoms)."""
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
            length = linear_branch(graph, root, atom)
            if length is None:
                return None
            lengths.append(length)
        substituents[position] = lengths
    return substituents


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
    locant_set = lowest_locant_set(
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
    validate_atoms_and_bonds(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported by this module (see chemonym._cyclic)"
        )

    if mol.GetNumAtoms() == 1:
        return alkane_name(1)

    graph = adjacency(mol)
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
