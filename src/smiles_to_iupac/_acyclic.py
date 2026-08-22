"""Naming of acyclic saturated hydrocarbons (alkanes), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-44.3 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the
  principal chain is the one with the greater number of skeletal atoms.
- P-45.2 (same chapter): remaining ties are broken, in order, by (1) the
  maximum number of substituent prefixes, (2) the lowest locant set for those
  prefixes, and (3) the lowest locants in the prefixes' order of citation.
- P-14.3.5 / P-14.4 / P-14.5 (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf):
  lowest-locant-set comparison, numbering, and alphanumerical order of prefixes.
- P-29.3.2.1 (Chapter P-2): unbranched substituent groups (methyl, ethyl, propyl, ...).
- P-29.4 / P-46 (Chapter P-2, P-4): branched ("compound") substituent groups,
  e.g. `(1-methylpropyl)` for a sec-butyl-like branch — see `_substituents.py`.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    lowest_locant_set,
    non_single_bonds,
    path_between,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch


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
    """Return {position (1-based) -> [(name, is_compound), ...]} for a
    candidate chain."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        substituents[position] = [name_branch(graph, root, atom) for root in branch_roots]
    return substituents


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, grouped):
    prefix = format_substituent_prefixes(grouped)
    return prefix + alkane_name(chain_length)


def _candidate_key(chain_length, substituents):
    """Sort key implementing P-45.2.1-P-45.2.3, most-preferred first."""
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, grouped)
    # Higher substituent count and lower locants are preferred, so negate the count
    # to sort every field in ascending "most preferred first" order.
    return (-total_count, locant_set, citation_locants, name), name


def name_acyclic_alkane(mol) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see "
            "smiles_to_iupac._unsaturated for alkenes/alkynes)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported by this module (see smiles_to_iupac._cyclic)"
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
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    return best_name
