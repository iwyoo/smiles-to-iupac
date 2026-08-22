"""Naming of acyclic hydrocarbons containing exactly one carbon-carbon double or
triple bond (alkenes and alkynes), per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-31.1.1.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): the
  presence of a double or triple bond in an otherwise saturated chain is
  denoted by changing the 'ane' ending of the parent hydride name to 'ene' or
  'yne'. Only the lower locant of the multiple bond is cited, e.g.
  'hex-2-ene (PIN)', 'but-2-yne (PIN) (not dimethylacetylene)'.
- P-31.1.2.1 (Chapter P-3): 'acetylene' is the retained preferred IUPAC name
  for unsubstituted HC#CH (substitution of any kind is disallowed for this
  retained name); once substituted or chain-extended, the systematic '-yne'
  name is used instead (e.g. 'but-2-yne', not 'dimethylacetylene').
- P-14.3.4.2(d) (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf):
  the locant '1' is omitted only for unsubstituted, two- or three-carbon
  alkenes/alkynes ('ethene (PIN)', 'acetylene (PIN)', 'propene (PIN)',
  'propyne (PIN)'). P-14.3.3 (same chapter): as soon as the name carries any
  other locant (i.e. any substituent is present) or the chain is four carbons
  or longer, the multiple bond's own locant becomes essential and must be
  cited (by analogy with '2-chloroethan-1-ol (PIN)', not '2-chloroethanol').
- P-44.3 / P-44.4.1 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf):
  chain length is chosen first (as in `_acyclic.py`); among chains tied for
  greatest length, the principal chain is the one that includes the
  double/triple bond.
- P-14.4(e) (Chapter P-1): the numbering direction is then chosen to give the
  lowest locant to the double/triple bond, ahead of substituent locants.
- P-45.2 (Chapter P-4) / P-14.4(f,g) (Chapter P-1): any remaining tie is
  broken exactly as in `_acyclic.py`: by substituent count, then lowest
  locant set, then lowest locants in citation order.
- P-29.4 / P-46 (Chapter P-2, P-4): substituents on the chain are named the
  same way as for alkanes/cycloalkanes, via `_substituents.py`, so simple
  *and* branched ("compound") substituents (e.g. isopropyl) are supported
  here too.

Two or more multiple bonds, a double bond together with a triple bond, a
multiple bond located in a substituent rather than the principal chain, and
rings are all out of scope for this module and raise `UnsupportedStructure`.
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

_SUFFIX_BY_ORDER = {2.0: "ene", 3.0: "yne"}


def _longest_chains(graph):
    """All maximum-length simple paths in the tree (P-44.3: greater number of
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


def _bond_locant(chain, bond_atoms):
    """1-based locant of the lower-numbered atom of the multiple bond under
    this chain ordering, or None if the bond isn't an edge of this chain."""
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def _name_from_substituents(chain_length, bond_locant, suffix, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-3]
    if chain_length <= 3 and not prefix:
        # P-14.3.4.2(d): the locant is omittable only when unsubstituted and
        # the chain is short enough that no other position is possible.
        if chain_length == 2 and suffix == "yne":
            return "acetylene"
        return stem + suffix
    return prefix + f"{stem}-{bond_locant}-{suffix}"


def _candidate_key(chain_length, bond_locant, suffix, substituents):
    """Sort key implementing P-14.4(e) then P-45.2.1-P-45.2.3, most-preferred
    first."""
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, bond_locant, suffix, grouped)
    # The multiple-bond locant outranks every substituent criterion (P-14.4e);
    # substituent count and lower locants are preferred, so negate the count
    # to sort every field in ascending "most preferred first" order.
    return (bond_locant, -total_count, locant_set, citation_locants, name), name


def name_acyclic_unsaturated(mol) -> str:
    validate_atoms_and_bonds(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported by this module (see smiles_to_iupac._cyclic)"
        )

    bonds = non_single_bonds(mol)
    if len(bonds) != 1 or bonds[0][2] not in _SUFFIX_BY_ORDER:
        raise UnsupportedStructure(
            "this module handles exactly one carbon-carbon double or triple "
            "bond; more than one multiple bond, or a double and a triple bond "
            "together, is not supported yet (see P-31.1.1.1)"
        )
    bond_atoms = bonds[0][0], bonds[0][1]
    suffix = _SUFFIX_BY_ORDER[bonds[0][2]]

    graph = adjacency(mol)
    chains = _longest_chains(graph)
    chain_length = len(chains[0])

    # P-44.3 / P-44.4.1: among the longest chains, only those containing the
    # multiple bond can be the principal chain.
    chains_with_bond = [c for c in chains if _bond_locant(c, bond_atoms) is not None]
    if not chains_with_bond:
        raise UnsupportedStructure(
            "the multiple bond does not lie on a longest chain; expressing it "
            "in a substituent (an alkenyl/alkynyl prefix) is not supported "
            "yet (see P-29.2, P-32.1)"
        )

    best_key = None
    best_name = None
    for chain in chains_with_bond:
        for candidate in (chain, list(reversed(chain))):
            bond_locant = _bond_locant(candidate, bond_atoms)
            substituents = _substituents_for_chain(graph, candidate)
            key, name = _candidate_key(chain_length, bond_locant, suffix, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    return best_name
