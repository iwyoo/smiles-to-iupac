"""Naming of linear (unbranched) polyspiro saturated hydrocarbon ring systems
(three or more saturated monocyclic carbocyclic rings connected in a chain by
spiro atoms, each internal ring sharing exactly one spiro atom with each of
its two chain neighbors, each terminal ring exactly one spiro atom), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-24.2.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  "Polyspiro parent hydrides consisting of unbranched assemblies of three or
  more saturated cycloalkane rings are named using the nondetachable
  prefixes 'dispiro', 'trispiro', etc., according to the number of spiro
  atoms present, cited in front of the name of the acyclic hydrocarbon that
  has the same number of skeletal atoms. The von Baeyer spiro descriptor
  indicates the number of carbon atoms linking the spiro atoms by arabic
  numbers that are cited in order, starting at the smaller terminal ring if
  one is smaller, and proceeding consecutively, always by the shorter path,
  to the other terminal ring through each spiro atom and then back to the
  first spiro atom; the numbers are separated by full stops and enclosed in
  square brackets. The compound is numbered in the order in which the
  numbers of the von Baeyer spiro descriptor are cited, including spiro
  atoms when encountered for the first time. Each time a spiro atom is
  reached for the second time its locant, which has already been assigned,
  is cited as a superscript number to the number of the preceding linking
  atoms." (e.g. dispiro[3.2.3^7.2^4]dodecane, trispiro[2.2.2.2^9.2^6.2^3]
  pentadecane).
- P-24.2.2.1: if there is a choice of numbers for the spiro descriptor, the
  smaller numbers are selected -- equivalently, the lowest possible locant
  set is allocated to the spiro atoms.
- P-24.2.2.2: if there is still a choice, the descriptor's own numbers
  (excluding superscripts) are compared in citation order, lowest at the
  first point of difference.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): only once the descriptor itself is
  uniquely fixed does the ordinary lowest-locants-to-substituents tie-break
  apply, exactly as `_spiro.py`/`_cyclic.py` use it.
- P-24.2.0 / P-31.1.5: this module covers only saturated systems; unsaturated
  alicyclic spiro ring systems are out of scope (reusing `_spiro.py`'s
  saturation check/message).
- P-24.2.3: branched polyspiro systems (a spiro atom shared by three or more
  rings, or a ring-adjacency graph that isn't a simple chain) are a
  structurally distinct case, out of scope here: `find_linear_polyspiro_chain`
  below simply returns None for that shape, and `core.py` falls through to
  its generic "not supported" error.
- P-24.2.4: heterocyclic spiro ring systems (skeletal replacement) are out of
  scope; `validate_atoms_and_bonds` (imported from `_common.py`) rejects any
  ring heteroatom other than a halogen substituent.
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the same
  way any other substituent does, exactly as in `_cyclic.py`/`_spiro.py`.

Fused, bridged, branched-polyspiro, and heterocyclic polyspiro ring systems
are out of scope for this module.
"""

from itertools import product

from ._cyclic import _group, _substituents_for_ring
from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, numerical_term
from ._spiro import _walk_ring_from_spiro
from ._substituents import alpha_sort_key, format_substituent_prefixes


def find_linear_polyspiro_chain(mol):
    """Return (ring_order, spiro_atoms) if `mol` is an unbranched chain of
    three or more saturated monocyclic rings connected by spiro atoms (each
    internal ring sharing exactly one spiro atom with each chain neighbor,
    each terminal ring exactly one spiro atom, no atom shared by three or
    more rings, no two rings sharing a bond), else None.

    `ring_order`: tuple of ring indices (into `mol.GetRingInfo().AtomRings()`)
    in chain order. `spiro_atoms`: tuple of atom indices, one shorter than
    `ring_order`; `spiro_atoms[i]` is the atom shared by `ring_order[i]` and
    `ring_order[i + 1]`.
    """
    ring_info = mol.GetRingInfo()
    n = ring_info.NumRings()
    if n < 3:
        return None
    atom_rings = [set(r) for r in ring_info.AtomRings()]
    bond_rings = [set(r) for r in ring_info.BondRings()]
    for i in range(n):
        for j in range(i + 1, n):
            if bond_rings[i] & bond_rings[j]:
                return None

    shared = {}
    atom_ring_count = {}
    for i in range(n):
        for j in range(i + 1, n):
            common = atom_rings[i] & atom_rings[j]
            if len(common) > 1:
                return None
            if len(common) == 1:
                atom = next(iter(common))
                shared[(i, j)] = atom
                atom_ring_count[atom] = atom_ring_count.get(atom, 0) + 1
    if any(count > 2 for count in atom_ring_count.values()):
        return None

    ring_graph = {i: [] for i in range(n)}
    for i, j in shared:
        ring_graph[i].append(j)
        ring_graph[j].append(i)
    degrees = [len(ring_graph[i]) for i in range(n)]
    if any(d > 2 for d in degrees) or degrees.count(1) != 2 or degrees.count(2) != n - 2:
        return None

    start_ring = degrees.index(1)
    order = [start_ring]
    previous = None
    current = start_ring
    while len(order) < n:
        candidates = [x for x in ring_graph[current] if x != previous]
        if len(candidates) != 1 or candidates[0] in order:
            return None
        previous, current = current, candidates[0]
        order.append(current)

    spiro_atoms = []
    for k in range(n - 1):
        i, j = order[k], order[k + 1]
        key = (i, j) if (i, j) in shared else (j, i)
        spiro_atoms.append(shared[key])
    return tuple(order), tuple(spiro_atoms)


def _ring_cyclic_order(graph, ring_atoms, start):
    ring_set = set(ring_atoms)
    order = [start]
    previous = None
    current = start
    while len(order) < len(ring_set):
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous, current = current, next_atom
    return order


def _internal_ring_arcs(graph, ring_atoms, s_a, s_b):
    """Both arcs between an internal ring's two spiro atoms: `arc_ab`
    (interior atoms ordered from next-to-`s_a` to next-to-`s_b`) and
    `arc_ba` (the other arc, ordered from next-to-`s_b` to next-to-`s_a`)."""
    cyc = _ring_cyclic_order(graph, ring_atoms, s_a)
    idx_b = cyc.index(s_b)
    return cyc[1:idx_b], cyc[idx_b + 1:]


def _chain_direction_candidates(ring_order, spiro_atoms, atom_rings):
    """(order, spiro_atoms) pairs for each valid overall chain direction
    (P-24.2.2: numbering starts at the smaller terminal ring if one is
    smaller; both directions are tried only when the terminal rings tie)."""
    reversed_order = tuple(reversed(ring_order))
    reversed_spiro = tuple(reversed(spiro_atoms))
    size_first = len(atom_rings[ring_order[0]])
    size_last = len(atom_rings[ring_order[-1]])
    if size_first < size_last:
        return [(ring_order, spiro_atoms)]
    if size_first > size_last:
        return [(reversed_order, reversed_spiro)]
    return [(ring_order, spiro_atoms), (reversed_order, reversed_spiro)]


def _arc_choice_options(graph, atom_rings, ring_order, spiro_atoms):
    """Per internal ring (chain positions 1..n-2), the forced or free set of
    {0, 1} choices for which physical arc is the outbound ('forward') one:
    0 keeps `arc_ab` outbound, 1 keeps (the reverse of) `arc_ba` outbound
    (P-24.2.2: 'always by the shorter path'; only a genuine tie in arc length
    leaves a choice, resolved later by P-24.2.2.1/P-24.2.2.2)."""
    options = []
    for k in range(1, len(ring_order) - 1):
        ring_atoms = atom_rings[ring_order[k]]
        s_a, s_b = spiro_atoms[k - 1], spiro_atoms[k]
        arc_ab, arc_ba = _internal_ring_arcs(graph, ring_atoms, s_a, s_b)
        if len(arc_ab) < len(arc_ba):
            options.append((0,))
        elif len(arc_ab) > len(arc_ba):
            options.append((1,))
        else:
            options.append((0, 1))
    return options


def _build_sequence(graph, atom_rings, ring_order, spiro_atoms, start_dir, end_dir, arc_choices):
    n = len(ring_order)
    s_first = spiro_atoms[0]
    first_ring_non_spiro = set(atom_rings[ring_order[0]]) - {s_first}
    order0 = _walk_ring_from_spiro(graph, first_ring_non_spiro, s_first, start_dir)
    seq = list(order0) + [s_first]
    descriptor = [len(order0)]
    superscripts = [None]

    backward_arcs = []
    for k in range(1, n - 1):
        ring_atoms = atom_rings[ring_order[k]]
        s_a, s_b = spiro_atoms[k - 1], spiro_atoms[k]
        arc_ab, arc_ba = _internal_ring_arcs(graph, ring_atoms, s_a, s_b)
        if arc_choices[k - 1] == 0:
            forward, backward = arc_ab, arc_ba
        else:
            forward, backward = list(reversed(arc_ba)), list(reversed(arc_ab))
        seq.extend(forward)
        seq.append(s_b)
        descriptor.append(len(forward))
        superscripts.append(None)
        backward_arcs.append(backward)

    s_last = spiro_atoms[-1]
    last_ring_non_spiro = set(atom_rings[ring_order[-1]]) - {s_last}
    order_last = _walk_ring_from_spiro(graph, last_ring_non_spiro, s_last, end_dir)
    seq.extend(order_last)
    descriptor.append(len(order_last))
    superscripts.append(s_last)

    for k in range(n - 2, 0, -1):
        seq.extend(backward_arcs[k - 1])
        descriptor.append(len(backward_arcs[k - 1]))
        superscripts.append(spiro_atoms[k - 1])

    return seq, descriptor, superscripts


def _candidate_key(parent, spiro_locants, descriptor, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = parent if not grouped else format_substituent_prefixes(grouped) + parent
    return (spiro_locants, descriptor, locant_set, citation_locants, name)


def name_linear_polyspiro(mol, chain) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated spiro ring systems are not supported yet (see "
            "P-31.1.5, unsaturated alicyclic spiro ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    atom_rings = mol.GetRingInfo().AtomRings()
    ring_order, spiro_atoms = chain
    spiro_count = len(spiro_atoms)
    prefix = numerical_term(spiro_count) + "spiro"

    best_key = None
    best_name = None
    for order, spiros in _chain_direction_candidates(ring_order, spiro_atoms, atom_rings):
        arc_options = _arc_choice_options(graph, atom_rings, order, spiros)
        first_ring_non_spiro = set(atom_rings[order[0]]) - {spiros[0]}
        last_ring_non_spiro = set(atom_rings[order[-1]]) - {spiros[-1]}
        start_dirs = [a for a in graph[spiros[0]] if a in first_ring_non_spiro]
        end_dirs = [a for a in graph[spiros[-1]] if a in last_ring_non_spiro]
        for start_dir in start_dirs:
            for end_dir in end_dirs:
                for arc_choices in product(*arc_options):
                    seq, descriptor, superscripts = _build_sequence(
                        graph, atom_rings, order, spiros, start_dir, end_dir, arc_choices
                    )
                    locants = {atom: pos for pos, atom in enumerate(seq, start=1)}
                    descriptor_str = ".".join(
                        str(num) if sup is None else f"{num}^{locants[sup]}"
                        for num, sup in zip(descriptor, superscripts)
                    )
                    parent = f"{prefix}[{descriptor_str}]{alkane_name(len(seq))}"
                    spiro_locants = tuple(sorted(locants[s] for s in spiros))
                    substituents = _substituents_for_ring(graph, seq, halogens)
                    key = _candidate_key(parent, spiro_locants, tuple(descriptor), substituents)
                    if best_key is None or key < best_key:
                        best_key, best_name = key, key[-1]

    return best_name
