"""Naming of ortho-fused mancude (mancude = maximum number of noncumulative
double bonds; here simply "aromatic") six-membered all-carbon ring chains,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-25.1.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): fused
  ring nomenclature for polycyclic mancude ring systems.
- P-25.1.2.1 / P-25.1.2.2: the retained names 'benzene' (one ring, PIN),
  'naphthalene', 'anthracene', 'tetracene', 'pentacene' (the linear
  "polyacene" family, P-25.1.2.1) and 'phenanthrene' (the first
  "polyaphene", n=3, P-25.1.2.2, an *angular* three-ring chain).
- P-25.3.1.3: "ortho-fusion" -- two rings sharing exactly one bond (two
  atoms); a ring atom shared by three or more rings is *peri*-fusion
  (e.g. pyrene, acenaphthylene) and out of scope here.
- P-25.3.2.3.1 / P-25.3.2.3.3 (orientation) and P-25.3.3.1.1 (numbering):
  "The numbering of peripheral atoms in the preferred orientation starts
  from the uppermost ring. If there is more than one uppermost ring, the
  ring furthest to the right is chosen. Numbering starts from the
  nonfused atom most counterclockwise in the ring selected and proceeds
  in a clockwise direction around the system... Each fusion carbon atom
  is given the same number as the immediately preceding nonfusion
  skeletal atom, modified by a Roman letter 'a', 'b', 'c', etc." Because a
  fusion carbon in an all-carbon mancude ring has no free valence, it is
  never a substituent locant, so this module tracks only the plain
  integer locants (1, 2, 3, ...) assigned to nonfusion atoms; the
  lettered fusion locants never need to be materialized as text.
- P-25.3.3: "Anthracene, phenanthrene, acridine, carbazole, xanthene and
  its chalcogen analogues, purine, and cyclopenta[a]phenanthrene are
  exceptions; traditional numberings are retained." Anthracene's
  traditional numbering does NOT follow the general peripheral-walk
  algorithm above -- e.g. positions 9 and 10 (the "meso" carbons) are the
  two ring-2 nonfusion atoms, but they are graph-theoretically nonadjacent
  (confirmed below and against PubChem's own round-tripped IUPAC names for
  1-bromo-, 1,8-dichloro- and 9,10-dibromoanthracene), so a simple
  clockwise walk could never number them consecutively; anthracene's
  numbering instead numbers both terminal rings first (1-4, 5-8), then the
  two meso carbons (9, 10) and remaining fusion carbons (4a, 8a, 9a, 10a).
  Phenanthrene's traditional numbering, by contrast, *does* coincide with
  a plain clockwise peripheral walk (verified the same way, via 1-bromo-,
  9-bromo- and 4,5-dimethylphenanthrene) -- it is an "exception" only in
  the sense that this module does not derive its starting orientation from
  the general P-25.3.2.3.3 criteria (out of scope, see below), not because
  its walk mechanics differ from naphthalene's.
- P-35.2.1 (Chapter P-3): halogen substituents on ring carbons are named
  exactly as for any other parent hydride -- see `_common.py`.

Scope, deliberately narrow (this is a hard, novel domain for this project):
only a simple *chain* of ortho-fused benzenoid rings (each ring fused to at
most two others, no ring-adjacency cycles, no atom shared by three or more
rings) is handled, and only when the whole system carries one of six
retained names:

- 1 ring: benzene.
- 2 rings: naphthalene (the only possible ortho-fused 2-ring shape).
- 3 rings, straight (both fusion bonds of the middle ring are opposite
  hexagon edges): anthracene, traditional numbering.
- 3 rings, angular (bent): phenanthrene, traditional numbering.
- 4 rings, straight: tetracene, general algorithm.
- 5 rings, straight: pentacene, general algorithm.
- Anything else -- 4+ rings that aren't a straight chain (angular
  "polyphenes" beyond phenanthrene, e.g. chrysene/benz[a]anthracene), 6+
  rings, peri-fusion, a branched or cyclic ring-adjacency graph, or any
  heteroatom in the ring system -- raises `UnsupportedStructure`. These
  would require genuine `benzo[x,y-z]fusion[...]` name construction
  (P-25.3.2) or the full P-25.3.2.3.3 orientation search, both out of
  scope for this PR.

Substituents (halogens and simple/compound alkyl groups) reuse the
existing `_substituents.py`/`_cyclic.py` machinery exactly as the other
ring modules do. A fusion carbon has no free valence in a mancude ring, so
substituents only ever attach to a plain-numbered (nonfusion) ring atom;
P-14.3.3's "locant not essential" simplification, however, applies only to
benzene (every ring position is equivalent before substitution) -- for
every other ring size here, ring positions are chemically distinct and a
single substituent's locant is always cited.
"""

from ._cyclic import _group
from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    specified_stereocenters,
    validate_atoms_and_bonds,
)
from ._substituents import alpha_sort_key, branch_atom_locant, format_substituent_prefixes, name_branch


def _ring_cycle(graph, ring_atoms):
    """Order a ring's atoms into a cyclic sequence by walking its bonds
    (works for aromatic rings too -- adjacency doesn't encode bond order)."""
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def _run_between(cycle, a, b):
    """`cycle`: a ring's 6-atom cyclic order, where `a` and `b` are the
    ring's two (adjacent, chord-bonded) fusion atoms. Return the other 4
    (nonfusion) atoms in order, starting with the one adjacent to `a` and
    ending with the one adjacent to `b`."""
    n = len(cycle)
    ia = cycle.index(a)
    step = -1 if cycle[(ia + 1) % n] == b else 1
    return [cycle[(ia + step * k) % n] for k in range(1, 5)]


def find_aromatic_fused_core(mol):
    """Return (atom_rings, ring_atom_sets, fusion_bond_idxs) if every SSSR
    ring in `mol` is a 6-membered all-carbon aromatic ring, else None (not
    an aromatic hydrocarbon ring system at all -- let other dispatch
    branches in core.py handle it, e.g. a heteroaromatic ring falls through
    to name_cycloalkane, which correctly rejects the heteroatom).

    Deliberately permissive about *shape* (peri-fused, branched, too many
    rings, etc. all return a core here) -- `name_aromatic_fused` does that
    validation and raises specific `UnsupportedStructure` messages, mirroring
    every other module's find_*/name_* split in this project. This also
    matters for correctness: an aromatic ring system's carbon skeleton can be
    graph-isomorphic to a *saturated* bicyclic through pentacyclic core
    (e.g. naphthalene's skeleton matches decahydronaphthalene's), so this
    check must run, and must succeed for any valid aromatic ring shape,
    before those saturated modules get a chance to misdetect it and raise a
    confusing "unsaturated ... not supported yet" error instead of this
    module's own clear message -- see core.py's dispatch-order comment.
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return None
    for ring in atom_rings:
        if len(ring) != 6:
            return None
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                return None

    bond_ring_count = {}
    for bonds in ring_info.BondRings():
        for b in bonds:
            bond_ring_count[b] = bond_ring_count.get(b, 0) + 1
    fusion_bond_idxs = {b for b, c in bond_ring_count.items() if c >= 2}
    return atom_rings, [set(r) for r in atom_rings], fusion_bond_idxs


def _validate_aromatic_bonds(mol):
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 1.0 or bond.GetIsAromatic():
            continue
        raise UnsupportedStructure(
            "a non-aromatic multiple bond is not supported within or "
            "attached to an aromatic ring system (see P-25.3.1.3)"
        )


def _ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol):
    """Ring-fusion graph (nodes = ring indices, edge = two rings sharing an
    aromatic bond) plus {frozenset({ring_i, ring_j}) -> (atom_a, atom_b)},
    the shared bond's two atoms, needed later to find each ring's fusion
    atoms without re-deriving them from scratch."""
    n = len(atom_rings)
    adj = {i: set() for i in range(n)}
    fusion_bonds_by_pair = {}
    for bidx in fusion_bond_idxs:
        bond = mol.GetBondWithIdx(bidx)
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        rings_with_bond = [i for i in range(n) if a in ring_atom_sets[i] and b in ring_atom_sets[i]]
        if len(rings_with_bond) != 2:
            raise UnsupportedStructure(
                "an aromatic bond shared by other than exactly two rings is "
                "not supported (see P-25.3.1.3)"
            )
        i, j = rings_with_bond
        adj[i].add(j)
        adj[j].add(i)
        fusion_bonds_by_pair[frozenset((i, j))] = (a, b)
    return adj, fusion_bonds_by_pair


def _atom_ring_membership(atom_rings):
    membership = {}
    for ring_idx, ring in enumerate(atom_rings):
        for atom in ring:
            membership.setdefault(atom, []).append(ring_idx)
    return membership


def _ring_path_order(adj, n):
    """Validate the ring-fusion graph is a simple path (P-25.3.1.3: only a
    chain of ortho-fused rings is in scope -- no branching, no cycle) and
    return the rings in path order."""
    if n == 1:
        return [0]
    if any(len(v) > 2 for v in adj.values()):
        raise UnsupportedStructure(
            "a branched ring-fusion arrangement (a ring ortho-fused to "
            "three or more others, e.g. triphenylene) is not supported "
            "(see P-25.3.1.3, only a simple chain is)"
        )
    endpoints = [i for i in range(n) if len(adj[i]) <= 1]
    total_edges = sum(len(v) for v in adj.values()) // 2
    if len(endpoints) != 2 or total_edges != n - 1:
        raise UnsupportedStructure(
            "a cyclic or disconnected ring-fusion arrangement is not "
            "supported (see P-25.3.1.3, only a simple chain of ortho-fused "
            "rings is)"
        )
    order = [endpoints[0]]
    previous, current = None, endpoints[0]
    while len(order) < n:
        next_ring = next(r for r in adj[current] if r != previous)
        order.append(next_ring)
        previous, current = current, next_ring
    return order


def _edge_index(ring_cycle, a, b):
    n = len(ring_cycle)
    for i in range(n):
        if {ring_cycle[i], ring_cycle[(i + 1) % n]} == {a, b}:
            return i
    raise UnsupportedStructure("could not locate a ring-fusion bond on its own ring")


def _classify_shape(graph, atom_rings, ring_order, fusion_bonds_by_pair):
    """For each internal ring (fused to two neighbors), classify whether its
    two fusion bonds are opposite hexagon edges ('straight', e.g.
    naphthalene/anthracene-style linear fusion) or not ('bent', e.g.
    phenanthrene's angular fusion) -- see the module docstring."""
    shape = []
    for pos in range(1, len(ring_order) - 1):
        ring_idx = ring_order[pos]
        cycle = _ring_cycle(graph, list(atom_rings[ring_idx]))
        left_edge = _edge_index(cycle, *fusion_bonds_by_pair[frozenset((ring_order[pos - 1], ring_idx))])
        right_edge = _edge_index(cycle, *fusion_bonds_by_pair[frozenset((ring_idx, ring_order[pos + 1]))])
        diff = (left_edge - right_edge) % 6
        shape.append("straight" if diff == 3 else "bent")
    return tuple(shape)


def _periphery_cycle(mol, ring_atom_sets, fusion_bond_idxs):
    """The fused system's outer boundary: a single Hamiltonian cycle over
    every ring atom, obtained by dropping the fusion (chord) bonds from the
    ring skeleton. Valid because this module's scope (ortho-fusion only, no
    atom in 3+ rings, no branching) guarantees no interior atoms."""
    ring_atoms = set()
    for s in ring_atom_sets:
        ring_atoms |= s
    sub_adj = {a: [] for a in ring_atoms}
    for bond in mol.GetBonds():
        if bond.GetIdx() in fusion_bond_idxs:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ring_atoms and b in ring_atoms:
            sub_adj[a].append(b)
            sub_adj[b].append(a)
    if any(len(v) != 2 for v in sub_adj.values()):
        raise UnsupportedStructure(
            "the aromatic ring system's outer boundary is not a simple "
            "cycle (see P-25.3.1.3)"
        )
    start = next(iter(ring_atoms))
    cycle = [start]
    previous, current = None, start
    while True:
        next_atom = next(a for a in sub_adj[current] if a != previous)
        if next_atom == start:
            break
        cycle.append(next_atom)
        previous, current = current, next_atom
    return cycle


def _walk_locants(cycle, fusion_atoms, start_index, step):
    n = len(cycle)
    locants = {}
    counter = 1
    for k in range(n):
        atom = cycle[(start_index + step * k) % n]
        if atom in fusion_atoms:
            continue
        locants[atom] = counter
        counter += 1
    return locants


def _terminal_run_ends(cycle, terminal_atoms):
    """`terminal_atoms`: a terminal ring's 4 nonfusion atoms. They form one
    contiguous run on `cycle`; return its two end indices."""
    n = len(cycle)
    run_indices = sorted(i for i, a in enumerate(cycle) if a in terminal_atoms)
    if len(run_indices) != 4:
        raise UnsupportedStructure(
            "a terminal ring of the aromatic chain does not have exactly "
            "four nonfusion atoms (see P-25.3.1.3)"
        )
    breaks = [i for i in run_indices if (i - 1) % n not in run_indices]
    if len(breaks) != 1:
        raise UnsupportedStructure(
            "a terminal ring's nonfusion atoms are not contiguous on the "
            "periphery (see P-25.3.1.3)"
        )
    first = breaks[0]
    ordered = [(first + k) % n for k in range(4)]
    return ordered[0], ordered[-1]


def _straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order):
    """Naphthalene/tetracene/pentacene (P-25.3.3.1.1's general peripheral
    walk): position 1 must start at one of the two terminal rings, at
    whichever end leads *into* that ring's own run first (not straight back
    out to the fusion atom); each of the 2 terminal rings x 2 ends is a
    distinct valid numbering (the straight chain's own symmetry), searched
    below like every other module's lowest-locant tie-break."""
    cycle = _periphery_cycle(mol, ring_atom_sets, fusion_bond_idxs)
    membership = _atom_ring_membership(atom_rings)
    fusion_atoms = {atom for atom, rings in membership.items() if len(rings) >= 2}

    candidates = []
    for terminal_ring in (ring_order[0], ring_order[-1]):
        terminal_atoms = ring_atom_sets[terminal_ring] - fusion_atoms
        end_a, end_b = _terminal_run_ends(cycle, terminal_atoms)
        n = len(cycle)
        if (end_b - end_a) % n == 3:
            candidates.append(_walk_locants(cycle, fusion_atoms, end_a, 1))
            candidates.append(_walk_locants(cycle, fusion_atoms, end_b, -1))
        else:
            candidates.append(_walk_locants(cycle, fusion_atoms, end_b, 1))
            candidates.append(_walk_locants(cycle, fusion_atoms, end_a, -1))
    return candidates


def _anthracene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order):
    """Anthracene's traditional (non-walk) numbering -- see module
    docstring. r1/r3 are the two terminal rings, r2 the middle ring; x1/x2
    are r1's fusion pair, y1/y2 are r3's fusion pair. Both terminal-ring
    choices (r1<->r3) and both run directions within r1 are tried, matching
    the numbering's real 4-fold symmetry."""
    candidates = []
    for r1, r2, r3 in ((ring_order[0], ring_order[1], ring_order[2]), (ring_order[2], ring_order[1], ring_order[0])):
        x1, x2 = fusion_bonds_by_pair[frozenset((r1, r2))]
        y1, y2 = fusion_bonds_by_pair[frozenset((r2, r3))]
        run1_cycle = _ring_cycle(graph, list(ring_atom_sets[r1]))
        run3_cycle = _ring_cycle(graph, list(ring_atom_sets[r3]))
        for near_x, far_x in ((x1, x2), (x2, x1)):
            # near_x ("9a") is adjacent to position 1; far_x ("4a") to position 4.
            pos1, pos2, pos3, pos4 = _run_between(run1_cycle, near_x, far_x)
            meso9 = next(a for a in graph[near_x] if a in ring_atom_sets[r2] and a not in (x1, x2))
            far_y = next(a for a in graph[meso9] if a in ring_atom_sets[r3])  # "8a"
            near_y = y1 if far_y == y2 else y2  # "10a"
            pos8, pos7, pos6, pos5 = _run_between(run3_cycle, far_y, near_y)
            meso10 = next(a for a in graph[near_y] if a in ring_atom_sets[r2] and a != far_y)

            locants = {
                pos1: 1, pos2: 2, pos3: 3, pos4: 4,
                pos5: 5, pos6: 6, pos7: 7, pos8: 8,
                meso9: 9, meso10: 10,
            }
            candidates.append(locants)
    return candidates


def _phenanthrene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order):
    """Phenanthrene's traditional numbering coincides with a plain
    peripheral walk starting at the terminal-ring end that leads to the bay
    bond first -- see module docstring. Both terminal-ring choices give
    valid, symmetric numberings (phenanthrene's real 2-fold symmetry); the
    other 2 combinatorial (ring, end) choices would number the K-region
    atoms too early and are not traditionally valid, so they're not tried."""
    candidates = []
    for r1, r2, r3 in ((ring_order[0], ring_order[1], ring_order[2]), (ring_order[2], ring_order[1], ring_order[0])):
        x1, x2 = fusion_bonds_by_pair[frozenset((r1, r2))]
        y1, y2 = fusion_bonds_by_pair[frozenset((r2, r3))]
        bay = None
        for xa in (x1, x2):
            for ya in (y1, y2):
                if ya in graph[xa]:
                    bay = (xa, ya)
        if bay is None:
            raise UnsupportedStructure(
                "could not locate phenanthrene's bay-region fusion bond "
                "(see P-25.3.1.3)"
            )
        bay_x, bay_y = bay
        other_x = x2 if bay_x == x1 else x1
        other_y = y2 if bay_y == y1 else y1

        pos4, pos3, pos2, pos1 = _run_between(_ring_cycle(graph, list(ring_atom_sets[r1])), bay_x, other_x)
        pos5, pos6, pos7, pos8 = _run_between(_ring_cycle(graph, list(ring_atom_sets[r3])), bay_y, other_y)

        k1 = next(a for a in graph[other_y] if a in ring_atom_sets[r2] and a != bay_y)
        k2 = next(a for a in graph[k1] if a in ring_atom_sets[r2] and a != other_y)

        locants = {
            pos1: 1, pos2: 2, pos3: 3, pos4: 4,
            pos5: 5, pos6: 6, pos7: 7, pos8: 8,
            k1: 9, k2: 10,
        }
        candidates.append(locants)
    return candidates


def _benzene_candidates(graph, ring_atoms):
    cycle = _ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    candidates = []
    for start in range(n):
        rotated = cycle[start:] + cycle[:start]
        for seq in (rotated, list(reversed(rotated))):
            candidates.append({atom: i + 1 for i, atom in enumerate(seq)})
    return candidates


def _ring_substituents(graph, locants, ring_atoms, halogens, mol=None):
    """Like `_cyclic._substituents_for_ring`, but `locants` (nonfusion atom
    -> integer position) covers only part of the ring skeleton -- a fusion
    carbon has no free valence and so is never itself a locant -- while
    `ring_atoms` (the full skeleton, fusion included) is what a neighboring
    atom must be excluded from to count as a substituent branch root."""
    substituents = {}
    for atom, position in locants.items():
        branch_roots = [n for n in graph[atom] if n not in ring_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _candidate_key(parent, locants, ring_atoms, graph, halogens, omit_single_locant, stereo_display=None, mol=None):
    substituents = _ring_substituents(graph, locants, ring_atoms, halogens, mol=mol)
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    total_count = sum(len(info["locants"]) for info in grouped.values())
    if not grouped:
        name = parent
    elif omit_single_locant and total_count == 1:
        # P-14.3.3: on benzene, every ring position is equivalent before
        # substitution, so a single substituent's locant is not essential.
        if stereo_display is not None:
            display = stereo_display
        else:
            (only_name,) = grouped
            display = f"({only_name})" if grouped[only_name]["compound"] else only_name
        name = f"{display}{parent}"
    else:
        if stereo_display is not None:
            # Not benzene, so the ring attachment locant is significant and
            # still needs citing; only the substituent's own display text is
            # replaced by the stereo-decorated one (already validated to be
            # the system's sole substituent -- see `_stereo_display`).
            (only_name,) = grouped
            grouped = {stereo_display: {"locants": grouped[only_name]["locants"], "compound": False}}
        name = format_substituent_prefixes(grouped) + parent
    return locant_set, citation_locants, name


_RETAINED_NAMES = {
    (1, ()): "benzene",
    (2, ()): "naphthalene",
    (3, ("straight",)): "anthracene",
    (3, ("bent",)): "phenanthrene",
    (4, ("straight", "straight")): "tetracene",
    (5, ("straight", "straight", "straight")): "pentacene",
}


def _stereo_display(mol, graph, n, ring_atoms, halogens):
    """If `mol` has one or more specified tetrahedral stereocenters
    (P-92), build the bracketed "[(<locant><R/S>)-<name>]" substituent
    display P-91.3 requires when the stereocenter sits on a substituent
    branch rather than the ring itself (the Blue Book's own worked
    example, since an all-carbon aromatic ring atom is never itself a
    stereocenter). Returns None if there's no specified stereocenter at
    all (the caller proceeds exactly as before). Deliberately narrow: only a plain benzene
    or naphthalene ring (n in (1, 2)) with exactly one substituent,
    carrying exactly one specified stereocenter, is supported; anything
    else raises `UnsupportedStructure`. For naphthalene the ring
    attachment locant is still significant (unlike benzene) and is cited
    by the caller's normal substituent-prefix machinery -- this function
    only builds the substituent's own decorated display text."""
    stereo = specified_stereocenters(mol)
    if stereo is None:
        return None
    if n not in (1, 2):
        raise UnsupportedStructure(
            "a specified stereocenter combined with anything other than a "
            "plain benzene or naphthalene ring parent is not supported yet "
            "(see P-91.3)"
        )
    if len(stereo) != 1:
        raise UnsupportedStructure(
            "more than one specified stereocenter on a substituent branch "
            "is not supported yet (see P-91.3)"
        )
    stereo_atom, r_or_s = stereo[0]
    if stereo_atom in ring_atoms:
        raise UnsupportedStructure(
            "a specified stereocenter on the aromatic ring itself is not "
            "supported (an all-carbon aromatic ring atom can't be a "
            "genuine stereocenter)"
        )
    branch_attachments = [
        (ring_atom, neighbor)
        for ring_atom in ring_atoms
        for neighbor in graph[ring_atom]
        if neighbor not in ring_atoms
    ]
    if len(branch_attachments) != 1:
        raise UnsupportedStructure(
            "a specified stereocenter combined with anything other than "
            "exactly one substituent on the ring is not supported yet "
            "(see P-91.3)"
        )
    ring_atom, branch_root = branch_attachments[0]
    branch_name, branch_compound = name_branch(graph, branch_root, ring_atom, halogens, mol=mol)
    site_locant = branch_atom_locant(graph, branch_root, ring_atom, stereo_atom, halogens, mol=mol)
    descriptor = f"({site_locant}{r_or_s})-{branch_name}"
    return f"[{descriptor}]" if branch_compound else f"({descriptor})"


def name_aromatic_fused(mol, core) -> str:
    validate_atoms_and_bonds(mol)
    _validate_aromatic_bonds(mol)

    atom_rings, ring_atom_sets, fusion_bond_idxs = core
    n = len(atom_rings)

    membership = _atom_ring_membership(atom_rings)
    if any(len(rings) >= 3 for rings in membership.values()):
        raise UnsupportedStructure(
            "peri-fused aromatic ring systems (an atom shared by three or "
            "more rings, e.g. pyrene, acenaphthylene) are not supported "
            "(see P-25.3.1.3, ortho-fusion only)"
        )

    adj, fusion_bonds_by_pair = _ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
    ring_order = _ring_path_order(adj, n)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set()
    for s in ring_atom_sets:
        ring_atoms |= s

    stereo_display = _stereo_display(mol, graph, n, ring_atoms, halogens)

    if n == 1:
        candidates = _benzene_candidates(graph, ring_atom_sets[0])
        parent = "benzene"
        omit_single_locant = True
    else:
        shape = _classify_shape(graph, atom_rings, ring_order, fusion_bonds_by_pair)
        parent = _RETAINED_NAMES.get((n, shape))
        if parent is None:
            raise UnsupportedStructure(
                "this ortho-fused aromatic ring chain does not have one of "
                "the retained names supported so far (benzene, naphthalene, "
                "anthracene, phenanthrene, tetracene, pentacene); longer or "
                "angular (polyphene, n>=4) chains require systematic "
                "benzo[x,y-z]fusion[...] name construction, which is out of "
                "scope for now (see P-25.3.1-P-25.3.3)"
            )
        omit_single_locant = False
        if parent == "anthracene":
            candidates = _anthracene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order)
        elif parent == "phenanthrene":
            candidates = _phenanthrene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order)
        else:
            candidates = _straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order)

    best_key = None
    best_name = None
    for locants in candidates:
        key = _candidate_key(parent, locants, ring_atoms, graph, halogens, omit_single_locant, stereo_display, mol=mol)
        if best_key is None or key < best_key:
            best_key, best_name = key, key[-1]
    return best_name
