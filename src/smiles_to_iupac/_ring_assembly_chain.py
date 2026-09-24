"""Naming of unbranched ring assemblies of 3-6 identical benzene rings
(terphenyl/quaterphenyl/quinquephenyl/sexiphenyl), identical monocyclic
saturated all-carbon rings (tercyclopropane, tercyclohexane, ...), or
identical pyridine rings (terpyridine, quaterpyridine, ...), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-28.3.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): an
  unbranched chain of N (>=3) identical cyclic parent hydrides, each
  consecutive pair joined by a single bond, is numbered with a new
  *composite-locant* scheme distinct from `_ring_assembly.py`'s own
  primed-locant numbering for the N=2 case: each ring gets a **primary
  locant** (1, 2, 3, ... for its position along the chain) and its own
  internal ring-atom locants are cited as **superscripts** on that primary
  locant (rendered here, matching the primary source's own plain-text
  fallback, as the two digits concatenated -- ring 2's atom 4 is "24").
  Locants indicating the junction (attachment) points are cited at the
  front of the name, each junction's pair of locants comma-separated, and
  junctions colon-separated in path order; the choice of which physical
  ring is "ring 1" and each ring's own internal numbering direction is
  made to give the lowest locant set to the set of all junction points
  taken together, then lowest order of citation if a tie remains --
  confirmed PIN worked examples `11,21:24,31-terphenyl` (p-terphenyl) and
  the tie-break pair `11,21:22,31:33,41-quatercyclobutane` (not
  `11,21:23,31:32,41-...`; `tmp/bluebook/P2.txt` ~7904-7915).
- P-28.3.2: for benzene rings specifically (not any other ring), the
  retained substituent name 'phenyl' is used with the Latin multiplying
  prefix (ter/quater/quinque/sexi, P-14.2.3 -- this project's existing
  `_numerals.py` only has the basic/Greek series di/tri/tetra/..., used
  for suffix/substituent counts, not this distinct Latin series, so a
  small local table is used here instead of extending that shared one)
  instead of the general 'ter'+parent-hydride-name construction --
  confirmed PIN worked examples `11,21:24,31-terphenyl` and
  `11,21:23,31:33,41-quaterphenyl` (`tmp/bluebook/P2.txt` ~7970-7978).
- P-28.3.1's *general* rule applies as-is to any other identical cyclic
  parent hydride: the Latin multiplying prefix directly in front of the
  plain parent-hydride name (not a substituent-group name) -- confirmed
  PIN worked examples `11,21:22,31-tercyclopropane` (not
  "tercyclopropyl") and `12,25:22,34-terpyridine`, `tmp/bluebook/P2.txt`
  ~7904, ~7933-7938.
- **Numbering for a ring whose own numbering isn't free** (pyridine):
  P-28.2.1's text, shared by P-28.2.1/P-28.2.2/P-28.3 alike, states "Each
  cyclic system is numbered in the traditional way... Lowest possible
  locants must be used to denote the positions of attachment" -- for a
  symmetric ring (benzo/cycloalkane) "traditional way" leaves every
  rotation/direction equally valid, which is what `_ring_numberings`
  already searches over. Pyridine's own numbering isn't free: its
  nitrogen is always locant 1 (confirmed by the worked example
  `2,2'-bipyridine`, `tmp/bluebook/P2.txt` ~7799), so only the
  *direction* from that fixed point is chosen, by the same
  lowest-junction-locant rule -- `name_ring_assembly_chain` reuses
  `_ring_numberings` unchanged for this too, just starting from each
  ring's own nitrogen atom instead of its junction atom(s).
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the
  same way as in every other ring module, cited under the same
  composite-locant scheme; a substituent locant only breaks a tie left
  after the junction-locant-set rule above is already satisfied (P-28.3.1
  gives the junction points numbering priority, not substituents).

Scope: an unbranched chain of 3-6 disjoint, identical (same kind, same
size) rings -- 6-membered all-carbon aromatic (benzo), monocyclic
saturated all-carbon (any one ring size, e.g. all cyclopropane or all
cyclohexane, not mixed), or 6-membered pyridine -- each consecutive pair
joined by exactly one single (non-aromatic) bond and no other inter-ring
bond (a branched or cyclic ring-assembly topology -- P-28.5/P-28.6 -- is
out of scope and must fall through to `UnsupportedStructure` elsewhere,
not be misnamed), each ring optionally bearing halogen substituents.
Explicitly out of scope (raise `UnsupportedStructure` via the generic
fallback in `core.py`, since `find_ring_assembly_chain_core` below simply
returns None for any of these): N=2 (stays `_ring_assembly.py`'s own
job, and doesn't yet handle pyridine either -- a separate gap), N>6, any
heteroaromatic ring other than pyridine (furan/thiophene/pyrrole/... --
each has its own role-sequence/symmetry and is its own follow-up), mixed
ring kinds/sizes, any branched/cyclic ring-assembly topology, and
indicated hydrogen (P-28.2.3 -- not reachable by any ring kind here
anyway).
"""

from itertools import product

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    ring_cycle,
    substituent_locant_set_and_citation,
    validate_allowed_atoms,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_MIN_RINGS, _MAX_RINGS = 3, 6
_MULTIPLIER = {3: "ter", 4: "quater", 5: "quinque", 6: "sexi"}


def _ring_kind(mol, ring):
    """("aromatic", 6) for an all-carbon benzo ring, ("pyridine", 6) for a
    mancude 6-ring with exactly one ring nitrogen (no H, no charge -- the
    pyridine-type "no spare valence" nitrogen, P-25.2.1) and five ring
    carbons (each bearing 0 or 1 H, degree 3 only when it's a junction/
    substituent-bearing atom), ("saturated", n) for an n-membered
    monocyclic all-carbon ring with only single ring bonds, else None."""
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if len(ring) == 6 and all(a.GetIsAromatic() for a in atoms):
        atomic_nums = [a.GetAtomicNum() for a in atoms]
        if all(n == 6 for n in atomic_nums):
            return "aromatic", 6
        nitrogens = [a for a in atoms if a.GetAtomicNum() == 7]
        carbons = [a for a in atoms if a.GetAtomicNum() == 6]
        if len(nitrogens) != 1 or len(carbons) != 5:
            return None
        (nitrogen,) = nitrogens
        if nitrogen.GetTotalNumHs() != 0 or nitrogen.GetFormalCharge() != 0 or nitrogen.GetIsotope() != 0:
            return None
        if any(c.GetTotalNumHs() not in (0, 1) or c.GetFormalCharge() != 0 for c in carbons):
            return None
        return "pyridine", 6
    if any(a.GetAtomicNum() != 6 or a.GetIsAromatic() for a in atoms):
        return None
    ring_set = set(ring)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ring_set and b in ring_set and bond.GetBondTypeAsDouble() != 1.0:
            return None
    return "saturated", len(ring)


def find_ring_assembly_chain_core(mol):
    """Return (path, connections, ring_kind) if `mol` is an unbranched
    chain of 3-6 disjoint identical-kind-and-size rings (all benzo, or all
    one saturated ring size), consecutive rings joined by exactly one
    single non-aromatic bond and no other inter-ring bond, else None.
    `path` is a list of ring atom-index tuples in one arbitrary path order
    (either end may become "ring 1" -- both are tried when numbering);
    `connections` is a list of (attach_in_ring_i, attach_in_ring_i+1)
    atom-idx pairs, one per consecutive pair along `path`; `ring_kind` is
    `_ring_kind`'s own ("aromatic", 6) or ("saturated", n)."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    n = len(atom_rings)
    if not (_MIN_RINGS <= n <= _MAX_RINGS):
        return None
    kinds = {_ring_kind(mol, ring) for ring in atom_rings}
    if len(kinds) != 1 or None in kinds:
        return None
    (ring_kind,) = kinds

    ring_sets = [set(r) for r in atom_rings]
    for i in range(n):
        for j in range(i + 1, n):
            if ring_sets[i] & ring_sets[j]:
                return None

    ring_of_atom = {}
    for idx, rs in enumerate(ring_sets):
        for a in rs:
            ring_of_atom[a] = idx

    ring_adj = {i: [] for i in range(n)}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ra, rb = ring_of_atom.get(a), ring_of_atom.get(b)
        if ra is None or rb is None or ra == rb:
            continue
        if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() != 1.0:
            return None
        ring_adj[ra].append((rb, a, b))
        ring_adj[rb].append((ra, b, a))

    if sum(len(v) for v in ring_adj.values()) != 2 * (n - 1):
        return None
    if any(len(v) > 2 for v in ring_adj.values()):
        return None
    if any(len({x[0] for x in v}) != len(v) for v in ring_adj.values()):
        return None

    endpoints = [i for i in range(n) if len(ring_adj[i]) == 1]
    if len(endpoints) != 2:
        return None

    order = [endpoints[0]]
    prev = None
    while len(order) < n:
        current = order[-1]
        nxt = [x for x in ring_adj[current] if x[0] != prev]
        if len(nxt) != 1:
            return None
        order.append(nxt[0][0])
        prev = current
    if order[-1] != endpoints[1]:
        return None

    connections = []
    for i in range(n - 1):
        ra, rb = order[i], order[i + 1]
        entry = next(x for x in ring_adj[ra] if x[0] == rb)
        connections.append((entry[1], entry[2]))

    path = [atom_rings[idx] for idx in order]
    return path, connections, ring_kind


def _ring_numberings(graph, ring_atoms, attach_atoms):
    """Yield {atom: local_position (1..len(ring_atoms))} for every valid
    numbering of one ring: a terminal ring (`attach_atoms` has 1 atom)
    always starts at that atom (2 direction choices); a middle ring
    (`attach_atoms` has 2) may start at either one (2 x 2 = 4 choices)."""
    cycle = ring_cycle(graph, list(ring_atoms))
    for start in attach_atoms:
        idx = cycle.index(start)
        rotated = cycle[idx:] + cycle[:idx]
        for seq in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
            yield {atom: position for position, atom in enumerate(seq, start=1)}


def name_ring_assembly_chain(mol, core) -> str:
    path, connections, ring_kind = core
    n = len(path)
    kind, ring_size = ring_kind

    ring_atoms_all = set()
    for ring in path:
        ring_atoms_all |= set(ring)

    if kind == "pyridine":
        nitrogen_idxs = {idx for idx in ring_atoms_all if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7}
        validate_allowed_atoms(
            mol,
            "heteroatoms other than the pyridine rings' own nitrogens and a "
            "halogen substituent are not supported yet",
            [
                (
                    7,
                    nitrogen_idxs,
                    "a nitrogen atom outside the pyridine rings' own is not supported yet",
                )
            ],
            aromatic_ring_atoms=ring_atoms_all,
        )
        ring_word = "pyridine"
    else:
        validate_atoms_and_bonds(mol)
        ring_word = "phenyl" if kind == "aromatic" else "cyclo" + alkane_name(ring_size)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms_all or idx in halogens:
            continue
        raise UnsupportedStructure(
            "a substituent other than a halogen is not supported yet by this module"
        )

    attach_atoms_per_ring = [set() for _ in range(n)]
    for i, (a, b) in enumerate(connections):
        attach_atoms_per_ring[i].add(a)
        attach_atoms_per_ring[i + 1].add(b)

    best_key = None
    best_name = None
    for order, conns in ((path, connections), (list(reversed(path)), None)):
        if conns is None:
            # Reversed path: recompute each ring's attach atoms and the
            # junction pairs in the new (reversed) path order.
            rev_attach = [set() for _ in range(n)]
            rev_connections = []
            for i in range(n - 1):
                a, b = connections[n - 2 - i]
                rev_connections.append((b, a))
                rev_attach[i].add(b)
                rev_attach[i + 1].add(a)
            attach_sets = rev_attach
            conns = rev_connections
        else:
            attach_sets = attach_atoms_per_ring

        per_ring_candidates = []
        for i in range(n):
            if kind == "pyridine":
                # P-28.2.1's "each cyclic system is numbered in the
                # traditional way" -- pyridine's own numbering always
                # fixes its nitrogen at locant 1 (confirmed by the
                # primary source's own '2,2'-bipyridine' worked example);
                # only the direction (which of `_ring_numberings`' two
                # candidates starting there) is chosen freely, by the
                # same lowest-junction-locant rule as every other kind.
                (nitrogen,) = [idx for idx in order[i] if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
                starts = {nitrogen}
            else:
                starts = attach_sets[i]
            per_ring_candidates.append(list(_ring_numberings(graph, order[i], starts)))
        for combo in product(*per_ring_candidates):
            locants = {}
            for i in range(n):
                for atom, pos in combo[i].items():
                    locants[atom] = (i + 1, pos)

            junction_pairs = []
            for i, (a, b) in enumerate(conns):
                junction_pairs.append((locants[a], locants[b]))
            all_junction_locants = [loc for pair in junction_pairs for loc in pair]
            junction_locant_set = tuple(sorted(all_junction_locants))
            junction_citation = tuple(all_junction_locants)

            substituents = {}
            for atom, composite in locants.items():
                branch_roots = [nb for nb in graph[atom] if nb not in ring_atoms_all]
                if branch_roots:
                    position = f"{composite[0]}{composite[1]}"
                    substituents[position] = [
                        name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots
                    ]
            grouped = group_substituents(substituents)
            sub_locant_set, _, sub_citation = substituent_locant_set_and_citation(grouped)
            prefix = format_substituent_prefixes(grouped) if grouped else ""

            junction_str = ":".join(
                f"{a[0]}{a[1]},{b[0]}{b[1]}" for a, b in junction_pairs
            )
            base = f"{junction_str}-{_MULTIPLIER[n]}{ring_word}"
            name = base if not prefix else f"{prefix}-{base}"

            key = (junction_locant_set, junction_citation, sub_locant_set, sub_citation, name)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    return best_name
