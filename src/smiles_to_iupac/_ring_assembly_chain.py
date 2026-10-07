"""Naming of unbranched ring assemblies of 3-6 identical benzene rings
(terphenyl/quaterphenyl/quinquephenyl/sexiphenyl), identical monocyclic
saturated all-carbon rings (tercyclopropane, tercyclohexane, ...), or
identical mancude heteroaromatic rings whose own numbering is fixed by a
role sequence with no N-H tautomer ambiguity and whose own parent-hydride
name doesn't start with a locant -- pyridine (terpyridine,
quaterpyridine, ...) and, since this module now reuses
`_hetero_monocyclic.py`'s own `_ROLE_SEQUENCES` table directly instead of
special-casing pyridine alone, furan, thiophene, selenophene,
tellurophene, pyridazine, pyrimidine, and pyrazine too (terthiophene,
terfuran, terpyrimidine, ... e.g. -- see `_NON_NH_ROLE_SEQUENCES` for why
the Hantzsch-Widman oxazole/thiazole/selenazole/tellurazole parents are
excluded for now), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-28.3.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): an
  unbranched chain of N (>=3) identical cyclic parent hydrides, each
  consecutive pair joined by a single bond, is numbered with a new
  *composite-locant* scheme distinct from `_ring_assembly.py`'s own
  primed-locant numbering for the N=2 case: each ring gets a **primary
  locant** (1, 2, 3, ... for its position along the chain) and its own
  internal ring-atom locants are cited as **superscripts** on that primary
  locant (ring 2's atom 4 is "2⁴").
  Locants indicating the junction (attachment) points are cited at the
  front of the name, each junction's pair of locants comma-separated, and
  junctions colon-separated in path order; the choice of which physical
  ring is "ring 1" and each ring's own internal numbering direction is
  made to give the lowest locant set to the set of all junction points
  taken together, then lowest order of citation if a tie remains --
  confirmed PIN worked examples `11,21:24,31-terphenyl` (p-terphenyl) and
  the tie-break pair `11,21:22,31:33,41-quatercyclobutane` (not
  `11,21:23,31:32,41-...`; the Blue Book).
- P-28.3.2: for benzene rings specifically (not any other ring), the
  retained substituent name 'phenyl' is used with the Latin multiplying
  prefix (ter/quater/quinque/sexi, P-14.2.3 -- this project's existing
  `_numerals.py` only has the basic/Greek series di/tri/tetra/..., used
  for suffix/substituent counts, not this distinct Latin series, so a
  small local table is used here instead of extending that shared one)
  instead of the general 'ter'+parent-hydride-name construction --
  confirmed PIN worked examples `11,21:24,31-terphenyl` and
  `11,21:23,31:33,41-quaterphenyl` (the Blue Book).
- P-28.3.1's *general* rule applies as-is to any other identical cyclic
  parent hydride: the Latin multiplying prefix directly in front of the
  plain parent-hydride name (not a substituent-group name) -- confirmed
  PIN worked examples `11,21:22,31-tercyclopropane` (not
  "tercyclopropyl") and `12,25:22,34-terpyridine`, the Blue Book.
- **Numbering for a ring whose own numbering isn't free** (any
  `_ROLE_SEQUENCES` parent other than the 3 N-H tautomer-ambiguous ones):
  P-28.2.1's text, shared by P-28.2.1/P-28.2.2/P-28.3 alike, states "Each
  cyclic system is numbered in the traditional way... Lowest possible
  locants must be used to denote the positions of attachment" -- for a
  symmetric ring (benzo/cycloalkane) "traditional way" leaves every
  rotation/direction equally valid, which is what `_ring_numberings`
  already searches over. A heteroaromatic parent's own numbering isn't
  free: its heteroatom(s) sit at fixed role-sequence positions (confirmed
  for pyridine by the worked example `2,2'-bipyridine`,
  the Blue Book, and the same "traditional way" text
  applies identically to every other role-sequence parent), so only the
  alignments (rotation + direction) whose element pattern actually
  matches that parent's role sequence are considered --
  `_hetero_ring_alignments` reuses `_hetero_monocyclic.py`'s own
  `_ring_alignments`/`_ROLE_SEQUENCES` matching machinery for this,
  rather than a per-parent special case, and `name_ring_assembly_chain`
  picks among the resulting candidates by the same lowest-junction-locant
  rule as every other ring kind. Excludes 1H-imidazole/1H-pyrazole (real
  N-H prototropic-tautomer ambiguity -- which ring nitrogen is "N1" isn't
  fixed by structure alone -- needs extra disambiguation this module
  doesn't do yet); 1H-pyrrole is handled separately below.
- P-28.2.3 (indicated hydrogen, extended to N>=3 chains): 1H-pyrrole is
  the one N-H tautomer-ambiguous parent supported here (see
  `_ring_assembly.py`'s own N=2 case, #949, for why imidazole/pyrazole
  stay excluded -- a second "pyridine-type" nitrogen with no spare
  valence of its own, unlike every one of pyrrole's five positions).
  Confirmed PIN worked examples `11H,33H-13,23:23,33-terindole` and
  `12H,22H,24H,32H-11,21:23,31-terpyrimidine` (the Blue Book) establish the composite-locant indicated-hydrogen citation
  format for N>=3 chains: front-of-name, comma-separated per occupied
  ring, each using that ring's own composite locant. For pyrrole
  specifically this reduces to the same per-ring check #949 already
  shipped for N=2 (every one of pyrrole's positions -- including its own
  N -- is a direct 1-for-1 substitution with no double-bond
  recalculation, unlike the pyrimidine case above): a ring needs
  indicated hydrogen at its own N-derived composite locant (always that
  ring's "...1", since 1H-pyrrole's N is fixed at role-sequence position
  1) if and only if none of that ring's own junction bonds sit at its own
  N.
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the
  same way as in every other ring module, cited under the same
  composite-locant scheme; a substituent locant only breaks a tie left
  after the junction-locant-set rule above is already satisfied (P-28.3.1
  gives the junction points numbering priority, not substituents).

Scope: an unbranched chain of 3-6 disjoint, identical (same kind, same
size) rings -- 6-membered all-carbon aromatic (benzo), monocyclic
saturated all-carbon (any one ring size, e.g. all cyclopropane or all
cyclohexane, not mixed), any mancude 5- or 6-membered ring matching a
`_NON_NH_ROLE_SEQUENCES` parent (pyridine, furan, thiophene, selenophene,
tellurophene, pyridazine, pyrimidine, pyrazine), or 1H-pyrrole -- each
consecutive pair joined by exactly one single (non-aromatic) bond and no
other inter-ring bond (a branched or cyclic ring-assembly topology --
P-28.5/P-28.6 -- is out of scope and must fall through to
`UnsupportedStructure` elsewhere, not be misnamed), each ring optionally
bearing halogen substituents.
Explicitly out of scope (raise `UnsupportedStructure` via the generic
fallback in `core.py`, since `find_ring_assembly_chain_core` below simply
returns None for any of these): N=2 (stays `_ring_assembly.py`'s own
job), N>6, 1H-imidazole/1H-pyrazole (real N-H tautomer ambiguity -- its
own follow-up), the locant-prefixed Hantzsch-Widman parents (1,3-/1,2-
oxazole/thiazole/selenazole/tellurazole -- unconfirmed composite-name
formatting, see `_NON_NH_ROLE_SEQUENCES`), mixed ring kinds/sizes, and
any branched/cyclic ring-assembly topology.
"""

from itertools import product

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    ring_cycle,
    superscript_locant,
    substituent_locant_set_and_citation,
    validate_allowed_atoms,
    validate_atoms_and_bonds,
)
from ._hetero_monocyclic import _ROLE_SEQUENCES, _TAUTOMER_AMBIGUOUS_UNLESS_N1, _ring_alignments
from ._numerals import alkane_name, numerical_term
from ._substituents import format_substituent_prefixes, name_branch

_MIN_RINGS, _MAX_RINGS = 3, 6
_MULTIPLIER = {3: "ter", 4: "quater", 5: "quinque", 6: "sexi"}

_NON_NH_ROLE_SEQUENCES = {
    name: seq
    for name, seq in _ROLE_SEQUENCES.items()
    # Excludes the N-H tautomer-ambiguous parents (pyrrole/imidazole/
    # pyrazole, see the module docstring) and, separately, the 8
    # Hantzsch-Widman parents whose own name starts with a locant
    # (1,3-/1,2-oxazole/thiazole/selenazole/tellurazole) -- P-28.3.1's own
    # worked examples only ever concatenate 'ter'/'quater'/... directly in
    # front of a plain-word parent name (e.g. 'terpyrimidine'), and this
    # project has no confirmed primary-source example of how that
    # concatenation is meant to read when the parent name itself starts
    # with a locant (bare 'ter1,3-thiazole' is ambiguous; general IUPAC
    # practice elsewhere in this same chapter encloses such a name in
    # square brackets when composing it into a larger name, e.g.
    # 'spiroter[[1,3,2]benzodioxathiole]', the Blue Book --
    # but that is P-24's dispiro/spiro construction, not P-28.3's, so it
    # is not assumed to carry over here without its own worked example).
    if name not in _TAUTOMER_AMBIGUOUS_UNLESS_N1 and not name[0].isdigit()
}


def _match_hetero_ring_parent(mol, graph, ring):
    """The `_ROLE_SEQUENCES` parent name (e.g. "pyridine", "thiophene",
    "pyrazine") this aromatic ring's element/H-count pattern matches, or
    None if it matches none of them -- reuses `_hetero_monocyclic.py`'s
    own role-sequence table and `_ring_alignments` (every rotation +
    direction of the ring's actual atom order) instead of a per-parent
    special case. Excludes the 3 N-H tautomer-ambiguous parents (pyrrole/
    imidazole/pyrazole -- see the module docstring)."""
    ring_order = ring_cycle(graph, list(ring))
    elements = {atom: mol.GetAtomWithIdx(atom).GetSymbol() for atom in ring}
    h_counts = {atom: mol.GetAtomWithIdx(atom).GetTotalNumHs() for atom in ring}
    if any(mol.GetAtomWithIdx(atom).GetFormalCharge() != 0 or mol.GetAtomWithIdx(atom).GetIsotope() != 0 for atom in ring):
        return None
    for parent_name, role_sequence in _NON_NH_ROLE_SEQUENCES.items():
        if len(role_sequence) != len(ring):
            continue
        for candidate in _ring_alignments(ring_order):
            valid = True
            for position, atom in enumerate(candidate, start=1):
                role_element, role_has_h = role_sequence[position - 1]
                if elements[atom] != role_element:
                    valid = False
                    break
                if role_has_h:
                    if h_counts[atom] not in (0, 1):
                        valid = False
                        break
                elif h_counts[atom] != 0:
                    valid = False
                    break
            if valid:
                return parent_name
    return None


def _hetero_ring_alignments(mol, graph, ring, parent_name):
    """Every valid {atom: local_position} numbering of one ring consistent
    with `parent_name`'s own fixed role sequence (P-28.2.1's "traditional
    numbering" applied to a ring whose heteroatom(s) fix its own
    numbering) -- unlike `_match_hetero_ring_parent`, returns every
    matching alignment (there can be more than one when the role
    sequence has its own rotational/reflective symmetry, e.g. pyrazine's
    N,C,C,N,C,C), since the caller picks among them by the whole
    assembly's junction-locant rule, not a per-ring choice."""
    ring_order = ring_cycle(graph, list(ring))
    role_sequence = _ROLE_SEQUENCES[parent_name]
    if len(role_sequence) != len(ring_order):
        return
    elements = {atom: mol.GetAtomWithIdx(atom).GetSymbol() for atom in ring}
    for candidate in _ring_alignments(ring_order):
        position_of = {atom: position for position, atom in enumerate(candidate, start=1)}
        if all(elements[atom] == role_sequence[position - 1][0] for atom, position in position_of.items()):
            yield position_of


def saturated_counterpart_kind(mol, graph, ring):
    """(parent_name, size) of the mancude parent whose fully saturated form `ring` is (piperidine -> pyridine,
    oxolane -> furan), for the 5- and 6-membered `_NON_NH_ROLE_SEQUENCES` parents, else None."""
    if len(ring) not in (5, 6):
        return None
    ring_set = set(ring)
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
        return None
    if any(
        b.GetBondTypeAsDouble() != 1.0
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
    ):
        return None
    if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ring):
        return None
    for parent_name, sequence in _NON_NH_ROLE_SEQUENCES.items():
        if len(sequence) == len(ring) and any(True for _ in _hetero_ring_alignments(mol, graph, ring, parent_name)):
            return parent_name, len(ring)
    return None


def hydro_sort_key(locant):
    """Order of a (possibly primed) locant: 1, 1', 2, 2' ..."""
    return (int(locant.rstrip("'")), locant.endswith("'"))


def hydro_locants(positions):
    """Sorted hydro locants for a saturated ring of a mancude assembly: every position of a six-membered ring,
    every one but the chalcogen (position 1) of a five-membered ring (P-54.3, P-31.1.4.2.4)."""
    positions = list(positions)
    if len(positions) == 5:
        positions = [p for p in positions if p.rstrip("'") != "1"]
    return sorted(positions, key=hydro_sort_key)


def hydro_prefix(locants):
    return f"{','.join(locants)}-{numerical_term(len(locants))}hydro-" if locants else ""


def _ring_kind(mol, ring):
    """("aromatic", 6) for an all-carbon benzo ring, (parent_name, size)
    for a mancude 5- or 6-ring matching one of `_ROLE_SEQUENCES`'s non-NH
    parents (e.g. ("pyridine", 6), ("thiophene", 5), ("pyrazine", 6)),
    ("saturated", n) for an n-membered monocyclic all-carbon ring with
    only single ring bonds, else None."""
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if len(ring) in (5, 6) and all(a.GetIsAromatic() for a in atoms):
        if len(ring) == 6 and all(a.GetAtomicNum() == 6 for a in atoms):
            return "aromatic", 6
        graph = adjacency(mol)
        parent_name = _match_hetero_ring_parent(mol, graph, ring)
        return (parent_name, len(ring)) if parent_name is not None else None
    if any(a.GetAtomicNum() != 6 or a.GetIsAromatic() for a in atoms):
        return None
    ring_set = set(ring)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ring_set and b in ring_set and bond.GetBondTypeAsDouble() != 1.0:
            return None
    return "saturated", len(ring)


def _pyrrole_ring_kind(mol, ring):
    """("1H-pyrrole", 5) if `ring` is a mancude 5-membered ring with one
    N and four C (element-only match against `_ROLE_SEQUENCES`'s own
    "1H-pyrrole" entry, via `_hetero_ring_alignments`), else None --
    `_ring_kind` itself excludes 1H-pyrrole (it's not in
    `_NON_NH_ROLE_SEQUENCES`, see the module docstring), so this is
    checked separately rather than widening that shared table (which
    would also, wrongly, enable 1H-imidazole/1H-pyrazole-style
    complications this module doesn't handle). Shared by this module's
    own 3-6-ring case and `_ring_assembly.py`'s N=2 case."""
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if len(ring) != 5 or not all(a.GetIsAromatic() for a in atoms):
        return None
    graph = adjacency(mol)
    if any(True for _ in _hetero_ring_alignments(mol, graph, ring, "1H-pyrrole")):
        return "1H-pyrrole", 5
    return None


def find_ring_assembly_chain_core(mol):
    """Return (path, connections, ring_kind) if `mol` is an unbranched
    chain of 3-6 disjoint identical-kind-and-size rings (all benzo, or all
    one saturated ring size), consecutive rings joined by exactly one
    single non-aromatic bond and no other inter-ring bond, else None.
    `path` is a list of ring atom-index tuples in one arbitrary path order
    (either end may become "ring 1" -- both are tried when numbering);
    `connections` is a list of (attach_in_ring_i, attach_in_ring_i+1)
    atom-idx pairs, one per consecutive pair along `path`; `ring_kind` is
    `_ring_kind`'s own ("aromatic", 6) or ("saturated", n), or
    `_pyrrole_ring_kind`'s own ("1H-pyrrole", 5)."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    n = len(atom_rings)
    if not (_MIN_RINGS <= n <= _MAX_RINGS):
        return None
    graph = adjacency(mol)
    kinds, hydro_rings = [], []
    for ring in atom_rings:
        kind = _ring_kind(mol, ring) or _pyrrole_ring_kind(mol, ring)
        if kind is None:
            kind = saturated_counterpart_kind(mol, graph, ring)
            if kind is not None:
                hydro_rings.append(frozenset(ring))
        kinds.append(kind)
    if ("aromatic", 6) in kinds:
        for i, kind in enumerate(kinds):
            if kind == ("saturated", 6):
                kinds[i] = ("aromatic", 6)
                hydro_rings.append(frozenset(atom_rings[i]))
    if len(set(kinds)) != 1 or None in kinds or len(hydro_rings) == n:
        return None
    (ring_kind,) = set(kinds)

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
    return path, connections, ring_kind, hydro_rings


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


def validate_hetero_ring_assembly_atoms(mol, ring_atoms_all, kind):
    """Reject any heteroatom outside `ring_atoms_all`'s own `kind`-role
    heteroatoms (plus a halogen substituent, per `validate_allowed_atoms`)
    -- shared by this module's own 3-6-ring case and `_ring_assembly.py`'s
    N=2 case for a `_NON_NH_ROLE_SEQUENCES` ring kind."""
    heteroatom_idxs = {idx for idx in ring_atoms_all if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6}
    atomic_nums_present = {mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in heteroatom_idxs}
    validate_allowed_atoms(
        mol,
        f"heteroatoms other than the {kind} rings' own heteroatoms and a "
        "halogen substituent are not supported yet",
        [
            (
                atomic_num,
                {idx for idx in heteroatom_idxs if mol.GetAtomWithIdx(idx).GetAtomicNum() == atomic_num},
                f"a heteroatom outside the {kind} rings' own is not supported yet",
            )
            for atomic_num in atomic_nums_present
        ],
        aromatic_ring_atoms=ring_atoms_all,
    )


def name_ring_assembly_chain(mol, core) -> str:
    path, connections, ring_kind, hydro_rings = core
    n = len(path)
    kind, ring_size = ring_kind

    ring_atoms_all = set()
    for ring in path:
        ring_atoms_all |= set(ring)

    if kind in _NON_NH_ROLE_SEQUENCES:
        validate_hetero_ring_assembly_atoms(mol, ring_atoms_all, kind)
        ring_word = kind
    elif kind == "1H-pyrrole":
        validate_hetero_ring_assembly_atoms(mol, ring_atoms_all, kind)
        ring_word = "pyrrole"
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

        indicated_hydrogen_needed = []
        if kind == "1H-pyrrole":
            # Structural (not numbering-dependent, see `_ring_assembly.py`'s
            # own N=2 case): a ring needs indicated hydrogen at its own
            # N-derived composite locant iff none of its own junction bonds
            # sit at its own N -- 1H-pyrrole's N is always role-sequence
            # position 1 (`_pyrrole_ring_kind`'s docstring), so that locant
            # is always "<ring-number>1".
            for i in range(n):
                n_atom = next(atom for atom in order[i] if mol.GetAtomWithIdx(atom).GetAtomicNum() == 7)
                if n_atom not in attach_sets[i]:
                    indicated_hydrogen_needed.append(f"{i + 1}1H")
        indicated_hydrogen_prefix = ",".join(indicated_hydrogen_needed) + "-" if indicated_hydrogen_needed else ""

        per_ring_candidates = []
        for i in range(n):
            if kind in _NON_NH_ROLE_SEQUENCES or kind == "1H-pyrrole":
                # P-28.2.1's "each cyclic system is numbered in the
                # traditional way" -- a heteroaromatic parent's own
                # numbering always fixes its heteroatom(s) at their role-
                # sequence positions (confirmed for pyridine by the
                # primary source's own '2,2'-bipyridine' worked example,
                # and the same text applies identically to every other
                # role-sequence parent); only the choice among the
                # resulting (possibly several, if the role sequence has
                # its own symmetry) alignments is made freely, by the
                # same lowest-junction-locant rule as every other kind.
                per_ring_candidates.append(list(_hetero_ring_alignments(mol, graph, order[i], kind)))
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
            hydro = sorted(
                (i + 1, pos)
                for i in range(n)
                if frozenset(order[i]) in hydro_rings
                for pos in combo[i].values()
                if not (ring_size == 5 and pos == 1)
            )
            all_junction_locants = [loc for pair in junction_pairs for loc in pair]
            junction_locant_set = tuple(sorted(all_junction_locants))
            junction_citation = tuple(all_junction_locants)

            substituents = {}
            for atom, composite in locants.items():
                branch_roots = [nb for nb in graph[atom] if nb not in ring_atoms_all]
                if branch_roots:
                    position = superscript_locant(*composite)
                    substituents[position] = [
                        name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots
                    ]
            grouped = group_substituents(substituents)
            sub_locant_set, _, sub_citation = substituent_locant_set_and_citation(grouped)
            prefix = format_substituent_prefixes(grouped) if grouped else ""

            junction_str = ":".join(
                f"{superscript_locant(*a)},{superscript_locant(*b)}" for a, b in junction_pairs
            )
            hydro_str = hydro_prefix([superscript_locant(*loc) for loc in hydro])
            base = f"{hydro_str}{junction_str}-{_MULTIPLIER[n]}{ring_word}"
            # P-28.2.3: indicated hydrogen, if any, is cited at the very
            # front of the name -- ahead of the substituent prefix too
            # (same placement `_ring_assembly.py`'s own N=2 case uses).
            name = indicated_hydrogen_prefix + (base if not prefix else f"{prefix}-{base}")

            key = (junction_locant_set, junction_citation, tuple(hydro), sub_locant_set, sub_citation, name)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    return best_name
