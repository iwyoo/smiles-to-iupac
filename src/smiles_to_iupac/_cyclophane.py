"""Naming of a symmetric phane ring assembly -- N identical, unsubstituted
benzene "superatoms" (N >= 2), each attached to its two ring neighbors by a
plain -CH2-...-CH2- bridge (bridges not necessarily the same length, #1030),
and every ring using the same para (1,4) or meta (1,3) local attachment
pattern -- per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-26 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): phane
  nomenclature names a ring assembly of a macrocycle ("...phane") built from
  component rings ("superatoms") joined by bridges, using composite locants
  and multiplying prefixes on the component-ring names (e.g. 'benzena' for a
  benzene component). Three worked examples confirmed directly against the
  Blue Book text (P-26.3.2.1/P-26.4.1.4, all three explicitly marked "(PIN)"):
  '1,4(1,4)-dibenzenacyclohexaphane' ("[2.2]paracyclophane"),
  '1,4(1,3)-dibenzenacyclohexaphane' ("[2.2]metacyclophane"), and
  '1,3,5,7(1,3)-tetrabenzenacyclooctaphane' ("[1.1.1.1]metacyclophane" /
  cyclotetrabenzylene) -- these three were previously recognized by exact
  whole-molecule hardcoded structure match, one at a time, because no
  general phane algorithm existed. This module instead detects the shared
  *shape* those three examples all happen to
  have (identical rings, symmetric attachment points, equal bridge length)
  and computes the name algorithmically for any N/L/pattern combination
  fitting that shape -- not just these three specific structures.
- Composite locant derivation: each superatom counts as ONE macrocycle
  position (not six; P-26's superatom simplification), so with N superatoms
  and N bridges, the macrocycle has N + sum(bridge lengths) total positions
  and superatom `i`'s locant is 1 plus the running total of (1 + bridge
  length) for every superatom before it in the walk. When every bridge
  happens to share the same length L this reduces to the simpler
  1 + i*(1+L) form, reproducing '1,4' for N=2/L=2 (macrocycle size 6,
  "cyclohexaphane") and '1,3,5,7' for N=4/L=1 (macrocycle size 8,
  "cyclooctaphane") exactly. Unlike the equal-bridge case, the locant *set*
  now genuinely depends on which ring is picked as the walk's start and
  which direction it goes (P-26.4.1.1: lowest set of locants, compared
  term by term in ascending order) -- `name_cyclophane` tries all N
  rotations of the bridge sequence in both directions (2N candidates) and
  keeps whichever gives the lexicographically lowest locant list; putting
  the shortest bridge first, ties broken by the next bridge, and so on,
  always wins by construction, so this reduces to the equal-bridge case's
  single fixed answer whenever every bridge length is the same. The
  parenthesized local locant ('1,4' for para, '1,3' for meta) is the same
  for every ring since this module requires a uniform pattern. The
  multiplying prefix on 'benzena' (P-26.2.3.1/P-26.2.2.1) is the plain
  basic numerical term ('di', 'tetra', ...), confirmed by both worked
  examples above -- not the 'bis'/'tetrakis' compound-substituent form.
- Since every ring is still required identical (plain benzene) with a
  uniform attachment pattern, this is deliberately narrower than the
  general phane algorithm (which would also need to handle different ring
  kinds -- see `_naphthalene_benzene_phane.py` for that separate,
  partially-built-out direction -- and a non-uniform attachment pattern).

Formulas/structures cross-checked (same three compounds as before this
generalization): [2.2]Paracyclophane (C16H16, PubChem CID 74210),
[2.2]Metacyclophane (C16H16, PubChem CID 137543), [1.1.1.1]Metacyclophane
(C28H24, matches PubChem CID 11740710 via InChIKey). PubChem's own
*computed* "IUPACName" property for any of the three is a von Baeyer
bridged-ring name, not a phane name, so it isn't usable for verifying any
name this module returns -- only the Blue Book's own worked examples are.

A single simple substituent (halogen, hydroxyl, amino, or a plain alkyl
branch) on exactly one component ring is also supported (#1038), cited as
an ordinary prefix before the phane name at that ring's own local locant
-- P-26.4.1.2: an amplificant's local numbering follows the numbering
rules of the parent name it's derived from (plain benzene: lowest
locants), independent of the macrocycle's own composite locants for that
ring. Confirmed directly against a Blue Book PIN worked example citing a
substituent this exact way (P-4, `4-substituted...phane` shape):
`3-chloro-2-(naphthalen-2-yl)-1(2)-naphthalena-3,5(1,4),7(1)-tribenzenaheptaphane
(PIN)`. Verified reachable against 3 real, named PubChem structures of
the simplest case this module handles (a single-substituent [2.2]para-
cyclophane): 4-methyl-, 4-bromo-, and 4-hydroxy-[2.2]paracyclophane --
PubChem's own computed name for these is again a von Baeyer name, not
usable for comparison (same caveat as the unsubstituted structures
above), so verification is by structure/connectivity match. The ring's
two bridge attachment points keep priority for the lowest available
local locants (P-29.2's free-valence-first principle, same reasoning
that already fixes the unsubstituted '1,4'/'1,3' convention); among the
numbering choices that preserve that, the substituent gets whichever
remaining locant is lowest.

Explicitly out of scope (the detector below simply returns None for any of
these, so `core.py`'s existing dispatch continues to raise
`UnsupportedStructure`, unchanged):
- Fewer than two rings, rings sharing atoms (fused, not phane), or a ring
  bridged back to itself.
- Any ring other than plain benzene (P-26 defines other "-ena" superatom
  names too, e.g. naphthalena, but none are supported here).
- A non-uniform attachment pattern (some rings para, others meta), an
  ortho (1,2) attachment pattern, or a ring with other than exactly two
  bridge attachment points.
- A bridge that isn't a plain unbranched -CH2-...-CH2- chain, or a
  superatom "supergraph" that isn't a single simple cycle covering every
  ring exactly once (e.g. a branched or multiply-connected assembly).
- More than one substituent total, or substituents on more than one
  component ring -- deferred to a follow-up (#1038's own scope note).
- A substituent whose own branch loops back into the phane's ring system.
"""

from ._multiplicative_text import enclose
from ._common import adjacency, halogen_substituents
from ._numerals import numerical_term
from ._polyspiro import _ring_cyclic_order
from ._substituents import name_branch

_LOCAL_LOCANTS = {"para": "1,4", "meta": "1,3"}


def _pattern_and_bridges(mol, graph, benzene_rings, ring_of, attachment_pairs):
    """(pattern, bridge_lengths) if every ring's two `attachment_pairs`
    entries share the same para/meta local pattern and the bridges between
    them form a single macrocycle spanning every ring in `benzene_rings`,
    else None. `attachment_pairs`: exactly 2 `(ring_atom, exo_atom)` tuples
    per ring, in `benzene_rings` order -- the two genuine bridge
    attachment points, already disambiguated from any substituent by the
    caller."""
    n = len(benzene_rings)
    patterns = []
    for ring, atts in zip(benzene_rings, attachment_pairs):
        cyc = _ring_cyclic_order(graph, list(ring), atts[0][0])
        dist = cyc.index(atts[1][0])
        dist = min(dist, len(cyc) - dist)
        if dist == 3:
            patterns.append("para")
        elif dist == 2:
            patterns.append("meta")
        else:
            return None
    if len(set(patterns)) != 1:
        return None
    pattern = patterns[0]

    # Each bridge is its own edge, identified by its own index below (not
    # by which pair of rings it connects) -- with only 2 rings, both
    # bridges connect the same ring pair, so a ring-pair-keyed lookup
    # can't tell two parallel bridges of different lengths apart.
    edges = []
    seen_starts = set()
    incident = {i: [] for i in range(n)}
    for ring_idx, atts in enumerate(attachment_pairs):
        for ring_atom, start in atts:
            if start in seen_starts:
                continue
            chain = [start]
            prev, cur = ring_atom, start
            while True:
                atom = mol.GetAtomWithIdx(cur)
                if (
                    atom.GetAtomicNum() != 6
                    or atom.GetIsAromatic()
                    or atom.GetFormalCharge() != 0
                    or atom.GetIsotope() != 0
                    or atom.GetTotalNumHs() != 2
                ):
                    return None
                nbrs = [nb for nb in graph[cur] if nb != prev]
                if len(nbrs) != 1:
                    return None
                nxt = nbrs[0]
                if nxt in ring_of:
                    other_ring = ring_of[nxt]
                    if other_ring == ring_idx:
                        return None
                    seen_starts.add(start)
                    seen_starts.add(chain[-1])
                    edge_idx = len(edges)
                    edges.append((ring_idx, other_ring, len(chain)))
                    incident[ring_idx].append(edge_idx)
                    incident[other_ring].append(edge_idx)
                    break
                if len(chain) > 20:
                    return None
                chain.append(nxt)
                prev, cur = cur, nxt

    if any(len(edge_idxs) != 2 for edge_idxs in incident.values()):
        return None
    # Every ring has exactly 2 super-edges; walk them one at a time
    # (always the incident edge not just arrived on) to confirm they form
    # one single cycle spanning all N rings, not several disjoint smaller
    # cycles, collecting each bridge's length in that same walk order.
    order = [0]
    bridge_lengths = []
    used_edge, current = None, 0
    while len(order) < n:
        next_edge = next((e for e in incident[current] if e != used_edge), None)
        if next_edge is None:
            return None
        a, b, length = edges[next_edge]
        next_ring = b if a == current else a
        if next_ring in order:
            return None
        bridge_lengths.append(length)
        order.append(next_ring)
        used_edge, current = next_edge, next_ring
    closing_edge = next((e for e in incident[current] if e != used_edge), None)
    if closing_edge is None:
        return None
    a, b, length = edges[closing_edge]
    closing_ring = b if a == current else a
    if closing_ring != 0:
        return None
    bridge_lengths.append(length)
    return pattern, bridge_lengths


def _branch_touches_rings(graph, root, coming_from, ring_of, limit=30):
    """True if the branch hanging off `root` (reached from `coming_from`)
    ever reaches one of the phane's own component-ring atoms -- ruling out
    a mis-picked "substituent" that's really another bridge, or a
    substituent fused back into the ring system, both out of scope."""
    seen = {coming_from}
    stack = [root]
    while stack:
        if len(seen) > limit:
            return True
        atom_idx = stack.pop()
        if atom_idx in seen:
            continue
        seen.add(atom_idx)
        if atom_idx in ring_of:
            return True
        stack.extend(n for n in graph[atom_idx] if n not in seen)
    return False


def _find_symmetric_phane(mol):
    """Return `(ring_count, bridge_lengths, pattern, substituent)` if `mol`
    fits the shape this module supports (see module docstring), else None.
    `bridge_lengths` is a list of `ring_count` bridge lengths, in the same
    cyclic ring order as the walk that discovers them (arbitrary starting
    ring/direction -- `name_cyclophane` tries every rotation/direction
    itself). `substituent` is None (the plain, fully symmetric case) or
    `(ring_atom, exo_atom, other_att_a, other_att_b)` for the one
    component ring carrying a substituent -- `ring_atom`/`exo_atom` the
    substituted ring atom and its exocyclic branch root, `other_att_a`/
    `other_att_b` that same ring's two genuine bridge attachment atoms
    (needed to re-derive the ring's own local numbering when naming it)."""
    ring_info = mol.GetRingInfo()
    benzene_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        if not all(
            mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring
        ):
            continue
        benzene_rings.append(set(ring))
    n = len(benzene_rings)
    if n < 2:
        return None
    for i in range(n):
        for j in range(i + 1, n):
            if benzene_rings[i] & benzene_rings[j]:
                return None

    graph = adjacency(mol)
    ring_of = {a: idx for idx, ring in enumerate(benzene_rings) for a in ring}

    attachments = []
    for ring in benzene_rings:
        atts = []
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None
            exo = [nb for nb in graph[a] if nb not in ring]
            if not exo:
                if atom.GetDegree() != 2 or atom.GetTotalNumHs() != 1:
                    return None
                continue
            if len(exo) != 1 or atom.GetDegree() != 3 or atom.GetTotalNumHs() != 0:
                return None
            atts.append((a, exo[0]))
        if len(atts) not in (2, 3):
            return None
        attachments.append(atts)

    substituted_rings = [i for i, atts in enumerate(attachments) if len(atts) == 3]
    if len(substituted_rings) > 1:
        return None

    if not substituted_rings:
        result = _pattern_and_bridges(mol, graph, benzene_rings, ring_of, attachments)
        if result is None:
            return None
        pattern, bridge_lengths = result
        return n, bridge_lengths, pattern, None

    (sub_ring_idx,) = substituted_rings
    candidates3 = attachments[sub_ring_idx]
    successes = []
    for drop in range(3):
        bridge_pair = [candidates3[i] for i in range(3) if i != drop]
        substituent_ring_atom, substituent_exo_atom = candidates3[drop]
        trial = list(attachments)
        trial[sub_ring_idx] = bridge_pair
        result = _pattern_and_bridges(mol, graph, benzene_rings, ring_of, trial)
        if result is None:
            continue
        if _branch_touches_rings(graph, substituent_exo_atom, substituent_ring_atom, ring_of):
            continue
        pattern, bridge_lengths = result
        successes.append((pattern, bridge_lengths, substituent_ring_atom, substituent_exo_atom, bridge_pair))

    if len(successes) != 1:
        return None
    pattern, bridge_lengths, substituent_ring_atom, substituent_exo_atom, bridge_pair = successes[0]
    (att_a, _), (att_b, _) = bridge_pair
    substituent = (substituent_ring_atom, substituent_exo_atom, att_a, att_b)
    return n, bridge_lengths, pattern, substituent


def has_cyclophane_name(mol) -> bool:
    return _find_symmetric_phane(mol) is not None


def _macro_locants(bridge_lengths):
    locants = []
    position = 1
    for length in bridge_lengths:
        locants.append(position)
        position += 1 + length
    return locants


def _substituent_ring_prefix(mol, graph, substituent):
    """`"<local-locant>-<name>"` for the one substituent `_find_symmetric_
    phane` found, per P-26.4.1.2: the substituted ring's own local
    numbering follows plain benzene numbering with the two bridge
    attachment points (the ring's "free valences" into the macrocycle)
    kept at the lowest available locants first (P-29.2), same principle
    that already fixes the unsubstituted case's '1,4'/'1,3' convention.
    For a para ring this alone leaves a genuine tie -- the other
    attachment sits exactly 3 steps away in *either* direction around a
    6-ring, so both directions from a given start atom already put it at
    local position 4 -- so all 4 combinations of (start atom, direction)
    are tried and whichever gives the lowest substituent locant among
    those tied for the lowest attachment locant wins (confirmed against
    two independently-built SMILES for the same real structure giving the
    same locant only once this tie was handled explicitly, not just
    picking `_ring_cyclic_order`'s own arbitrary walk direction)."""
    ring_atom, exo_atom, att_a, att_b = substituent
    ring = next(
        r
        for r in mol.GetRingInfo().AtomRings()
        if ring_atom in r and len(r) == 6 and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)
    )
    ring_set = set(ring)
    candidates = []
    for start, other in ((att_a, att_b), (att_b, att_a)):
        cyc = _ring_cyclic_order(graph, list(ring_set), start)
        idx = cyc.index(other)
        candidates.append((idx, cyc))
        candidates.append((len(cyc) - idx, [cyc[0]] + list(reversed(cyc[1:]))))
    best_other_locant = min(idx for idx, _ in candidates)
    tied = [cyc for idx, cyc in candidates if idx == best_other_locant]
    best_order = min(tied, key=lambda order: order.index(ring_atom))
    local_locant = best_order.index(ring_atom) + 1

    prefixes = dict(halogen_substituents(mol))
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1:
            prefixes[atom.GetIdx()] = "hydroxy"
        elif atom.GetAtomicNum() == 7 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2:
            prefixes[atom.GetIdx()] = "amino"

    name, is_compound = name_branch(graph, exo_atom, ring_atom, prefixes, mol=mol)
    display = enclose(name) if is_compound else name
    return f"{local_locant}-{display}"


def name_cyclophane(mol) -> str:
    n, bridge_lengths, pattern, substituent = _find_symmetric_phane(mol)

    # P-26.4.1.1: the lowest set of superatom locants wins, compared term
    # by term in ascending order -- try every rotation of the bridge
    # sequence in both directions (2N candidates; harmless duplicates when
    # the sequence has its own symmetry) and keep the smallest resulting
    # locant list. Putting the shortest bridge first always wins the first
    # point of difference, so this reduces to the single fixed answer the
    # equal-bridge case already had when every bridge length is the same.
    reversed_lengths = list(reversed(bridge_lengths))
    candidates = [
        seq[start:] + seq[:start] for seq in (bridge_lengths, reversed_lengths) for start in range(n)
    ]
    best_locants = min(_macro_locants(candidate) for candidate in candidates)

    macro_str = ",".join(str(loc) for loc in best_locants)
    local_str = _LOCAL_LOCANTS[pattern]
    macrocycle_size = n + sum(bridge_lengths)
    multiplier = numerical_term(n)
    phane_name = f"{macro_str}({local_str})-{multiplier}benzenacyclo{numerical_term(macrocycle_size)}phane"
    if substituent is None:
        return phane_name
    graph = adjacency(mol)
    return f"{_substituent_ring_prefix(mol, graph, substituent)}-{phane_name}"
