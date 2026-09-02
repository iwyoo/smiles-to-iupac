"""Naming of a symmetric phane ring assembly -- N identical, unsubstituted
benzene "superatoms" (N >= 2), each attached to its two ring neighbors by a
bridge of the same length L, all bridges the same length, and every ring
using the same para (1,4) or meta (1,3) local attachment pattern -- per the
IUPAC 2013 Recommendations ("the Blue Book"):

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
- Composite locant derivation: since every ring is identical and every
  bridge the same length L, the N superatoms sit at evenly-spaced positions
  around the macrocycle -- each superatom counts as ONE macrocycle position
  (not six; P-26's superatom simplification), so with N superatoms and N
  bridges of L atoms each, the macrocycle has N*(1+L) total positions and
  the superatoms sit at 1, 1+(1+L), 1+2*(1+L), .... This reproduces
  '1,4' for N=2/L=2 (macrocycle size 6, "cyclohexaphane") and '1,3,5,7' for
  N=4/L=1 (macrocycle size 8, "cyclooctaphane") exactly, without needing to
  actually trace a specific starting ring/direction in the input molecule
  (rotational symmetry makes the locant *set* independent of which ring is
  "first"). The parenthesized local locant ('1,4' for para, '1,3' for meta)
  is the same for every ring since this module requires a uniform pattern.
  The multiplying prefix on 'benzena' (P-26.2.3.1/P-26.2.2.1) is the plain
  basic numerical term ('di', 'tetra', ...), confirmed by both worked
  examples above -- not the 'bis'/'tetrakis' compound-substituent form.
- Since every ring/bridge is required identical, this is deliberately a
  narrower shape than the general phane algorithm (which would also need
  to handle different ring kinds, asymmetric attachment points, and a
  non-trivial attachment-locant-ordering rule for those, still open for
  that broader case).

Formulas/structures cross-checked (same three compounds as before this
generalization): [2.2]Paracyclophane (C16H16, PubChem CID 74210),
[2.2]Metacyclophane (C16H16, PubChem CID 137543), [1.1.1.1]Metacyclophane
(C28H24, matches PubChem CID 11740710 via InChIKey). PubChem's own
*computed* "IUPACName" property for any of the three is a von Baeyer
bridged-ring name, not a phane name, so it isn't usable for verifying any
name this module returns -- only the Blue Book's own worked examples are.

Explicitly out of scope (the detector below simply returns None for any of
these, so `core.py`'s existing dispatch continues to raise
`UnsupportedStructure`, unchanged):
- Fewer than two rings, rings sharing atoms (fused, not phane), or a ring
  bridged back to itself.
- Any ring other than plain benzene (P-26 defines other "-ena" superatom
  names too, e.g. naphthalena, but none are supported here).
- A non-uniform attachment pattern (some rings para, others meta), an
  ortho (1,2) attachment pattern, or a ring with other than exactly two
  attachment points.
- Bridges of unequal length, a bridge that isn't a plain unbranched
  -CH2-...-CH2- chain, or a superatom "supergraph" that isn't a single
  simple cycle covering every ring exactly once (e.g. a branched or
  multiply-connected assembly).
- Any substituent anywhere.
"""

from ._common import adjacency
from ._numerals import numerical_term
from ._polyspiro import _ring_cyclic_order

_LOCAL_LOCANTS = {"para": "1,4", "meta": "1,3"}


def _find_symmetric_phane(mol):
    """Return (ring_count, bridge_length, pattern) if `mol` fits the shape
    this module supports (see module docstring), else None."""
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
        if len(atts) != 2:
            return None
        attachments.append(atts)

    patterns = []
    for ring, atts in zip(benzene_rings, attachments):
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

    bridge_lengths = set()
    seen_starts = set()
    super_adj = {i: [] for i in range(n)}
    for ring_idx, atts in enumerate(attachments):
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
                    bridge_lengths.add(len(chain))
                    super_adj[ring_idx].append(other_ring)
                    super_adj[other_ring].append(ring_idx)
                    break
                if len(chain) > 20:
                    return None
                chain.append(nxt)
                prev, cur = cur, nxt

    if len(bridge_lengths) != 1 or any(len(neighbors) != 2 for neighbors in super_adj.values()):
        return None
    # Every ring has exactly 2 super-edges; confirm they form one single
    # cycle spanning all N rings, not several disjoint smaller cycles.
    order = [0]
    previous, current = None, 0
    while len(order) < n:
        next_ring = next((r for r in super_adj[current] if r != previous), None)
        if next_ring is None or next_ring == 0:
            return None
        order.append(next_ring)
        previous, current = current, next_ring
    if 0 not in super_adj[current]:
        return None
    return n, next(iter(bridge_lengths)), pattern


def has_cyclophane_name(mol) -> bool:
    return _find_symmetric_phane(mol) is not None


def name_cyclophane(mol) -> str:
    n, bridge_length, pattern = _find_symmetric_phane(mol)
    macro_locants = [1 + i * (1 + bridge_length) for i in range(n)]
    macro_str = ",".join(str(loc) for loc in macro_locants)
    local_str = _LOCAL_LOCANTS[pattern]
    macrocycle_size = n * (1 + bridge_length)
    multiplier = numerical_term(n)
    return f"{macro_str}({local_str})-{multiplier}benzenacyclo{numerical_term(macrocycle_size)}phane"
