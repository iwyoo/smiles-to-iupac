"""Fusion-locant-letter naming for two five-membered heteromonocycles
(thiophene, furan -- identical or mixed) ortho-fused to each other, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-25.3.1.3 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  each peripheral bond of the parent (base) component is lettered a, b, c...
  starting at the 1,2-bond; the earliest possible letter is used for the
  fusion bond, and the attached component's own locants at the two fusion
  atoms are cited in the order that matches the parent's lettering
  direction. Worked example given verbatim in this section for the exact
  same shape (two identical five-membered chalcogen heterocycles
  self-fused): 'selenopheno[2,3-b]selenophene (PIN)',
  'selenopheno[3,4-b]selenophene (PIN)', 'selenopheno[3,2-b]selenophene
  (PIN)' -- confirming both that fusion letters are computed (not a
  numeral-prefix retained-name shortcut, unlike benzofuran/benzothiophene
  in `_heteroaromatic_fused.py`) and that a third isomer (letter 'c', the
  fusion bond not touching either ring's own heteroatom) exists and must be
  told apart from the two 'b'-lettered ones handled here.
- P-25.3.2.4, criterion (a): when the two components are *different* heteromonocycles,
  the parent (base) component is the one containing the heteroatom earlier
  in the seniority order N > F > Cl > Br > I > O > S > Se > Te > ...; the
  Blue Book's own worked example for this exact O-vs-S case ('2H-[1,4]dithiepino
  [2,3-c]furan (PIN)') is annotated directly in the text as "furan is
  senior to dithiepine; O > S" -- confirmed against the primary source
  before implementing, not assumed. That ordering (O > S > Se > Te) covers
  all four chalcogens this module recognizes: whichever ring's heteroatom
  comes first in O/S/Se/Te is always the base component when the two
  rings differ. When they're the same, the ordering has no effect (either
  ring could be "base" and the result is identical since both keep the
  same ring name), so ring index 0 is simply used as a tiebreak.
- P-25.2.1 / Table 2.2: thiophene and furan both fix their own heteroatom
  at locant 1 (independent of fusion), with a mirror symmetry across the
  heteroatom that makes both ring-numbering directions equally valid until
  a fusion (or substituent) breaks the tie via "lowest locants".
- P-25.3.3.1.1 (whole-fused-system peripheral numbering): unlike the
  6-5 pyridine+heterocycle shape in `_pyridine_heterocycle_fusion.py`
  (where the indicated-hydrogen locant is always 1, independent of
  orientation, because the attached 5-ring is always numbered first), a
  plain 5-5 bicyclic has no such fixed shortcut -- the two rings are the
  same size, so which one is numbered first, and which of its own two
  ways of being drawn (P-25.3.2.3.1's pentagon shape has two mirror
  orientations) is preferred, is a genuine choice resolved by
  `_whole_system_locants` below: enumerate the four geometrically valid
  embeddings for a plain 2-ring ortho-fusion (which physical ring sits on
  the right x which of the two shared fusion atoms is drawn "on top" --
  both freedoms P-25.3.2.3.3 leaves open for a symmetric 2-ring row),
  number each one per P-25.3.3.1.1 (start at the
  nonfused atom immediately clockwise from the top fusion atom, proceed
  clockwise), then pick the winner via P-25.3.3.1.2's tie-break cascade
  (a) low locants to heteroatoms as a set, (b) then by element seniority
  F > Cl > Br > I > O > S > Se > Te > N > P > ... (a *different* order
  from P-25.3.2.4(a)'s base-component seniority below -- confirmed
  against the primary source, not assumed), (c) low locants to fusion
  atoms, (f) low locants to indicated hydrogen. Verified against three
  real PubChem structures spanning all three outcomes this shape can
  have: furo[3,2-b]pyrrole (CID 21679450) -> '4H-', furo[2,3-b]pyrrole
  (CID 21868890) -> '6H-', pyrrolo[3,2-b]pyrrole (CID 20203864) -> no
  indicated hydrogen at all (both ring nitrogens end up pyridine-type).

Scope, deliberately narrow:
- Exactly two rings, both aromatic, both 5-membered, both with exactly one
  ring heteroatom, each one of N/O/S/Se/Te (identical or any mixed pair).
- Ortho-fusion only (the two rings share exactly one bond).
- The fusion bond must touch an atom adjacent to *each* ring's own
  heteroatom (letter 'b', per the worked examples above) -- a fusion bond
  positioned elsewhere on the ring (e.g. letter 'c', the
  'selenopheno[3,4-b]selenophene'/'thieno[2,3-c]thiophene' shape,
  independently confirmed against PubChem CID 520171) is out of scope and
  explicitly rejected, not silently misnamed.
- Completely unsubstituted (no atoms beyond the 8 ring atoms). Indicated
  hydrogen *is* handled (pyrrole's N-H), computed via the whole-system
  numbering above; furan/thiophene/selenophene/tellurophene's own
  heteroatoms never carry one.

Anything else (3+ components, other heteroatoms, non-'b' fusion bonds, any
substituent) raises `UnsupportedStructure` and falls through to other
dispatch branches in core.py, exactly like every other retained/computed-
name module in this project.
"""

from ._common import UnsupportedStructure, ring_cycle

_RING_NAMES = {
    7: ("pyrrolo", "pyrrole"),
    8: ("furo", "furan"),
    16: ("thieno", "thiophene"),
    34: ("selenopheno", "selenophene"),
    52: ("telluropheno", "tellurophene"),
}
# P-25.3.2.4(a)'s base-*component*-selection seniority order (N > F > Cl >
# Br > I > O > S > Se > Te > ...), restricted to the elements this module's
# five rings use -- lower value is more senior (becomes the base component).
_SENIORITY_ORDER = {7: 0, 8: 1, 16: 2, 34: 3, 52: 4}
# P-25.3.3.1.2(b)'s *separate* numbering-locant tie-break seniority order
# (low locants to heteroatoms in the order F > Cl > Br > I > O > S > Se >
# Te > N > P > ...) -- deliberately different from `_SENIORITY_ORDER`
# above: N is the most senior element for choosing the base *component*,
# but the *least* senior of this module's elements when only choosing which
# of two tied numbering candidates gets the lower *locant*. Confirmed
# against the primary source text directly (both rule numbers give an
# explicit, different heteroatom order), not assumed from a naming
# resemblance to P-25.3.2.4(a). Restricted the same way as `_SENIORITY_ORDER`.
_NUMBERING_SENIORITY_ORDER = {8: 0, 16: 1, 34: 2, 52: 3, 7: 4}


def _local_numbering(graph, ring_atoms, heteroatom, fusion_atoms):
    """Number a ring 1 (heteroatom) .. 5, walking in whichever direction
    makes a fusion atom appear at position 2 (the 'lowest locants' choice
    forced by P-25.3.1.3, see module docstring). Returns {atom: locant}, or
    None if neither ring-neighbor of the heteroatom is a fusion atom (the
    fusion bond doesn't touch this ring's own heteroatom -- letter 'c',
    out of scope)."""
    cycle = ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    start = cycle.index(heteroatom)
    for step in (1, -1):
        neighbor = cycle[(start + step) % n]
        if neighbor in fusion_atoms:
            return {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
    return None


def _ring_arc(graph, ring_atoms, start, other_fusion_atom):
    """Walk `ring_atoms`'s cycle from `start`, in whichever direction leads
    away from `other_fusion_atom` first, collecting every atom up to but not
    including `other_fusion_atom`. This is one ring's contribution to a
    whole-system peripheral walk: `start` is a fusion atom being treated as
    already-numbered, and the returned list is that ring's own non-fusion
    atoms in P-25.3.3.1.1 walk order."""
    cycle = ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    start_idx = cycle.index(start)
    for step in (1, -1):
        if cycle[(start_idx + step) % n] == other_fusion_atom:
            continue
        arc = []
        i = start_idx
        while True:
            i = (i + step) % n
            atom = cycle[i]
            if atom == other_fusion_atom:
                return arc
            arc.append(atom)
    raise AssertionError("a ring atom always has exactly one non-fusion-partner neighbor at a fusion atom")


def _embedding_locants(graph, right_ring, left_ring, top, bottom):
    """P-25.3.3.1.1's numbering for one specific geometrically valid
    embedding: `right_ring` drawn on the right with `top` as the upper of
    its two shared fusion atoms (`bottom` the lower). Numbering starts at
    the nonfused atom immediately clockwise from `top` (i.e. `top`'s
    neighbor in `right_ring` other than `bottom`), proceeds clockwise
    through the rest of `right_ring`, across `bottom` (a fusion atom), on
    through `left_ring`, and back to `top` (the other fusion atom) last.
    Returns {atom: (int_locant, is_fusion)}."""
    right_arc = _ring_arc(graph, right_ring, top, bottom)
    left_arc = _ring_arc(graph, left_ring, bottom, top)
    sequence = right_arc + [bottom] + left_arc + [top]
    locants = {}
    counter = 0
    for atom in sequence:
        if atom in (top, bottom):
            locants[atom] = (counter, True)
        else:
            counter += 1
            locants[atom] = (counter, False)
    return locants


def _whole_system_locants(graph, ring_atom_sets, fusion_atoms, mol):
    """The winning P-25.3.3.1.1/.1.2 whole-system numbering for a plain
    2-ring ortho-fused bicyclic: among the four geometrically valid
    embeddings (which ring is on the right x which fusion atom is on top,
    see module docstring), pick the one P-25.3.3.1.2's tie-break cascade
    prefers. Returns {atom: (int_locant, is_fusion)} for the winner."""
    p, q = tuple(fusion_atoms)
    ring_a, ring_b = ring_atom_sets
    candidates = [
        _embedding_locants(graph, right, left, top, bottom)
        for right, left in ((ring_a, ring_b), (ring_b, ring_a))
        for top, bottom in ((p, q), (q, p))
    ]

    def heteroatoms_by_locant(locants):
        return sorted(
            (locant, atom) for atom, (locant, _fusion) in locants.items() if mol.GetAtomWithIdx(atom).GetAtomicNum() != 6
        )

    # (a) low locants to heteroatoms as a set.
    candidates = _min_by(candidates, lambda loc: [n for n, _a in heteroatoms_by_locant(loc)])
    # (b) among ties, low locants in P-25.3.3.1.2(b)'s element order.
    candidates = _min_by(
        candidates,
        lambda loc: [
            n
            for n, a in sorted(
                heteroatoms_by_locant(loc), key=lambda pair: _NUMBERING_SENIORITY_ORDER[mol.GetAtomWithIdx(pair[1]).GetAtomicNum()]
            )
        ],
    )
    # (c) low locants to fusion atoms.
    candidates = _min_by(candidates, lambda loc: sorted(n for _a, (n, fusion) in loc.items() if fusion))
    # (f) low locants to indicated hydrogen atoms -- restricted to
    # heteroatoms: a ring CH carbon's implicit hydrogen is part of the
    # ordinary mancude pattern, not an "indicated hydrogen" (P-14.7's own
    # sense, e.g. pyrrole's N-H vs pyridine-type N with none).
    candidates = _min_by(
        candidates,
        lambda loc: sorted(
            n for a, (n, _fusion) in loc.items() if mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and mol.GetAtomWithIdx(a).GetTotalNumHs() > 0
        ),
    )
    return candidates[0]


def _min_by(candidates, key):
    keyed = [(key(c), c) for c in candidates]
    best = min(k for k, _c in keyed)
    return [c for k, c in keyed if k == best]


def _indicated_hydrogen_prefix(locants, mol):
    """"<n>H-" for whichever atom carries an explicit hydrogen in the
    winning whole-system numbering (pyrrole's N-H), else ""."""
    for atom, (locant, is_fusion) in locants.items():
        if is_fusion:
            continue
        heavy_atom = mol.GetAtomWithIdx(atom)
        if heavy_atom.GetAtomicNum() != 6 and heavy_atom.GetTotalNumHs() > 0:
            return f"{locant}H-"
    return ""


def find_two_component_heterocycle_fusion_core(mol):
    """Return (ring_atom_sets, fusion_atoms, heteroatoms, elements) if `mol`
    is exactly two ortho-fused 5-membered aromatic rings each with one
    O/S/Se/Te heteroatom (identical or mixed) and no other atoms, else
    None."""
    if mol.GetNumAtoms() != 8:
        return None
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2 or any(len(r) != 5 for r in atom_rings):
        return None
    ring_atom_sets = [set(r) for r in atom_rings]
    for ring in ring_atom_sets:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                return None

    heteroatoms = []
    for ring in ring_atom_sets:
        hetero = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        if len(hetero) != 1 or mol.GetAtomWithIdx(hetero[0]).GetAtomicNum() not in _RING_NAMES:
            return None
        heteroatoms.append(hetero[0])
    elements = tuple(mol.GetAtomWithIdx(h).GetAtomicNum() for h in heteroatoms)

    shared = ring_atom_sets[0] & ring_atom_sets[1]
    if len(shared) != 2:
        return None
    a, b = shared
    if b not in {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()}:
        return None
    if heteroatoms[0] in shared or heteroatoms[1] in shared:
        return None

    return ring_atom_sets, shared, heteroatoms, elements


def has_two_component_heterocycle_fusion_name(mol) -> bool:
    return find_two_component_heterocycle_fusion_core(mol) is not None


def name_two_component_heterocycle_fusion(mol) -> str:
    core = find_two_component_heterocycle_fusion_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "this two-ring heteroaromatic system is not a supported "
            "two-component ortho-fusion (see P-25.3.1.3)"
        )
    ring_atom_sets, fusion_atoms, heteroatoms, elements = core
    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}

    # P-25.3.2.4(a): the more senior heteroatom's ring is the base component;
    # a tie (identical rings) is broken by ring index, which doesn't affect
    # the resulting name string (see module docstring).
    if _SENIORITY_ORDER[elements[0]] <= _SENIORITY_ORDER[elements[1]]:
        base_idx, attached_idx = 0, 1
    else:
        base_idx, attached_idx = 1, 0

    base_numbering = _local_numbering(graph, ring_atom_sets[base_idx], heteroatoms[base_idx], fusion_atoms)
    attached_numbering = _local_numbering(graph, ring_atom_sets[attached_idx], heteroatoms[attached_idx], fusion_atoms)
    if base_numbering is None or attached_numbering is None:
        raise UnsupportedStructure(
            "the fusion bond does not touch either ring's own heteroatom "
            "(only the 'b'-lettered fusion is supported, see P-25.3.1.3)"
        )

    base_low, base_high = sorted(fusion_atoms, key=lambda atom: base_numbering[atom])
    citation = f"{attached_numbering[base_low]},{attached_numbering[base_high]}"

    attached_prefix, _ = _RING_NAMES[elements[attached_idx]]
    _, base_name = _RING_NAMES[elements[base_idx]]

    whole_system_locants = _whole_system_locants(graph, ring_atom_sets, fusion_atoms, mol)
    indicated_h = _indicated_hydrogen_prefix(whole_system_locants, mol)

    return f"{indicated_h}{attached_prefix}[{citation}-b]{base_name}"
