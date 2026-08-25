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
- P-25.3.2.4, criterion (a) (as of `tasks/hetero-two-component-fusion-naming.md`,
  2026-08-25): when the two components are *different* heteromonocycles,
  the parent (base) component is the one containing the heteroatom earlier
  in the seniority order N > F > Cl > Br > I > O > S > Se > Te > ...; the
  Blue Book's own worked example for this exact O-vs-S case ('2H-[1,4]dithiepino
  [2,3-c]furan (PIN)') is annotated directly in the text as "furan is
  senior to dithiepine; O > S" -- confirmed against the primary source
  before implementing, not assumed. Within this module's O/S-only scope,
  that reduces to: a furan component is always the base, a thiophene
  component is always the attached ('thieno'-prefixed) one, whenever the
  two rings' heteroatoms differ. When they're the same (S+S or O+O), the
  ordering has no effect (either ring could be "base" and the result is
  identical since both keep the same ring name), so ring index 0 is simply
  used as a tiebreak.
- P-25.2.1 / Table 2.2: thiophene and furan both fix their own heteroatom
  at locant 1 (independent of fusion), with a mirror symmetry across the
  heteroatom that makes both ring-numbering directions equally valid until
  a fusion (or substituent) breaks the tie via "lowest locants".

Scope, deliberately narrow (see tasks/homo-heterocycle-fusion-naming.md and
tasks/hetero-two-component-fusion-naming.md):
- Exactly two rings, both aromatic, both 5-membered, both with exactly one
  ring heteroatom, each either O or S (identical or mixed) -- any other
  heteroatom (Se, Te, N, ...) is out of scope; P-25.3.2.4(a)'s full
  seniority order is not implemented, only the O-vs-S slice of it that's
  independently verified above.
- Ortho-fusion only (the two rings share exactly one bond).
- The fusion bond must touch an atom adjacent to *each* ring's own
  heteroatom (letter 'b', per the worked examples above) -- a fusion bond
  positioned elsewhere on the ring (e.g. letter 'c', the
  'selenopheno[3,4-b]selenophene'/'thieno[2,3-c]thiophene' shape,
  independently confirmed against PubChem CID 520171) is out of scope and
  explicitly rejected, not silently misnamed.
- Completely unsubstituted (no atoms beyond the 8 ring atoms, no indicated
  hydrogen -- both rings' own heteroatoms are pyridine/thiophene/furan-type
  divalent atoms with no N-H case to consider, since only O/S are in
  scope).

Anything else (3+ components, heteroatoms other than O/S, non-'b' fusion
bonds, any substituent) raises `UnsupportedStructure` and falls through to
other dispatch branches in core.py, exactly like every other
retained/computed-name module in this project.
"""

from ._common import UnsupportedStructure

_RING_NAMES = {8: ("furo", "furan"), 16: ("thieno", "thiophene")}
# P-25.3.2.4(a): heteroatom seniority order, restricted to this module's
# O/S-only scope -- lower value is more senior (becomes the base component).
_SENIORITY_ORDER = {8: 0, 16: 1}


def _ring_cycle(graph, ring_atoms):
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def _local_numbering(graph, ring_atoms, heteroatom, fusion_atoms):
    """Number a ring 1 (heteroatom) .. 5, walking in whichever direction
    makes a fusion atom appear at position 2 (the 'lowest locants' choice
    forced by P-25.3.1.3, see module docstring). Returns {atom: locant}, or
    None if neither ring-neighbor of the heteroatom is a fusion atom (the
    fusion bond doesn't touch this ring's own heteroatom -- letter 'c',
    out of scope)."""
    cycle = _ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    start = cycle.index(heteroatom)
    for step in (1, -1):
        neighbor = cycle[(start + step) % n]
        if neighbor in fusion_atoms:
            return {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
    return None


def find_two_component_heterocycle_fusion_core(mol):
    """Return (ring_atom_sets, fusion_atoms, heteroatoms, elements) if `mol`
    is exactly two ortho-fused 5-membered aromatic rings each with one O/S
    heteroatom (identical or mixed) and no other atoms, else None."""
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
    return f"{attached_prefix}[{citation}-b]{base_name}"
