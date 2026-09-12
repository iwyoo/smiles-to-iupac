"""Fusion-locant-letter naming for a named five-membered heteromonocycle
(furan, thiophene) ortho-fused onto pyridine, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-25.3.2.4(a): when the two components have different heteroatoms, the
  one with the heteroatom earlier in the order N > F > Cl > Br > I > O >
  S > Se > Te > ... is the base (parent) component. Pyridine's N always
  outranks furan/thiophene's O/S, so pyridine is always the base
  component here, never the attached one - this mirrors
  `_two_component_heterocycle_fusion.py`'s own seniority check, just with
  a fixed winner instead of a per-pair comparison.
- P-25.3.1.3: the base component's peripheral bonds are lettered a, b,
  c... starting at the 1,2-bond, walking around the ring; the fusion bond
  gets the earliest possible letter. Pyridine's own numbering fixes N at
  locant 1 independent of fusion (same convention furan/thiophene already
  use as the *attached* component in
  `_two_component_heterocycle_fusion.py`), leaving exactly one free
  choice - which of the two ring-walk directions from N is "1, 2, 3..." -
  and that choice is made to minimize the fusion bond's letter
  (`_base_bond_letter`). Unlike the plain-benzo case, this letter is not
  fixed to any one value: it depends on which pyridine bond the attached
  ring is fused at (confirmed against real registered isomers below,
  spanning letters 'b' and 'c'). The attached component's own locants are
  then cited in the order matching the base's lettering direction, the
  same rule already implemented for the plain-benzo case in
  `_fusion_locant_letter.py` and for the two-heteromonocycle case in
  `_two_component_heterocycle_fusion.py`.
- The attached ring's own fusion locants must include one atom adjacent
  to its own heteroatom (reusing `_two_component_heterocycle_fusion.py`'s
  `_local_numbering` convention unchanged) - this still allows citations
  like "2,3" (touching the heteroatom's immediate neighbor) as well as
  "3,2", just not a fusion bond entirely on the far side of the ring from
  the heteroatom.
- Fusion bonds touching pyridine's own nitrogen are out of scope (would
  require citing an indicated-nitrogen/bridgehead special case this
  module doesn't attempt).
- Pyrrole is deliberately excluded from this step's scope even though its
  seniority relative to pyridine works out the same way (N ties on
  element, pyridine's larger ring wins via P-25.3.2.4(c)): every real
  pyrrolo-pyridine isomer needs an indicated-hydrogen citation ("1H-") for
  its N-H, a distinct P-14.7 mechanism this module doesn't implement -
  silently omitting it would produce a real IUPAC name with a token
  missing rather than a rejection, so pyrrole is excluded outright
  instead (`UnsupportedStructure` via the heteroatom-not-recognized path
  below) until indicated hydrogen is handled.

Real registered examples confirming this family and scope: furo[3,2-b]
pyridine (PubChem CID 12210217) and furo[2,3-b]pyridine (CID 12421098),
thieno[3,2-b]pyridine (CID 12210218) and thieno[2,3-b]pyridine (CID
289928) - each pair being the two ways the attached ring can straddle a
given pyridine bond - plus furo[2,3-c]pyridine (CID 12826108), confirming
the base letter isn't fixed at 'b'.

Scope, deliberately narrow, matching the sibling fusion modules:
- Exactly two rings: pyridine (base) plus one furan/thiophene (attached),
  ortho-fused, unsubstituted, no extra atoms.
- The fusion bond must touch (or be adjacent to) the attached ring's own
  heteroatom per `_local_numbering`'s rule, and must not touch pyridine's
  own nitrogen.

Anything else (3+ rings, other heteroatoms including pyrrole, fusion at
pyridine's own nitrogen, fusion not reachable from the attached ring's
heteroatom, any other substituent) raises `UnsupportedStructure` and
falls through to other dispatch branches in core.py, exactly like every
other retained/computed-name module in this project.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, ring_cycle

_ATTACHED_RING_NAMES = {
    8: "furo",
    16: "thieno",
}


def _local_numbering(graph, ring_atoms, heteroatom, fusion_atoms):
    """Number a ring 1 (heteroatom) .. n, walking in whichever direction
    makes a fusion atom appear at position 2 (P-25.3.1.3's own
    lowest-locants choice). Returns {atom: locant}, or None if neither
    ring-neighbor of the heteroatom is a fusion atom (the fusion bond
    doesn't touch this ring's own heteroatom)."""
    cycle = ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    start = cycle.index(heteroatom)
    for step in (1, -1):
        neighbor = cycle[(start + step) % n]
        if neighbor in fusion_atoms:
            return {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
    return None


def _base_bond_letter(graph, ring_atoms, nitrogen, fusion_atoms):
    """Letter (a, b, c...) of pyridine's own fusion bond: N is fixed at
    locant 1, and the ring-walk direction (the only free choice once N is
    fixed) is picked to give the fusion bond the earliest possible
    letter. Returns (letter, {atom: locant}), or None if the fusion bond
    touches N itself (out of scope) or isn't a single peripheral bond."""
    cycle = ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    start = cycle.index(nitrogen)
    best = None
    for step in (1, -1):
        numbering = {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
        if numbering[nitrogen] in (numbering[a] for a in fusion_atoms):
            continue
        locants = sorted(numbering[a] for a in fusion_atoms)
        if locants[1] - locants[0] == 1:
            position = locants[0]
        elif locants == [1, n]:
            position = n
        else:
            continue
        letter = chr(ord("a") + position - 1)
        if best is None or letter < best[0]:
            best = (letter, numbering)
    return best


def find_pyridine_heterocycle_fusion_core(mol):
    """Return (pyridine_ring, attached_ring, nitrogen, attached_hetero,
    fusion_atoms) if `mol` is exactly pyridine ortho-fused to one furan/
    thiophene/pyrrole, else None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    sizes = sorted(len(r) for r in atom_rings)
    if sizes != [5, 6]:
        return None
    pyridine_ring = next(r for r in atom_rings if len(r) == 6)
    attached_ring = next(r for r in atom_rings if len(r) == 5)

    for ring in (pyridine_ring, attached_ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                return None
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None

    pyridine_hetero = [idx for idx in pyridine_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(pyridine_hetero) != 1 or mol.GetAtomWithIdx(pyridine_hetero[0]).GetAtomicNum() != 7:
        return None
    nitrogen = pyridine_hetero[0]

    attached_hetero = [idx for idx in attached_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(attached_hetero) != 1:
        return None
    attached_heteroatom = attached_hetero[0]
    if mol.GetAtomWithIdx(attached_heteroatom).GetAtomicNum() not in _ATTACHED_RING_NAMES:
        return None

    shared = set(pyridine_ring) & set(attached_ring)
    if len(shared) != 2:
        return None
    a, b = shared
    if b not in {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()}:
        return None

    if mol.GetNumAtoms() != len(set(pyridine_ring) | set(attached_ring)):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    return pyridine_ring, attached_ring, nitrogen, attached_heteroatom, shared


def has_pyridine_heterocycle_fusion_name(mol) -> bool:
    return find_pyridine_heterocycle_fusion_core(mol) is not None


def name_pyridine_heterocycle_fusion(mol) -> str:
    core = find_pyridine_heterocycle_fusion_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "this two-ring system is not a supported pyridine + named "
            "five-membered heterocycle ortho-fusion (see P-25.3.1.3)"
        )
    pyridine_ring, attached_ring, nitrogen, attached_heteroatom, fusion_atoms = core
    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}

    base = _base_bond_letter(graph, pyridine_ring, nitrogen, fusion_atoms)
    if base is None:
        raise UnsupportedStructure(
            "fusion at pyridine's own nitrogen is not supported (see P-25.3.1.3)"
        )
    letter, base_numbering = base

    attached_numbering = _local_numbering(graph, attached_ring, attached_heteroatom, fusion_atoms)
    if attached_numbering is None:
        raise UnsupportedStructure(
            "the fusion bond does not touch the attached ring's own "
            "heteroatom (only the 'b'-lettered fusion is supported, see "
            "P-25.3.1.3)"
        )

    base_low, base_high = sorted(fusion_atoms, key=lambda atom: base_numbering[atom])
    citation = f"{attached_numbering[base_low]},{attached_numbering[base_high]}"

    attached_atomic_num = mol.GetAtomWithIdx(attached_heteroatom).GetAtomicNum()
    attached_prefix = _ATTACHED_RING_NAMES[attached_atomic_num]
    return f"{attached_prefix}[{citation}-{letter}]pyridine"
