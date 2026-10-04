"""Fusion-locant-letter naming for a two-ring, ortho-fused aromatic system
whose shared fusion atom is itself a heteroatom common to both components
(P-25.3.2.5.1: "A heteroatom common to two components must be indicated in
the name of each component"), per the IUPAC 2013 Recommendations ("the
Blue Book").

Every other 2-component fusion module in this project explicitly rejects
this shape: `_two_component_heterocycle_fusion.py` has `if heteroatoms[0]
in shared or heteroatoms[1] in shared: return None`, and
`_pyridine_heterocycle_fusion.py` rejects any fusion touching pyridine's
own nitrogen the same way. The Blue Book's own P-25.3.2.5.1 worked example
(imidazo[2,1-b][1,3]thiazole, PubChem CID 817024) is exactly this shape -
a bridgehead nitrogen shared by imidazole and thiazole - so this is a
real, previously out-of-scope family, not a one-off.

Two mechanisms, reused for both the base and attached component:

- `_numbering_for_pattern`: each supported ring type has a fixed cyclic
  heteroatom pattern (e.g. imidazole = N,C,N,C,C starting at locant 1;
  pyridine = N,C,C,C,C,C). For a ring matching that pattern, this tries
  every (start, direction) combination consistent with the ring's own
  fixed numbering convention and returns whichever gives the *lowest*
  locant set for the fusion atoms - the same "lowest locants" principle
  `_fusion_locant_letter.py` already applies via automorphism search, and
  `_local_numbering` applies via direction choice, generalized here to a
  ring with more than one heteroatom whose relative positions are fixed
  by the ring's own retained name (so the free choice is which physical
  atom plays which numbered role, not the pattern itself).
- Seniority (`_ring_senior`, P-25.3.2.4(a)-(e)): compares the two rings'
  (size, heteroatom-multiset) to decide which is the base component. This
  needs a real multi-heteroatom evaluator, not the existing sibling
  module's single-heteroatom `_SENIORITY_ORDER` shortcut - imidazole (N,
  N) and thiazole (N, S) tie on (a) (both contain N), (b) (1 ring each),
  (c) (both same ring size when compared - see below), and (d) (2
  heteroatoms each); only (e), "the greater variety of heteroatoms,"
  decides it (thiazole's N+S beats imidazole's N+N) - confirmed against
  the real name (thiazole is the base). For a 5+6 pair like imidazo[1,2-a]
  pyridine, (c) ("the larger ring at the first point of difference")
  already decides it before (d)/(e) are reached (pyridine's 6-ring beats
  imidazole's 5-ring) - matching that real name's choice of pyridine as
  base too.

Scope, deliberately narrow:
- Exactly two rings, both aromatic, unsubstituted, ortho-fused (sharing
  exactly one bond), with exactly one of the two shared atoms being a
  heteroatom (the bridgehead) - the other shared atom is plain carbon.
- Each ring must match one of the fixed patterns in `_RING_TYPES` below
  (imidazole/pyrazole as the always-5-membered attached component;
  thiazole/pyridine/pyrimidine as the base, 5- or 6-membered) - any other
  ring type, or any ring with more than 2 heteroatoms (e.g. a
  1,3,4-thiadiazole base, which also matches the Blue Book's own
  P-25.3.2.5.1 example family but needs a 3-heteroatom pattern this step
  doesn't add) is out of scope and correctly rejected, not guessed at.

Anything else (3+ rings, unsupported ring types, more than one bridgehead
heteroatom, any substituent) raises `UnsupportedStructure` and falls
through to other dispatch branches in core.py, exactly like every other
retained/computed-name module in this project.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, ring_cycle

# (pattern of atomic numbers for locants 1..n, base name, attached-form prefix)
_RING_TYPES = (
    ((7, 6, 7, 6, 6), "imidazole", "imidazo"),
    ((7, 7, 6, 6, 6), "pyrazole", "pyrazolo"),
    ((16, 6, 7, 6, 6), "thiazole", "thiazolo"),
    ((7, 6, 6, 6, 6, 6), "pyridine", "pyrido"),
    ((7, 6, 7, 6, 6, 6), "pyrimidine", "pyrimido"),
    ((7, 7, 6, 6, 6, 6), "pyridazine", "pyridazino"),
    ((7, 6, 6, 7, 6, 6), "pyrazine", "pyrazino"),
    ((7, 7, 6, 7, 6, 6), "[1,2,4]triazine", "[1,2,4]triazino"),
)

# P-25.3.2.4(a): heteroatom seniority order (abbreviated to the elements
# `_RING_TYPES` above actually uses).
_HETEROATOM_SENIORITY = {7: 0, 8: 1, 16: 2}


def _bond_letter_position(numbering, fusion_atoms, n):
    """Peripheral-bond position (1-indexed, matching a=1, b=2, ...) of the
    bond between `fusion_atoms` under `numbering` - the *wraparound* bond
    (locants {1, n}) is position n (the *last* letter), not position 1,
    even though "1" is numerically the smaller locant. Returns None if
    `fusion_atoms` aren't consecutive in this numbering at all."""
    low, high = sorted(numbering[atom] for atom in fusion_atoms)
    if high - low == 1:
        return low
    if (low, high) == (1, n):
        return n
    return None


def _numberings_for_pattern(graph, ring_atoms, pattern, mol):
    """Every (start, direction) numbering of this ring consistent with
    `pattern` (atomic numbers for locants 1..n) - there are 1 or 2 for
    every `_RING_TYPES` entry here (2 exactly when the ring's own pattern
    has a mirror symmetry, e.g. imidazole/pyrazole/pyridazine's two
    heteroatoms swapping under reflection); which one is "correct" depends
    on whether this ring ends up playing the base or attached role (see
    `_numbering_for_base`/`_numbering_for_attached` below), so this
    doesn't pick one itself."""
    cycle = ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    matches = []
    for start in range(n):
        for step in (1, -1):
            numbering = {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
            by_locant = sorted(numbering.items(), key=lambda kv: kv[1])
            actual = tuple(mol.GetAtomWithIdx(atom).GetAtomicNum() for atom, _ in by_locant)
            if actual == pattern:
                matches.append(numbering)
    return matches


def _numbering_for_base(numberings, fusion_atoms):
    """Among candidate numberings, pick whichever gives the *earliest bond
    letter* for `fusion_atoms` (P-25.3.1.3: "the letter as early in the
    alphabet as possible") - not the numerically lowest locant set, which
    disagrees with the earliest letter exactly for the wraparound bond
    (locants {1, n} is the *last* letter, not the first, even though 1 is
    numerically smallest). Returns None if no candidate has `fusion_atoms`
    as a single peripheral bond at all."""
    n = len(numberings[0])
    best = None
    for numbering in numberings:
        position = _bond_letter_position(numbering, fusion_atoms, n)
        if position is None:
            continue
        if best is None or position < best[1]:
            best = (numbering, position)
    return best[0] if best else None


def _numbering_for_attached(numberings, fusion_atoms):
    """Among candidate numberings, pick whichever gives the numerically
    lowest citation locants for `fusion_atoms` - the attached component's
    own locants have no "letter" of their own, so P-25.3.1.3's "as low as
    is consistent with the numbering" applies directly to the locant
    values themselves. Returns None if no candidate has `fusion_atoms` as
    a single peripheral bond at all (mirrors `_numbering_for_base`)."""
    n = len(numberings[0])
    best = None
    for numbering in numberings:
        if _bond_letter_position(numbering, fusion_atoms, n) is None:
            continue
        locants = tuple(sorted(numbering[atom] for atom in fusion_atoms))
        if best is None or locants < best[1]:
            best = (numbering, locants)
    return best[0] if best else None


def _match_ring_type(graph, ring_atoms, mol):
    """Return (pattern, name, prefix, numberings) for whichever
    `_RING_TYPES` entry matches this ring (`numberings`: every valid
    (start, direction) candidate, see `_numberings_for_pattern`), or
    None."""
    if len(ring_atoms) not in (5, 6):
        return None
    for pattern, name, prefix in _RING_TYPES:
        if len(pattern) != len(ring_atoms):
            continue
        numberings = _numberings_for_pattern(graph, ring_atoms, pattern, mol)
        if numberings:
            return pattern, name, prefix, numberings
    return None


def _ring_senior(pattern_a, size_a, pattern_b, size_b):
    """P-25.3.2.4(a)-(e): return -1 if ring a is senior (becomes the base
    component), 1 if ring b is senior, 0 if still tied after (e)."""
    hetero_a = [z for z in pattern_a if z != 6]
    hetero_b = [z for z in pattern_b if z != 6]
    best_a = min(_HETEROATOM_SENIORITY[z] for z in hetero_a)
    best_b = min(_HETEROATOM_SENIORITY[z] for z in hetero_b)
    if best_a != best_b:
        return -1 if best_a < best_b else 1
    # (b): ring count -- always 1 for the monocyclic components this step supports.
    if size_a != size_b:
        return -1 if size_a > size_b else 1
    if len(hetero_a) != len(hetero_b):
        return -1 if len(hetero_a) > len(hetero_b) else 1
    variety_a, variety_b = len(set(hetero_a)), len(set(hetero_b))
    if variety_a != variety_b:
        return -1 if variety_a > variety_b else 1
    return 0


def find_bridgehead_heteroatom_fusion_core(mol):
    """Return ((base_name, base_numbering), (attached_prefix,
    attached_numbering), fusion_atoms) if `mol` is exactly two ortho-fused
    rings from `_RING_TYPES` sharing one bridgehead heteroatom, else
    None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic() or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None

    shared = set(atom_rings[0]) & set(atom_rings[1])
    if len(shared) != 2:
        return None
    a, b = shared
    if b not in {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()}:
        return None
    shared_heteroatoms = [idx for idx in shared if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(shared_heteroatoms) != 1:
        return None

    if mol.GetNumAtoms() != len(set(atom_rings[0]) | set(atom_rings[1])):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}

    matches = []
    for ring in atom_rings:
        match = _match_ring_type(graph, ring, mol)
        if match is None:
            return None
        matches.append(match)

    (pattern_0, name_0, prefix_0, numberings_0), (pattern_1, name_1, prefix_1, numberings_1) = matches
    verdict = _ring_senior(pattern_0, len(atom_rings[0]), pattern_1, len(atom_rings[1]))
    if verdict == 0:
        return None
    if verdict < 0:
        base_name, base_numberings = name_0, numberings_0
        attached_prefix, attached_numberings = prefix_1, numberings_1
    else:
        base_name, base_numberings = name_1, numberings_1
        attached_prefix, attached_numberings = prefix_0, numberings_0

    base_numbering = _numbering_for_base(base_numberings, shared)
    attached_numbering = _numbering_for_attached(attached_numberings, shared)
    if base_numbering is None or attached_numbering is None:
        return None

    return (base_name, base_numbering), (attached_prefix, attached_numbering), shared


def has_bridgehead_heteroatom_fusion_name(mol) -> bool:
    return find_bridgehead_heteroatom_fusion_core(mol) is not None


def name_bridgehead_heteroatom_fusion(mol) -> str:
    core = find_bridgehead_heteroatom_fusion_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "this two-ring system is not a supported bridgehead-heteroatom "
            "ortho-fusion (see P-25.3.2.5.1)"
        )
    (base_name, base_numbering), (attached_prefix, attached_numbering), fusion_atoms = core

    n = len(base_numbering)
    position = _bond_letter_position(base_numbering, fusion_atoms, n)
    letter = chr(ord("a") + position - 1)

    low, high = sorted(fusion_atoms, key=lambda atom: base_numbering[atom])
    citation = f"{attached_numbering[low]},{attached_numbering[high]}"
    return f"{attached_prefix}[{citation}-{letter}]{base_name}"
