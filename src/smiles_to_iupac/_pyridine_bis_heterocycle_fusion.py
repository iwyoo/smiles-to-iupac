"""Fusion naming (P-25.3.6.1's "identical attached components") for a
three-ring system: pyridine as the sole (senior) parent component, with
two *identical* five-membered heteromonocycles (furan or thiophene) as
first-order attached components at two separate bonds, per the IUPAC 2013
Recommendations ("the Blue Book").

Unlike `_benzo_bis_heterocycle_fusion.py`'s multiparent mechanism (#604),
no multiparent construction is needed here: pyridine's nitrogen already
outranks furan/thiophene's O/S outright per P-25.3.2.4(a), so pyridine
stays the single parent and the "di" multiplying prefix goes on the
*attached* component instead (e.g. "dithieno[2,3-b:3',2'-e]pyridine",
contrast with #604's "benzo[...]dithiophene", where "di" was on the
*parent* because the bridging ring there wasn't senior to its neighbors).

Two numbering problems, reusing existing machinery:

- Each attached ring's own citation locants: `_local_numbering`, reused
  unchanged from `_pyridine_heterocycle_fusion.py`.
- Pyridine's own numbering must satisfy *two* simultaneous fusion bonds at
  once, like `_benzo_bis_heterocycle_fusion.py`'s `_benzo_numbering` - but
  pyridine's nitrogen is fixed at locant 1, so there are only 2 candidate
  numberings (one per ring-walk direction) instead of benzo's 12. For
  each, both fusion bonds must avoid touching the nitrogen (reusing
  `_base_bond_letter`'s own check), then the bond with the earlier letter
  is cited first (unprimed), the other second (primed) - `_two_bond_letters`
  picks whichever of the 2 directions makes the first-cited letter as
  early as possible.

Real registered example confirming this family: dithieno[2,3-b:3',2'-e]
pyridine (PubChem CID 129867025) - pyridine as parent (N outranks S
outright, no multiparent needed at all), two identical thiophenes as
first-order attached components at bonds 'b' and 'e'. Only one real
example was found for this exact shape despite a reasonable search - the
generalization here is the *mechanism* (extending a fixed-locant-1 ring's
numbering search to two simultaneous fusion bonds, reusing the
per-connection letter/locant computation unchanged), not a large
named-molecule sample.

Two related shapes checked and explicitly left out of scope (see the
GitHub issue this module closes for the detail): two *different* first-
order attached components (P-25.3.4.2.3.1's own worked example wasn't
found registered on PubChem, so not used as a verification anchor here),
and true second-order attachment (a component fused onto another attached
component, not the parent - the Blue Book's own simplest example needs a
retained bicyclic name, a separate prerequisite).

Scope, deliberately narrow, matching the sibling fusion modules:
- Exactly three rings: pyridine (parent) plus two identical five-membered
  heteromonocycles (both furan, or both thiophene), each ortho-fused to
  pyridine at a separate bond not touching pyridine's own nitrogen, not
  fused to each other.
- Each heteromonocycle's own fusion bond must touch its own heteroatom's
  neighbor (the same `_local_numbering` restriction the sibling modules
  already carry).
- Completely unsubstituted, no extra atoms, single molecular fragment.

Anything else (different attached rings, 3+ occurrences, attachment at
pyridine's own nitrogen, a heteroatom other than O/S, any substituent)
raises `UnsupportedStructure` and falls through to other dispatch
branches in core.py, exactly like every other retained/computed-name
module in this project.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, ring_cycle
from ._pyridine_heterocycle_fusion import _local_numbering

_ATTACHED_RING_NAMES = {
    8: "furo",
    16: "thieno",
}


def _bond_letter_position(numbering, fusion_atoms, n):
    """Same convention as the sibling fusion modules: the wraparound bond
    (locants {1, n}) is position n (the last letter), not position 1."""
    low, high = sorted(numbering[atom] for atom in fusion_atoms)
    if high - low == 1:
        return low
    if (low, high) == (1, n):
        return n
    return None


def _two_bond_letters(graph, pyridine_ring, nitrogen, fusion_pairs):
    """Pyridine's nitrogen is fixed at locant 1, leaving only 2 candidate
    numberings (one per ring-walk direction). For each, both fusion bonds
    must avoid touching the nitrogen; the bond with the earlier letter is
    cited first (unprimed), the other second (primed). Returns
    (numbering, first_pair, first_letter, second_pair, second_letter), or
    None if no direction makes both bonds simple peripheral bonds clear
    of the nitrogen."""
    cycle = ring_cycle(graph, list(pyridine_ring))
    n = len(cycle)
    start = cycle.index(nitrogen)
    best = None
    for step in (1, -1):
        numbering = {cycle[(start + step * k) % n]: k + 1 for k in range(n)}
        if numbering[nitrogen] in (numbering[a] for pair in fusion_pairs for a in pair):
            continue
        positions = []
        ok = True
        for pair in fusion_pairs:
            position = _bond_letter_position(numbering, pair, n)
            if position is None:
                ok = False
                break
            positions.append(position)
        if not ok:
            continue
        order = sorted(range(2), key=lambda i: positions[i])
        first_position = positions[order[0]]
        if best is None or first_position < best[0]:
            best = (first_position, numbering, order, positions)
    if best is None:
        return None
    _first_position, numbering, order, positions = best
    first_pair, second_pair = fusion_pairs[order[0]], fusion_pairs[order[1]]
    first_letter = chr(ord("a") + positions[order[0]] - 1)
    second_letter = chr(ord("a") + positions[order[1]] - 1)
    return numbering, first_pair, first_letter, second_pair, second_letter


def find_pyridine_bis_heterocycle_fusion_core(mol):
    """Return (pyridine_ring, nitrogen, element, [(ring, heteroatom,
    fusion_atoms), (ring, heteroatom, fusion_atoms)]) if `mol` is exactly
    pyridine ortho-fused to two identical furan/thiophene rings at
    separate bonds, else None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 3:
        return None
    sizes = sorted(len(r) for r in atom_rings)
    if sizes != [5, 5, 6]:
        return None
    pyridine_ring = next(r for r in atom_rings if len(r) == 6)
    attached_rings = [r for r in atom_rings if len(r) == 5]

    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic() or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None

    pyridine_hetero = [idx for idx in pyridine_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(pyridine_hetero) != 1 or mol.GetAtomWithIdx(pyridine_hetero[0]).GetAtomicNum() != 7:
        return None
    nitrogen = pyridine_hetero[0]

    heteroatoms = []
    for ring in attached_rings:
        hetero = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        if len(hetero) != 1 or mol.GetAtomWithIdx(hetero[0]).GetAtomicNum() not in _ATTACHED_RING_NAMES:
            return None
        heteroatoms.append(hetero[0])
    if mol.GetAtomWithIdx(heteroatoms[0]).GetAtomicNum() != mol.GetAtomWithIdx(heteroatoms[1]).GetAtomicNum():
        return None

    if set(attached_rings[0]) & set(attached_rings[1]):
        return None

    fusion_pairs = []
    for ring in attached_rings:
        shared = set(pyridine_ring) & set(ring)
        if len(shared) != 2:
            return None
        a, b = shared
        if b not in {n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()}:
            return None
        fusion_pairs.append(shared)
    if fusion_pairs[0] & fusion_pairs[1]:
        return None

    if mol.GetNumAtoms() != len(set(pyridine_ring) | set(attached_rings[0]) | set(attached_rings[1])):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    element = mol.GetAtomWithIdx(heteroatoms[0]).GetAtomicNum()
    connections = list(zip(attached_rings, heteroatoms, fusion_pairs))
    return pyridine_ring, nitrogen, element, connections


def has_pyridine_bis_heterocycle_fusion_name(mol) -> bool:
    return find_pyridine_bis_heterocycle_fusion_core(mol) is not None


def name_pyridine_bis_heterocycle_fusion(mol) -> str:
    core = find_pyridine_bis_heterocycle_fusion_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "this three-ring system is not a supported pyridine + two "
            "identical five-membered heterocycle ortho-fusion (see "
            "P-25.3.6.1)"
        )
    pyridine_ring, nitrogen, element, connections = core
    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}

    fusion_pairs = [pair for _ring, _hetero, pair in connections]
    result = _two_bond_letters(graph, pyridine_ring, nitrogen, fusion_pairs)
    if result is None:
        raise UnsupportedStructure(
            "fusion at pyridine's own nitrogen is not supported (see P-25.3.1.3)"
        )
    pyridine_numbering, first_pair, first_letter, second_pair, second_letter = result
    by_pair = {frozenset(pair): (ring, hetero) for ring, hetero, pair in connections}

    parts = []
    for pair, letter, primed in ((first_pair, first_letter, False), (second_pair, second_letter, True)):
        ring, hetero = by_pair[frozenset(pair)]
        attached_numbering = _local_numbering(graph, ring, hetero, pair)
        if attached_numbering is None:
            raise UnsupportedStructure(
                "the fusion bond does not touch an attached ring's own "
                "heteroatom (only the 'b'-lettered fusion is supported, "
                "see P-25.3.1.3)"
            )
        low, high = sorted(pair, key=lambda atom: pyridine_numbering[atom])
        # Unlike `_benzo_bis_heterocycle_fusion.py` (#604), where the
        # *parent* is the doubled component and its own letters get
        # primed, here the *attached* component is doubled - so its own
        # numbers get primed instead, not the (singular) parent's letter.
        prime = "'" if primed else ""
        parts.append(
            f"{attached_numbering[low]}{prime},{attached_numbering[high]}{prime}-{letter}"
        )

    attached_prefix = _ATTACHED_RING_NAMES[element]
    return f"di{attached_prefix}[{parts[0]}:{parts[1]}]pyridine"
