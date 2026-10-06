"""Naming of a saturated single-heteroatom monocyclic ring (P-22.2.1's
Hantzsch-Widman/retained names, `_hetero_monocyclic.py`) carrying one extra
ring C=C double bond, via the 'didehydro' prefix, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-31.2.2/P-31.2.4.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'dehydro' prefixes (multiplied 'di-'/'tetra-' etc. for an even H-atom
  count, since each extra double bond removes 2 H) cite the creation of a
  double or triple bond relative to the fully-saturated parent hydride
  named, placed immediately before the parent name. Verified against the
  primary source's own worked example (the Blue Book):
  'oxepane' (PIN) -> '2,3-didehydrooxepane'. This is the mirror-image
  mechanism of the 'hydro' prefix (which *removes* a double bond from a
  mancude parent); here a double bond is *added* to a saturated parent.
- Numbering: the heteroatom keeps its already-fixed locant 1 (same
  convention as every saturated monocyclic ring in `_hetero_monocyclic.py`
  -- there being only one heteroatom, no alternative numbering start
  exists), and of the two ring-walk directions from it, whichever gives
  the lower locant pair for the double bond wins (P-14.4's general
  lowest-locants rule -- the only choice left once the heteroatom's own
  locant is fixed), mirroring `_match_pyran_indicated_hydrogen`'s
  forward/backward locant computation.

Scope, deliberately narrow: exactly one ring, one heteroatom (O, S, or N,
matching `_hetero_monocyclic.py`'s saturated-ring table), unsubstituted,
with exactly one ring C=C double bond between two ring *carbons* -- a
double bond directly at the heteroatom (e.g. a cyclic imine) is a
different, out-of-scope shape (no spare valence question the same way,
and no real worked example checked here). Anything else returns None, so
`core.py`'s dispatch falls through unchanged.
"""

from rdkit import Chem

from ._common import adjacency, ring_cycle
from ._hetero_monocyclic import saturated_ring_name

_HETEROATOM_ELEMENTS = {"O", "S", "N"}


def find_didehydro_ring_core(mol):
    """Return (parent_name, (locant_a, locant_b)) for a `_hetero_monocyclic.py`
    saturated single-heteroatom ring with one extra ring C=C double bond
    between two ring carbons, else None (see module docstring)."""
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    size = len(ring_atoms)
    if mol.GetNumAtoms() != size:
        return None

    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 1:
        return None
    hetero = heteroatoms[0]
    hetero_atom = mol.GetAtomWithIdx(hetero)
    element = hetero_atom.GetSymbol()
    if element not in _HETEROATOM_ELEMENTS:
        return None
    if hetero_atom.GetFormalCharge() != 0 or hetero_atom.GetIsotope() != 0:
        return None
    if hetero_atom.GetIsAromatic() or hetero_atom.GetDegree() != 2:
        return None
    expected_hetero_h = 1 if element == "N" else 0
    if hetero_atom.GetTotalNumHs() != expected_hetero_h:
        return None

    carbons = [a for a in ring_atoms if a != hetero]
    for a in carbons:
        atom = mol.GetAtomWithIdx(a)
        if (
            atom.GetIsAromatic()
            or atom.GetAtomicNum() != 6
            or atom.GetDegree() != 2
            or atom.GetFormalCharge() != 0
            or atom.GetIsotope() != 0
        ):
            return None

    double_bonds = []
    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order == 1.0:
            continue
        if order != 2.0:
            return None
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if hetero in (a, b):
            return None
        double_bonds.append((a, b))
    if len(double_bonds) != 1:
        return None
    db_a, db_b = double_bonds[0]

    for a in carbons:
        expected_h = 1 if a in (db_a, db_b) else 2
        if mol.GetAtomWithIdx(a).GetTotalNumHs() != expected_h:
            return None

    name = saturated_ring_name(element, size)
    if name is None:
        return None

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, ring_atoms)
    start = ring_order.index(hetero)
    ring_order = ring_order[start:] + ring_order[:start]

    forward = {atom: i + 1 for i, atom in enumerate(ring_order)}
    backward_order = [ring_order[0]] + list(reversed(ring_order[1:]))
    backward = {atom: i + 1 for i, atom in enumerate(backward_order)}

    forward_pair = tuple(sorted((forward[db_a], forward[db_b])))
    backward_pair = tuple(sorted((backward[db_a], backward[db_b])))
    locants = min(forward_pair, backward_pair)

    return name, locants


def has_didehydro_ring_name(mol) -> bool:
    return find_didehydro_ring_core(mol) is not None


def name_didehydro_ring(mol) -> str:
    name, locants = find_didehydro_ring_core(mol)
    return f"{locants[0]},{locants[1]}-didehydro{name}"
