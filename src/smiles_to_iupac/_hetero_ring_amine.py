"""Naming of a plain, saturated, single-heteroatom-nitrogen monocyclic
ring (piperidine/pyrrolidine/azetidine/aziridine/azepane, P-22.2.1) whose
own nitrogen carries no substituent (a plain ring N-H) but exactly one
ring CARBON bears a single primary exocyclic -NH2, cited as the '-amine'
suffix (P-33.1) the same way `_amine.py`'s `_name_cyclic_amine` cites it
on a plain carbocyclic ring -- e.g. 'piperidin-4-amine' (PubChem CID
424361), 'piperidin-3-amine' (CID 148119), 'piperidin-2-amine' (CID
421842), 'pyrrolidin-2-amine' (CID 14298876), 'pyrrolidin-3-amine' (CID
164401), 'azetidin-2-amine' (CID 21661087), 'aziridin-2-amine' (CID
20071107), 'azepan-2-amine' (CID 12040473), 'azepan-3-amine' (CID
2756445) -- all PubChem structure matches, one worked example per ring
size (3-7) and, for piperidine/pyrrolidine, two locants each.

`_amine.py` explicitly rejects any ring bearing a nitrogen heteroatom
(see that module's own scope docstring -- the plain ring N-H here still
looks like a "secondary amine nitrogen on a ring" to its naive
neighbor-count model, and the coexisting exocyclic amine trips its
multi-nitrogen guard); `_ring_amine.py` handles the opposite shape (the
ring nitrogen itself carries a substituent, no exocyclic suffix group).
This module fills the gap in between: ring nitrogen plain, exocyclic
-NH2 suffix on a ring carbon.

The ring's own nitrogen is always locant 1 (Hantzsch-Widman/retained-name
numbering, P-22.2.1); only the traversal *direction* around the ring is
free, chosen (P-31.1.4.2.4, same "lowest locant" tie-break as any other
ring) to give the exocyclic amine's own locant the lower of its two
possible values.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one exocyclic amine, or any other substituent (including a
  halogen) anywhere on the ring -- separate follow-ups.
- The ring nitrogen bearing a substituent -- `_ring_amine.py`'s shape.
- A two-heteroatom ring (the morpholine/piperazine/thiomorpholine
  siblings), ring unsaturation, or a ring size outside 3-7 --
  `saturated_ring_name`'s own scope.
- A secondary/tertiary exocyclic amine (the -NH2 substituent's own
  nitrogen bearing more than one carbon), or a specified stereocenter --
  mirror `_amine.py`'s identical restrictions, deferred as separate
  follow-ups.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, ring_cycle, specified_stereocenters
from ._hetero_monocyclic import saturated_ring_name

_RING_SIZES = (3, 4, 5, 6, 7)


def _hetero_ring_amine_shape(mol):
    """(ring_atoms, ring_n, amine_carbon) if `mol` is a single plain
    saturated N-heterocycle (unsubstituted ring N-H) with exactly one
    exocyclic primary amine on a ring carbon and nothing else; else
    None."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    if len(ring_atoms) not in _RING_SIZES:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None

    ring_set = set(ring_atoms)
    ring_n_atoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    if len(ring_n_atoms) != 1:
        return None
    (ring_n,) = ring_n_atoms
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring_atoms if a != ring_n):
        return None

    for atom in ring_atoms:
        rd_atom = mol.GetAtomWithIdx(atom)
        if rd_atom.GetFormalCharge() != 0 or rd_atom.GetIsotope() != 0:
            return None

    ring_n_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(ring_n).GetNeighbors()]
    if len(ring_n_neighbors) != 2 or any(n not in ring_set for n in ring_n_neighbors):
        return None

    exo_amines = []
    for atom in ring_atoms:
        if atom == ring_n:
            continue
        exo = [n for n in mol.GetAtomWithIdx(atom).GetNeighbors() if n.GetIdx() not in ring_set]
        if not exo:
            continue
        if len(exo) != 1:
            return None
        (exo_atom,) = exo
        if exo_atom.GetAtomicNum() != 7:
            return None
        if exo_atom.GetFormalCharge() != 0 or exo_atom.GetIsotope() != 0:
            return None
        if len(list(exo_atom.GetNeighbors())) != 1:
            return None
        exo_amines.append(atom)

    if len(exo_amines) != 1:
        return None
    (amine_carbon,) = exo_amines

    if any(bond.GetBondTypeAsDouble() != 1.0 for bond in mol.GetBonds()):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    return ring_atoms, ring_n, amine_carbon


def has_hetero_ring_amine_shape(mol) -> bool:
    return _hetero_ring_amine_shape(mol) is not None


def name_hetero_ring_amine(mol) -> str:
    shape = _hetero_ring_amine_shape(mol)
    if shape is None:
        raise UnsupportedStructure(
            "not a plain saturated single-nitrogen heterocycle with "
            "exactly one exocyclic primary amine on a ring carbon"
        )
    ring_atoms, ring_n, amine_carbon = shape
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a hetero-ring amine is "
            "not supported yet"
        )

    parent = saturated_ring_name("N", len(ring_atoms))
    if parent is None:
        raise UnsupportedStructure(
            "this ring size has no saturated single-nitrogen "
            "Hantzsch-Widman/retained name (P-22.2.1)"
        )

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, ring_atoms)
    start = ring_order.index(ring_n)
    rotated = ring_order[start:] + ring_order[:start]
    # Reversing direction must keep the ring nitrogen fixed at locant 1
    # (Hantzsch-Widman numbering never starts anywhere else) -- reversing
    # `rotated` outright would instead move it to the last position.
    other_direction = [rotated[0]] + list(reversed(rotated[1:]))

    best_locant = None
    for candidate in (rotated, other_direction):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        locant = position_of[amine_carbon]
        if best_locant is None or locant < best_locant:
            best_locant = locant

    stem = parent[:-1] if parent.endswith("e") else parent
    return f"{stem}-{best_locant}-amine"
