"""Naming of the pyridinone tautomer -- a 6-membered, single-nitrogen
mancude ring bearing one ring-carbon oxo group -- per the IUPAC 2013
Recommendations:

- P-31.1.4.3.4 (indicated hydrogen for a mancude ring's oxo tautomer):
  the ring nitrogen's own hydrogen is cited as front indicated hydrogen
  ('1H-', the nitrogen always at locant 1 for this single-heteroatom
  ring, P-22.2.2.1.2), and the oxo group is cited as an ordinary '-one'
  suffix at its own ring locant -- e.g. '1H-pyridin-2-one', not
  'pyridin-2(1H)-one' (confirmed as PubChem's own form for CID 8871/
  12290). This is the same tautomer-locant mechanism `_hetero_monocyclic
  .py`'s already-shipped `name_pyran_indicated_hydrogen` (2H-/4H-pyran)
  uses, with a real oxo substituent standing in for the second hydrogen.
- Ring traversal direction is chosen to give the oxo group the lowest
  locant, mirroring `_match_pyran_indicated_hydrogen`'s own direction
  choice.
- Confirmed against two real structures: `O=c1cccc[nH]1` -> '1H-pyridin-
  2-one' (PubChem CID 8871) and `O=c1cc[nH]cc1` -> '1H-pyridin-4-one'
  (CID 12290).

Out of scope: a second ring heteroatom (the pyrimidinedione shape needed
for uracil/thymine/cytosine), any other exocyclic substituent, and any
non-6-membered or non-mancude ring.
"""

from ._common import adjacency, ring_cycle


def _match_pyridinone(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    if len(ring_atoms) != 6:
        return None
    ring_set = set(ring_atoms)
    if any(not mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms):
        return None

    nitrogens = [atom for atom in ring_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() == 7]
    if len(nitrogens) != 1:
        return None
    (n_atom,) = nitrogens
    n_data = mol.GetAtomWithIdx(n_atom)
    if n_data.GetTotalNumHs() != 1 or n_data.GetFormalCharge() != 0:
        return None
    if any(n.GetIdx() not in ring_set for n in n_data.GetNeighbors()):
        return None

    oxo_atoms = []
    for atom in ring_atoms:
        if atom == n_atom:
            continue
        data = mol.GetAtomWithIdx(atom)
        if data.GetAtomicNum() != 6:
            return None
        exocyclic = [n for n in data.GetNeighbors() if n.GetIdx() not in ring_set]
        if not exocyclic:
            if data.GetTotalNumHs() != 1:
                return None
            continue
        if len(exocyclic) != 1 or data.GetTotalNumHs() != 0:
            return None
        (oxygen,) = exocyclic
        if oxygen.GetAtomicNum() != 8 or oxygen.GetFormalCharge() != 0 or oxygen.GetDegree() != 1:
            return None
        if mol.GetBondBetweenAtoms(atom, oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            return None
        oxo_atoms.append(atom)
    if len(oxo_atoms) != 1:
        return None
    (oxo_atom,) = oxo_atoms

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, ring_atoms)
    start = ring_order.index(n_atom)
    ring_order = ring_order[start:] + ring_order[:start]
    forward_locant = ring_order.index(oxo_atom) + 1
    backward_locant = len(ring_order) - forward_locant + 2
    return min(forward_locant, backward_locant)


def has_pyridinone_shape(mol) -> bool:
    return _match_pyridinone(mol) is not None


def name_pyridinone(mol) -> str:
    locant = _match_pyridinone(mol)
    return f"1H-pyridin-{locant}-one"
