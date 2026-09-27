"""Naming of the pyrimidinone tautomer -- a 6-membered mancude ring with
two nitrogens in a 1,3 relationship, one bearing the ring's indicated
hydrogen, the other a plain unsubstituted ring nitrogen, plus one
ring-carbon oxo group -- per the IUPAC 2013 Recommendations:

- P-31.1.4.3.4, the same indicated-hydrogen mechanism `_pyridinone.py`
  already ships for the single-nitrogen case (P-22.2.2.1.2 fixes the
  indicated-hydrogen nitrogen at locant 1); here the second nitrogen's own
  lowest-locant-set requirement (P-22.2.2.1.2: heteroatoms as a set get
  the lowest locants, {1,3} beating {1,5}) fixes the ring traversal
  direction outright, leaving no separate freedom to also optimize the
  oxo group's own locant the way `_pyridinone.py`'s single-heteroatom
  case does.
- Confirmed against a real structure: `O=c1nccc[nH]1` -> '1H-pyrimidin-
  2-one' (PubChem CID 68401).

Out of scope: a second oxo group (the uracil/thymine diketo shape, which
needs the Blue Book's own parenthesized multi-locant added-hydrogen
convention -- a different, harder construction), any other exocyclic
substituent, and any non-1,3 relationship between the two nitrogens.
"""

from ._common import adjacency, ring_cycle


def _match_pyrimidinone(mol):
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
    if len(nitrogens) != 2:
        return None
    nh_candidates = [n for n in nitrogens if mol.GetAtomWithIdx(n).GetTotalNumHs() == 1]
    plain_candidates = [n for n in nitrogens if mol.GetAtomWithIdx(n).GetTotalNumHs() == 0]
    if len(nh_candidates) != 1 or len(plain_candidates) != 1:
        return None
    (n_h_atom,) = nh_candidates
    (n_plain_atom,) = plain_candidates
    for n_atom in nitrogens:
        n_data = mol.GetAtomWithIdx(n_atom)
        if n_data.GetFormalCharge() != 0:
            return None
        if any(n.GetIdx() not in ring_set for n in n_data.GetNeighbors()):
            return None

    oxo_atoms = []
    for atom in ring_atoms:
        if atom in nitrogens:
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
    start = ring_order.index(n_h_atom)
    ring_order = ring_order[start:] + ring_order[:start]
    other_n_locant = ring_order.index(n_plain_atom) + 1
    if other_n_locant not in (3, 5):
        return None
    if other_n_locant == 5:
        ring_order = [ring_order[0]] + list(reversed(ring_order[1:]))
    return ring_order.index(oxo_atom) + 1


def has_pyrimidinone_shape(mol) -> bool:
    return _match_pyrimidinone(mol) is not None


def name_pyrimidinone(mol) -> str:
    locant = _match_pyrimidinone(mol)
    return f"1H-pyrimidin-{locant}-one"
