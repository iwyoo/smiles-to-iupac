"""Naming of the phenanthroline/naphthyridine retained diaza names, per
P-2 Table 2.8 (https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): phenanthrene/
naphthalene with two ring CH positions replaced by pyridine-type N, cited
with the N locants -- not generic replacement ("aza") nomenclature.
Naphthyridine reuses naphthalene's own numbering; phenanthroline uses its
own traditional numbering, a fixed rotation of phenanthrene's own walk
(confirmed against PubChem CIDs 67473/1318/67472/72812).
"""

from ._aromatic import (
    _classify_shape,
    _phenanthrene_candidates,
    _ring_adjacency,
    _ring_path_order,
    _straight_chain_candidates,
)
from ._common import UnsupportedStructure, adjacency

_NAME_BY_PARENT = {"naphthalene": "naphthyridine", "phenanthrene": "phenanthroline"}
_BENZODIAZINES = {(1, 2): "cinnoline", (1, 3): "quinazoline", (1, 4): "quinoxaline", (2, 3): "phthalazine"}


def find_phenanthroline_naphthyridine_core(mol):
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) not in (2, 3):
        return None
    if mol.GetNumAtoms() != (10 if len(atom_rings) == 2 else 14):
        return None
    for ring in atom_rings:
        if len(ring) != 6:
            return None

    nitrogens = []
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic():
                return None
            if atom.GetAtomicNum() == 7:
                if atom.GetTotalNumHs() != 0 or atom.GetFormalCharge() != 0 or atom.GetDegree() != 2:
                    return None
                nitrogens.append(idx)
            elif atom.GetAtomicNum() != 6:
                return None
    if len(set(nitrogens)) != 2:
        return None

    ring_atom_sets = [set(r) for r in atom_rings]
    bond_ring_count = {}
    for bonds in ring_info.BondRings():
        for b in bonds:
            bond_ring_count[b] = bond_ring_count.get(b, 0) + 1
    fusion_bond_idxs = {b for b, c in bond_ring_count.items() if c >= 2}

    try:
        adj, fusion_bonds_by_pair = _ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
        n = len(atom_rings)
        ring_order = _ring_path_order(adj, n)
        graph = adjacency(mol)
        shape = () if n == 2 else _classify_shape(graph, atom_rings, ring_order, fusion_bonds_by_pair)
    except UnsupportedStructure:
        return None

    parent = {(2, ()): "naphthalene", (3, ("bent",)): "phenanthrene"}.get((n, shape))
    if parent is None:
        return None

    return parent, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order, fusion_bonds_by_pair, tuple(nitrogens)


def name_phenanthroline_naphthyridine(mol, core) -> str:
    parent, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order, fusion_bonds_by_pair, nitrogens = core
    graph = adjacency(mol)

    if parent == "phenanthrene":
        candidates = [
            {atom: ((locant - 1 + 6) % 10) + 1 for atom, locant in c.items()}
            for c in _phenanthrene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order)
        ]
    else:
        candidates = _straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order)

    ring_of = [next(i for i, ring in enumerate(ring_atom_sets) if n in ring) for n in nitrogens]
    same_ring = ring_of[0] == ring_of[1]
    if parent == "phenanthrene" and (same_ring or set(ring_of) != {ring_order[0], ring_order[-1]}):
        raise UnsupportedStructure(
            "a diazaphenanthrene without one nitrogen in each terminal ring is not a phenanthroline"
        )

    best = None
    for locants in candidates:
        pair = tuple(sorted(locants[a] for a in nitrogens))
        if best is None or pair < best:
            best = pair

    if parent == "naphthalene" and same_ring:
        # Both nitrogens in one ring: the benzodiazines keep their own retained names (P-25.2.2.4).
        retained = _BENZODIAZINES.get(best)
        if retained is None:
            raise UnsupportedStructure("this benzodiazine locant pattern has no retained name")
        return retained
    return f"{best[0]},{best[1]}-{_NAME_BY_PARENT[parent]}"
