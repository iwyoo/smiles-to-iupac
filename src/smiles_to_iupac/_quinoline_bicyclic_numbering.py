"""Full P-25.3.3.1 whole-system peripheral numbering for the linear
(anthracene-shaped) benzo[g]quinoline/benzo[g]isoquinoline tricyclics
that `_pyridine_bicyclic_fusion.py` already names via their fusion
letter alone. Numbering starts in whichever terminal ring gives the
heteroatom its lowest locant, verified against the Blue Book's own
`10,5-[2,3]furanobenzo[g]quinoline` diagram (`tmp/bluebook/P2.pdf` p.106):
the two non-fusion middle-ring atoms land on locants 5 and 10 exactly as
drawn. The angular benzo[f]/[h] fusions are a structurally different
shape (no invariant meso-atom pair) and are deliberately out of scope.
"""

from collections import Counter, defaultdict

from rdkit import Chem

_QUINOLINE_REF = Chem.MolFromSmiles("c1ccc2ncccc2c1")
_ISOQUINOLINE_REF = Chem.MolFromSmiles("c1ccc2cnccc2c1")
_ROLE_TO_IDX = {"C5": 9, "C6": 0, "C7": 1, "C8": 2}
_REFERENCES = (_QUINOLINE_REF, _ISOQUINOLINE_REF)


def _periphery_cycle(mol, fusion_atoms, interior_bonds):
    adjacency = defaultdict(list)
    for bond in mol.GetBonds():
        if bond.GetIdx() in interior_bonds:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        adjacency[a].append(b)
        adjacency[b].append(a)
    start = next(iter(adjacency))
    cycle = [start]
    prev, cur = None, start
    while True:
        nxt = next(x for x in adjacency[cur] if x != prev)
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt
    return cycle


def _numbering_from(cycle, fusion_atoms, start_i, direction):
    n = len(cycle)
    locants = {}
    counter = 0
    for step in range(n):
        atom = cycle[(start_i + direction * step) % n]
        if atom in fusion_atoms:
            locants[atom] = f"{counter}a"
        else:
            counter += 1
            locants[atom] = str(counter)
    return locants


def _find_g_fusion_extra_atoms(mol):
    for reference in _REFERENCES:
        matches = mol.GetSubstructMatches(reference, useChirality=False)
        if len(matches) != 1:
            continue
        (match,) = matches
        idx_to_role = {v: k for k, v in _ROLE_TO_IDX.items()}
        target_to_role = {match[i]: role for i, role in idx_to_role.items()}
        extra_atoms = set(range(mol.GetNumAtoms())) - set(match)
        if len(extra_atoms) != 4:
            continue
        fusion_roles = set()
        ok = True
        for a in extra_atoms:
            for nb in mol.GetAtomWithIdx(a).GetNeighbors():
                if nb.GetIdx() in target_to_role:
                    fusion_roles.add(target_to_role[nb.GetIdx()])
                elif nb.GetIdx() not in extra_atoms:
                    ok = False
        if ok and fusion_roles == {"C6", "C7"}:
            return True
    return False


def peripheral_numbering(mol):
    """Return {atom_idx: locant_str} for a benzo[g]quinoline/isoquinoline
    shape, or None if this isn't that shape (including angular f/h fusion)."""
    if not _find_g_fusion_extra_atoms(mol):
        return None

    ring_info = mol.GetRingInfo()
    atom_rings = [set(r) for r in ring_info.AtomRings()]
    bond_count = Counter()
    for bond_ring in ring_info.BondRings():
        for bond_idx in bond_ring:
            bond_count[bond_idx] += 1
    interior_bonds = {b for b, c in bond_count.items() if c >= 2}
    fusion_atoms = set()
    for bond_idx in interior_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        fusion_atoms.add(bond.GetBeginAtomIdx())
        fusion_atoms.add(bond.GetEndAtomIdx())

    terminal_rings = [
        (ring & fusion_atoms, ring - fusion_atoms)
        for ring in atom_rings
        if len(ring & fusion_atoms) == 2 and len(ring - fusion_atoms) == 4
    ]
    if len(terminal_rings) != 2:
        return None

    cycle = _periphery_cycle(mol, fusion_atoms, interior_bonds)
    idx_of = {a: i for i, a in enumerate(cycle)}
    hetero = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 6]
    if len(hetero) != 1:
        return None

    candidates = []
    for fus, nonfus in terminal_rings:
        for f_atom in fus:
            i = idx_of[f_atom]
            for step in (1, -1):
                neighbor = cycle[(i + step) % len(cycle)]
                if neighbor in nonfus:
                    candidates.append(_numbering_from(cycle, fusion_atoms, (i + step) % len(cycle), step))

    return min(candidates, key=lambda locants: int(locants[hetero[0]]))
