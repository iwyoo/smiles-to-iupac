"""Fusion-locant-letter naming (P-25.3.1.3) for a plain cyclopenta ring
ortho-fused onto naphthalene at a bond not touching its own 4a/8a ring-
fusion carbons, plus P-25.3.3.1's whole-system numbering for the resulting
indicated hydrogen. Unlike `_pyridine_bicyclic_fusion.py`'s benzo-fusion
case, the attached 5-ring is never fully aromatic, so its numbering (and
therefore the indicated-hydrogen locant) isn't fixed by the base alone:
when one of the 5-ring's two fusion atoms is itself peripherally adjacent
to a third (naphthalene-internal) fusion atom, numbering starts next to
that atom (matching literature numbering: 1 is adjacent to '9b', not
'3a'); otherwise the two fusion atoms are symmetric and the tie-break is
the lowest locant to the indicated-hydrogen atom itself (P-25.3.3.1.2(f)).
Verified against PubChem CIDs 11745004 (1H-, letter a), 11105721 (3H-,
letter a), 6451436 (1H-, letter b).
"""

from collections import Counter, defaultdict

from rdkit import Chem

from ._common import UnsupportedStructure

_NAPHTHALENE_REF = Chem.MolFromSmiles("c1ccc2ccccc2c1")
# Role -> ref atom index, derived from the reference SMILES's own atom
# order (see module history/PR description for the periphery walk).
_ROLE_TO_IDX = {"C1": 9, "C2": 0, "C3": 1, "C4": 2, "C4a": 3, "C5": 4, "C6": 5, "C7": 6, "C8": 7, "C8a": 8}
_PERIPHERY_ROLES = ["C1", "C2", "C3", "C4", "C4a", "C5", "C6", "C7", "C8", "C8a"]
_LETTERS = "abcdefghij"


def _periphery_cycle(mol):
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
    return cycle, fusion_atoms, atom_rings


def _find_fusion_letter(mol):
    matches = mol.GetSubstructMatches(_NAPHTHALENE_REF, useChirality=False)
    if len(matches) != 1:
        return None, None
    (match,) = matches
    idx_to_role = {v: k for k, v in _ROLE_TO_IDX.items()}
    target_to_role = {match[i]: role for i, role in idx_to_role.items()}
    core_atoms = set(match)
    extra_atoms = set(range(mol.GetNumAtoms())) - core_atoms
    if len(extra_atoms) != 3:
        return None, None
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in extra_atoms):
        return None, None

    fusion_roles = set()
    for a in extra_atoms:
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            if n.GetIdx() in target_to_role:
                fusion_roles.add(target_to_role[n.GetIdx()])
            elif n.GetIdx() not in extra_atoms:
                return None, None
    if len(fusion_roles) != 2:
        return None, None
    i1 = _PERIPHERY_ROLES.index(next(iter(fusion_roles)))
    i2 = _PERIPHERY_ROLES.index(next(r for r in fusion_roles if _PERIPHERY_ROLES.index(r) != i1))
    lo, hi = sorted((i1, i2))
    if (lo, hi) != (0, 9) and hi - lo != 1:
        return None, None
    letter_index = 9 if (lo, hi) == (0, 9) else lo
    if _PERIPHERY_ROLES[letter_index].endswith("a"):
        return None, None
    return _LETTERS[letter_index], extra_atoms


def _indicated_hydrogen_locant(mol, fusion_atoms_naphthalene, cyclopenta_atoms):
    cycle, fusion_atoms, atom_rings = _periphery_cycle(mol)
    n = len(cycle)
    five_ring = next((r for r in atom_rings if len(r) == 5 and r & cyclopenta_atoms), None)
    if five_ring is None:
        return None
    fusion5 = five_ring & fusion_atoms
    if len(fusion5) != 2:
        return None
    idx_of = {atom: i for i, atom in enumerate(cycle)}

    def numbering_from(start_i, direction):
        locants = {}
        counter = 0
        for step in range(n):
            atom = cycle[(start_i + direction * step) % n]
            if atom in fusion_atoms:
                locants[atom] = counter
            else:
                counter += 1
                locants[atom] = counter
        return locants

    def start_away_from(f_atom):
        i = idx_of[f_atom]
        for step in (1, -1):
            neighbor = cycle[(i + step) % n]
            if neighbor in five_ring and neighbor not in fusion5:
                return (i + step) % n, step
        return None

    sp3_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() in five_ring and a.GetTotalNumHs() == 2]
    if len(sp3_atoms) != 1:
        return None
    sp3 = sp3_atoms[0]

    f_doubles = []
    for f in fusion5:
        i = idx_of[f]
        for step in (1, -1):
            neighbor = cycle[(i + step) % n]
            if neighbor in fusion_atoms and neighbor not in five_ring:
                f_doubles.append(f)
                break

    if len(f_doubles) == 1:
        start_i, direction = start_away_from(f_doubles[0])
        return numbering_from(start_i, direction)[sp3]
    candidates = [start_away_from(f) for f in fusion5]
    candidates = [c for c in candidates if c is not None]
    if not candidates:
        return None
    return min(numbering_from(s, d)[sp3] for s, d in candidates)


def _find_core(mol):
    if mol.GetNumAtoms() != 13:
        return None
    if mol.GetRingInfo().NumRings() != 3:
        return None
    if any(atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()):
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    letter, cyclopenta_atoms = _find_fusion_letter(mol)
    if letter is None:
        return None
    naphthalene_fusion_atoms = set(range(mol.GetNumAtoms())) - cyclopenta_atoms
    indicated_h = _indicated_hydrogen_locant(mol, naphthalene_fusion_atoms, cyclopenta_atoms)
    if indicated_h is None:
        return None
    return f"{indicated_h}H-cyclopenta[{letter}]naphthalene"


def has_cyclopenta_naphthalene_name(mol) -> bool:
    return _find_core(mol) is not None


def name_cyclopenta_naphthalene(mol) -> str:
    name = _find_core(mol)
    if name is None:
        raise UnsupportedStructure(
            "this tricyclic system is not a supported cyclopenta-fused "
            "naphthalene shape (see P-25.3.1.3/P-25.3.3.1)"
        )
    return name
