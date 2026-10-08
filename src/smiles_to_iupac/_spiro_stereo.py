"""Spiro atoms of the type Xabab (P-93.5.3.2): RDKit drops the tetrahedral mark of a spiro atom whose two rings are
constitutionally identical, so it is read from the raw SMILES parse.

The ligands form two equivalent pairs a/a' and b/b' with a > b. The analysis starts at the lowest-numbered ring, so a
is a ligand of that ring and b the other-class ligand in the same ring; the sequence a > a' > b > b' gives R or S.
Swapping both pairs is an even permutation, so the descriptor does not depend on which equivalent ligand is taken as a.
"""

import re

from rdkit import Chem

from ._axial_stereo import _merge

_DEPTH = 8
_SPIRO_LOCANT = re.compile(r"(?<!di)(?<!tri)(?<!tetra)spiro\[(\d+)\.(\d+)\]")


def _signature(mol, root, center):
    """Atomic numbers of the hierarchical digraph of the branch at `root`, sphere by sphere; a multiple bond is
    duplicated and a ring closure ends the branch with a duplicate atom."""
    nodes = [(root, center, frozenset({center, root}), mol.GetAtomWithIdx(root).GetAtomicNum())]
    spheres = []
    for _ in range(_DEPTH):
        spheres.append(tuple(sorted((z for *_, z in nodes), reverse=True)))
        following = []
        for atom_index, parent, path, _z in nodes:
            if atom_index is None:
                continue
            atom = mol.GetAtomWithIdx(atom_index)
            for bond in atom.GetBonds():
                other = bond.GetOtherAtomIdx(atom_index)
                if other == parent:
                    continue
                z = mol.GetAtomWithIdx(other).GetAtomicNum()
                if other in path:
                    following.append((None, atom_index, path, z))
                else:
                    following.append((other, atom_index, path | {other}, z))
                for _ in range(max(round(bond.GetBondTypeAsDouble()) - 1, 0)):
                    following.append((None, atom_index, path, z))
            following.extend((None, atom_index, path, 1) for _ in range(atom.GetTotalNumHs()))
        nodes = following
        if not nodes:
            break
    return spheres


def _compare(first, second):
    for a, b in zip(first, second):
        if a != b:
            return 1 if a > b else -1
    return (len(first) > len(second)) - (len(first) < len(second))


def _classes(mol, center):
    neighbours = [n.GetIdx() for n in center.GetNeighbors()]
    signatures = {n: _signature(mol, n, center.GetIdx()) for n in neighbours}
    groups = []
    for n in neighbours:
        for group in groups:
            if _compare(signatures[group[0]], signatures[n]) == 0:
                group.append(n)
                break
        else:
            groups.append([n])
    if sorted(len(g) for g in groups) != [2, 2]:
        return None
    first, second = groups
    return (first, second) if _compare(signatures[first[0]], signatures[second[0]]) > 0 else (second, first)


def _parity(sequence, target):
    order = {atom: i for i, atom in enumerate(sequence)}
    permutation = [order[atom] for atom in target]
    swaps = 0
    for i in range(len(permutation)):
        while permutation[i] != i:
            j = permutation[i]
            permutation[i], permutation[j] = permutation[j], permutation[i]
            swaps += 1
    return swaps % 2


def _descriptor(raw, center):
    classes = _classes(raw, center)
    if classes is None:
        return None
    high, low = classes
    rings = [set(r) for r in raw.GetRingInfo().AtomRings() if center.GetIdx() in r]
    ring = next((r for r in rings if high[0] in r), None)
    partner = next((n for n in low if ring is not None and n in ring), None)
    if partner is None:
        return None
    other_low = next(n for n in low if n != partner)
    priority = [high[0], high[1], partner, other_low]
    sequence = [n.GetIdx() for n in center.GetNeighbors()]
    tag = center.GetChiralTag()
    if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW):
        return None
    counterclockwise = (tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW) != bool(
        _parity(sequence, [priority[3], priority[0], priority[1], priority[2]])
    )
    return "R" if counterclockwise else "S"


def cite_spiro_stereo(smiles, name):
    if "@" not in smiles or name.count("spiro") != 1:
        return name
    locants = _SPIRO_LOCANT.search(name)
    if locants is None:
        return name
    raw = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = Chem.MolFromSmiles(smiles)
    if raw is None or mol is None or raw.GetNumAtoms() != mol.GetNumAtoms():
        return name
    raw.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(raw)
    seen = {e.centeredOn for e in Chem.FindPotentialStereo(mol)}
    for atom in raw.GetAtoms():
        if (
            atom.GetIdx() in seen
            or atom.GetChiralTag() in (Chem.ChiralType.CHI_UNSPECIFIED, Chem.ChiralType.CHI_OTHER)
            or atom.GetAtomicNum() != 6
            or atom.GetDegree() != 4
            or raw.GetRingInfo().NumAtomRings(atom.GetIdx()) != 2
        ):
            continue
        code = _descriptor(raw, atom)
        if code is not None:
            return _merge(name, f"{int(locants.group(1)) + 1}{code}")
    return name
