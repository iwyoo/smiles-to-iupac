"""Substituted steroids named on the retained steroid parent hydrides (P-101.2, P-101.4): the 1989 steroid numbering
fixes every locant, unsaturation is cited as 'ene' endings (estra-1,3,5(10)-triene), groups as suffixes or prefixes and
configuration with α/β/ξ locants (cholest-5-en-3β-ol, 17β-hydroxy-5α-androstan-3-one). Configuration is read from a 3D
embedding: β is the face from which the numbering of ring A runs anticlockwise."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from ._common import UnsupportedStructure, adjacency, halogen_substituents, multiplied_word
from ._substituents import format_substituent_prefixes, name_branch

_RING_BONDS = [
    (1, 2), (2, 3), (3, 4), (4, 5), (5, 10), (10, 1), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
    (9, 11), (11, 12), (12, 13), (13, 14), (14, 8), (14, 15), (15, 16), (16, 17), (17, 13),
]
_ATTACHED = {
    "gon": [],
    "estr": [(13, 18)],
    "androst": [(13, 18), (10, 19)],
    "pregn": [(13, 18), (10, 19), (17, 20), (20, 21)],
    "chol": [(13, 18), (10, 19), (17, 20), (20, 21), (20, 22), (22, 23), (23, 24)],
    "cholest": [(13, 18), (10, 19), (17, 20), (20, 21), (20, 22), (22, 23), (23, 24), (24, 25), (25, 26), (25, 27)],
}
_ORDER = ["cholest", "chol", "pregn", "androst", "estr", "gon"]
_RING_POSITIONS = {p for bond in _RING_BONDS for p in bond}
_NATURAL_FACE = {8: "b", 9: "a", 10: "b", 13: "b", 14: "a"}
_SUFFIX = {"acid": "oic acid", "ketone": "one", "alcohol": "ol", "amine": "amine"}
_PREFIX = {"alcohol": "hydroxy", "ketone": "oxo", "amine": "amino"}
_SENIORITY = ["acid", "ketone", "alcohol", "amine"]
_FACES = {"a": "α", "b": "β", "x": "ξ"}


class Loc(int):
    """A steroid locant that may carry a configuration symbol: 17β."""

    def __new__(cls, number, face=""):
        obj = super().__new__(cls, number)
        obj.face = face
        return obj

    def __str__(self):
        return f"{int(self)}{_FACES.get(self.face, '')}"

    __repr__ = __str__
    __format__ = lambda self, spec: str(self)


def _template(stem):
    bonds = list(_RING_BONDS) + _ATTACHED[stem]
    adjacency_of = {}
    for a, b in bonds:
        adjacency_of.setdefault(a, set()).add(b)
        adjacency_of.setdefault(b, set()).add(a)
    return adjacency_of


def _embeddings(stem, graph, ring_atoms):
    """Every assignment {steroid position: atom} of the parent skeleton: ring positions on ring atoms, chain ones on
    acyclic carbons."""
    wanted = _template(stem)
    order, seen, queue = [], set(), [1]
    while queue:
        position = queue.pop(0)
        if position in seen:
            continue
        seen.add(position)
        order.append(position)
        queue.extend(sorted(wanted[position] - seen))
    results = []

    def extend(assigned, used):
        if len(assigned) == len(order):
            results.append(dict(assigned))
            return
        position = order[len(assigned)]
        anchors = [assigned[p] for p in wanted[position] if p in assigned]
        pool = graph[anchors[0]] if anchors else set(graph)
        for atom in sorted(pool):
            if atom in used or any(atom not in graph[x] for x in anchors):
                continue
            if (position in _RING_POSITIONS) != (atom in ring_atoms):
                continue
            assigned[position] = atom
            used.add(atom)
            extend(assigned, used)
            del assigned[position]
            used.discard(atom)

    extend({}, set())
    return results


def _core_atoms(mol, graph):
    ring_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring_atoms) or len(ring_atoms) != 17:
        return None
    info = mol.GetRingInfo()
    if info.NumRings() != 4:
        return None
    return ring_atoms


def _ene_locant(pair):
    low, high = pair
    return f"{low}({high})" if high - low != 1 else f"{low}"


def _plane_faces(mol, mapping):
    """{atom: 'a'|'b'} for every atom attached to a ring atom, from a 3D embedding; None when it cannot be embedded."""
    hydrogens = Chem.AddHs(mol)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20240607
    if AllChem.EmbedMolecule(hydrogens, parameters) != 0:
        parameters.useRandomCoords = True
        if AllChem.EmbedMolecule(hydrogens, parameters) != 0:
            return None
    coordinates = hydrogens.GetConformer().GetPositions()
    ring = np.array([coordinates[mapping[p]] for p in sorted(_RING_POSITIONS)])
    centroid = ring.mean(axis=0)
    normal = np.linalg.svd(ring - centroid)[2][2]
    turn = np.zeros(3)
    cycle = [1, 2, 3, 4, 5, 10]
    for i in range(6):
        a, b, c = (coordinates[mapping[cycle[(i + k) % 6]]] for k in range(3))
        turn += np.cross(b - a, c - b)
    if np.dot(turn, normal) < 0:
        normal = -normal
    return hydrogens, coordinates, normal


def _exocyclic_neighbors(hydrogens, atom, ring_atoms):
    return [n.GetIdx() for n in hydrogens.GetAtomWithIdx(atom).GetNeighbors() if n.GetIdx() not in ring_atoms]


def _face_of(coordinates, normal, center, others, atom):
    """'b' when `atom` lies on the β side of `center` relative to the other exocyclic neighbours."""
    height = {x: float(np.dot(coordinates[x] - coordinates[center], normal)) for x in [atom] + others}
    if others:
        return "b" if height[atom] > max(height[x] for x in others) else "a"
    return "b" if height[atom] > 0 else "a"


def _chirality_signature(hydrogens, atom, key):
    tag = hydrogens.GetAtomWithIdx(atom).GetChiralTag()
    if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
        return None
    neighbors = [n.GetIdx() for n in hydrogens.GetAtomWithIdx(atom).GetNeighbors()]
    ranked = sorted(neighbors, key=key)
    swaps = 0
    permutation = [neighbors.index(x) for x in ranked]
    for i in range(len(permutation)):
        for j in range(i + 1, len(permutation)):
            if permutation[i] > permutation[j]:
                swaps += 1
    base = 1 if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW else 2
    return base if swaps % 2 == 0 else 3 - base


_NATURAL_C20 = {}


def _natural_c20_signature():
    if "value" not in _NATURAL_C20:
        mol = Chem.MolFromSmiles("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        graph = {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 6} for a in mol.GetAtoms()}
        ring_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
        mapping = _embeddings("cholest", graph, ring_atoms)[0]
        hydrogens = Chem.AddHs(mol)
        position_of = {atom: p for p, atom in mapping.items()}
        key = lambda n: position_of.get(n, 0)
        _NATURAL_C20["value"] = _chirality_signature(hydrogens, mapping[20], key)
    return _NATURAL_C20["value"]


def name_steroid(mol):
    """Name of a substituted steroid on a retained parent, or None when `mol` is not one this module handles."""
    if mol.GetRingInfo().NumRings() != 4 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    if any(a.GetAtomicNum() not in (1, 6, 7, 8, 9, 16, 17, 35, 53) for a in mol.GetAtoms()):
        return None
    graph = {
        a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 6}
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
    }
    ring_atoms = _core_atoms(mol, graph)
    if ring_atoms is None:
        return None
    for stem in _ORDER:
        found = _embeddings(stem, graph, ring_atoms)
        if found:
            try:
                return _assemble(mol, stem, found, ring_atoms)
            except UnsupportedStructure:
                return None
    return None


def _groups(mol, mapping):
    """({class: [(position, atoms)]}, [(position, branch root)]) for the groups on the skeleton atoms."""
    mapped = set(mapping.values())
    classes, branches = {}, []
    for position, atom in mapping.items():
        center = mol.GetAtomWithIdx(atom)
        outside = [n for n in center.GetNeighbors() if n.GetIdx() not in mapped]
        oxo = [
            n.GetIdx()
            for n in outside
            if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl = [n.GetIdx() for n in outside if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1]
        amino = [n.GetIdx() for n in outside if n.GetAtomicNum() == 7 and n.GetDegree() == 1 and n.GetTotalNumHs() == 2]
        if any(mol.GetBondBetweenAtoms(atom, n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in outside if n.GetIdx() not in oxo):
            raise UnsupportedStructure("an exocyclic multiple bond on a steroid is not supported")
        if oxo and hydroxyl and len(outside) == 2 and position == 24:
            classes.setdefault("acid", []).append((position, {atom, oxo[0], hydroxyl[0]}))
            continue
        if oxo:
            if len(outside) != 1 or (position not in _RING_POSITIONS and center.GetTotalNumHs() > 0):
                raise UnsupportedStructure("an aldehyde or acyl group on a steroid chain is not supported")
            classes.setdefault("ketone", []).append((position, {oxo[0]}))
            continue
        for neighbor in outside:
            idx = neighbor.GetIdx()
            if idx in hydroxyl:
                classes.setdefault("alcohol", []).append((position, {idx}))
            elif idx in amino:
                classes.setdefault("amine", []).append((position, {idx}))
            else:
                if neighbor.GetAtomicNum() == 6 and position >= 22:
                    raise UnsupportedStructure("a carbon substituent on the side chain needs a larger steroid parent")
                branches.append((position, idx))
    return classes, branches


def _words(base_stem, ene, yne, principal_word, principal_locants):
    segments = []
    if ene:
        segments.append((ene, multiplied_word(len(ene), "ene")))
    if yne:
        segments.append((yne, multiplied_word(len(yne), "yne")))
    if not segments:
        segments.append(([], "ane"))
    if principal_word:
        segments.append((principal_locants, principal_word))
    text = base_stem + ("a" if segments[0][1].startswith(("di", "tri", "tetra")) else "")
    for i, (locants, word) in enumerate(segments):
        if word.endswith("e") and i + 1 < len(segments) and segments[i + 1][1][0] in "aeiouy":
            word = word[:-1]
        text += (f"-{','.join(locants)}-" if locants else "") + word
    return text


def _assemble(mol, stem, mappings, ring_atoms):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    best = None
    for mapping in mappings:
        classes, branches = _groups(mol, mapping)
        principal = next((c for c in _SENIORITY if c in classes), None)
        key = (
            tuple(sorted(p for p, _ in classes.get(principal, []))),
            tuple(mapping[p] for p in sorted(mapping)),
        )
        if best is None or key < best[0]:
            best = (key, mapping, classes, branches, principal)
    _, mapping, classes, branches, principal = best
    position_of = {atom: p for p, atom in mapping.items()}
    double, triple = _bonds(mol, mapping, position_of)
    if not classes and not branches and not double and not triple:
        raise UnsupportedStructure("a bare steroid parent hydride is named by its retained-name module")
    hydrocarbon = all(a.GetAtomicNum() in (1, 6) for a in mol.GetAtoms())
    if hydrocarbon and branches:
        raise UnsupportedStructure("a methylated steroid hydrocarbon may need a nor/homo/seco parent")
    if stem == "gon" and any(p in (10, 13) for p, _ in branches):
        raise UnsupportedStructure("a carbon on C-10 or C-13 of gonane is a nor-steroid")
    faces, parent_descriptors = _configuration(mol, stem, mapping, classes, branches)

    prefix_groups = {}
    for cls, members in classes.items():
        if cls == principal:
            continue
        for position, atoms in members:
            prefix_groups.setdefault(_PREFIX[cls], {"locants": [], "compound": False})["locants"].append(
                Loc(position, faces.get(next(iter(atoms)), ""))
            )
    for position, root in branches:
        name, compound = name_branch(graph, root, mapping[position], halogens, frozenset(), mol=mol, unsaturated=True)
        prefix_groups.setdefault(name, {"locants": [], "compound": compound})["locants"].append(
            Loc(position, faces.get(root, ""))
        )
    prefix = format_substituent_prefixes(prefix_groups) if prefix_groups else ""

    if principal is None:
        word, locants = "", []
    else:
        members = sorted(classes[principal], key=lambda m: m[0])
        word = multiplied_word(len(members), _SUFFIX[principal])
        locants = [str(Loc(p, faces.get(next(iter(sorted(atoms))), ""))) for p, atoms in members]
    core = _words(stem, [_ene_locant(b) for b in double], [_ene_locant(b) for b in triple], word, locants)
    descriptor = ",".join(parent_descriptors) + "-" if parent_descriptors else ""
    if prefix and descriptor:
        return f"{prefix}-{descriptor}{core}"
    return f"{prefix}{core}" if prefix else f"{descriptor}{core}"


def _bonds(mol, mapping, position_of):
    ring_a = [mapping[p] for p in (1, 2, 3, 4, 5, 10)]
    aromatic_a = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_a)
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for p, a in mapping.items() if p not in (1, 2, 3, 4, 5, 10)):
        raise UnsupportedStructure("an aromatic ring other than ring A of a steroid is not supported")
    double, triple = set(), set()
    if aromatic_a:
        double |= {(1, 2), (3, 4), (5, 10)}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()
        if bond.GetIsAromatic() or order == 1.0:
            continue
        if a not in position_of and b not in position_of:
            continue
        if a not in position_of or b not in position_of:
            if order == 2.0 and (mol.GetAtomWithIdx(a).GetAtomicNum() == 8 or mol.GetAtomWithIdx(b).GetAtomicNum() == 8):
                continue
            raise UnsupportedStructure("an exocyclic multiple bond on a steroid is not supported")
        pair = tuple(sorted((position_of[a], position_of[b])))
        (double if order == 2.0 else triple).add(pair)
    return sorted(double), sorted(triple)


def _configuration(mol, stem, mapping, classes, branches):
    """({exocyclic atom: face letter}, parent descriptors like '5α') from the specified centres."""
    potential = {e.centeredOn: e for e in Chem.FindPotentialStereo(mol) if e.type == Chem.StereoType.Atom_Tetrahedral}
    specified = {a for a, e in potential.items() if e.specified == Chem.StereoSpecified.Specified}
    if not specified & set(mapping.values()):
        return {}, []
    embedded = _plane_faces(mol, mapping)
    if embedded is None:
        raise UnsupportedStructure("the steroid could not be embedded in three dimensions")
    hydrogens, coordinates, normal = embedded
    ring_atoms = {mapping[p] for p in _RING_POSITIONS}
    faces, descriptors = {}, []
    for position, atom in sorted(mapping.items(), key=lambda item: item[0]):
        if atom not in potential:
            continue
        if position not in _RING_POSITIONS:
            if position != 20 or stem not in ("chol", "cholest"):
                if atom in specified:
                    raise UnsupportedStructure("a stereocentre on a steroid chain beyond C-20 is not supported")
                continue
            if atom in specified:
                if _c20_signature(hydrogens, mapping) != _natural_c20_signature():
                    raise UnsupportedStructure("an unnatural C-20 configuration is not supported")
            else:
                descriptors.append((20, "20ξ"))
            continue
        exocyclic = _exocyclic_neighbors(hydrogens, atom, ring_atoms)
        if atom in specified:
            face = {x: _face_of(coordinates, normal, atom, [y for y in exocyclic if y != x], x) for x in exocyclic}
        else:
            face = {x: "x" for x in exocyclic}
        heavy = [x for x in exocyclic if hydrogens.GetAtomWithIdx(x).GetAtomicNum() != 1]
        single = exocyclic[0] if len(exocyclic) == 1 else None
        if position in _NATURAL_FACE:
            if face[single] != _NATURAL_FACE[position]:
                descriptors.append((position, f"{position}{_FACES[face[single]]}"))
            continue
        if position == 5:
            descriptors.append((5, f"5{_FACES[face[single]]}"))
            continue
        if position == 17 and stem in ("pregn", "chol", "cholest"):
            chain = mapping[20]
            if face[chain] != "b":
                descriptors.append((17, f"17{_FACES[face[chain]]}"))
            continue
        for x in heavy:
            faces[x] = face[x]
    return faces, [text for _, text in sorted(descriptors)]


def _c20_signature(hydrogens, mapping):
    position_of = {atom: p for p, atom in mapping.items()}
    return _chirality_signature(hydrogens, mapping[20], lambda n: position_of.get(n, 0))
