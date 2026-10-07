"""Substituted steroids named on the retained steroid parent hydrides (P-101.2, P-101.4): the 1989 steroid numbering
fixes every locant, unsaturation is cited as 'ene' endings (estra-1,3,5(10)-triene), groups as suffixes or prefixes and
configuration with α/β/ξ locants (cholest-5-en-3β-ol, 17β-hydroxy-5α-androstan-3-one). Configuration is read from a 3D
embedding: β is the face from which the numbering of ring A runs anticlockwise."""

from collections import Counter

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, halogen_substituents, multiplied_word
from ._numerals import multiplying_prefix
from ._substituents import cited_branch_stereo, format_substituent_prefixes, name_branch

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
_SUFFIX = {"ketone": "one", "alcohol": "ol", "amine": "amine", "ester_o": "yl"}
_ACYL_SUFFIX = {
    ("acid", "c"): "carboxylic acid",
    ("acid", "o"): "oic acid",
    ("ester", "c"): "carboxylate",
    ("ester", "o"): "oate",
    ("amide", "c"): "carboxamide",
    ("amide", "o"): "amide",
}
_PREFIX = {"alcohol": "hydroxy", "ketone": "oxo", "amine": "amino"}
_SENIORITY = ["acid", "ester", "ester_o", "amide", "ketone", "alcohol", "amine"]
_ACYL_CLASSES = ("acid", "ester", "ester_o", "amide")
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


def _carboxyl(mol, carbon, parent, mapped):
    """('acid'|'ester'|'amide', ester alkyl carbon or None) when `carbon` is a C(=O)X group attached to `parent`."""
    others = [n for n in mol.GetAtomWithIdx(carbon).GetNeighbors() if n.GetIdx() != parent]
    oxo = [n for n in others if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    rest = [n for n in others if n not in oxo]
    if len(oxo) != 1 or len(rest) != 1 or mol.GetBondBetweenAtoms(carbon, rest[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    x = rest[0]
    if x.GetAtomicNum() == 8 and x.GetDegree() == 1 and x.GetTotalNumHs() == 1:
        return "acid", None
    if x.GetAtomicNum() == 8 and x.GetDegree() == 2:
        alkyl = next(n for n in x.GetNeighbors() if n.GetIdx() != carbon)
        if alkyl.GetAtomicNum() != 6 or alkyl.GetIdx() in mapped:
            raise UnsupportedStructure("a lactone or a non-carbon ester group on a steroid is not supported")
        return "ester", (x.GetIdx(), alkyl.GetIdx())
    if x.GetAtomicNum() == 7 and x.GetDegree() == 1 and x.GetTotalNumHs() == 2:
        return "amide", None
    return None


def _acyl_oxygen(mol, oxygen, parent, mapped):
    """The acyl carbon of an O-acyl group (R-C(=O)-O-) on `parent`, else None."""
    acyl = next((n for n in mol.GetAtomWithIdx(oxygen).GetNeighbors() if n.GetIdx() != parent), None)
    if acyl is None or acyl.GetAtomicNum() != 6 or acyl.GetIdx() in mapped:
        return None
    oxo = [
        n for n in acyl.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(acyl.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if len(oxo) != 1 or acyl.GetDegree() > 3:
        return None
    if any(n.GetIdx() not in (oxygen, oxo[0].GetIdx()) and n.GetAtomicNum() != 6 for n in acyl.GetNeighbors()):
        return None
    return acyl.GetIdx()


def _groups(mol, mapping, terminals, parent_of):
    """({class: [(position, atoms, extra)]}, [(position, branch root)]) for the groups on the skeleton atoms."""
    mapped = set(mapping.values())
    classes, branches = {}, []
    for position, atom in mapping.items():
        center = mol.GetAtomWithIdx(atom)
        outside = [n for n in center.GetNeighbors() if n.GetIdx() not in mapped]
        if position in terminals:
            carboxyl = _carboxyl(mol, atom, mapping[parent_of[position]], mapped)
            if carboxyl:
                kind, alkyl = carboxyl
                extra = {"kind": "o", "ester": alkyl}
                classes.setdefault(kind, []).append((position, {atom}, extra))
                continue
        oxo = [
            n.GetIdx()
            for n in outside
            if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl = [n.GetIdx() for n in outside if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1]
        amino = [n.GetIdx() for n in outside if n.GetAtomicNum() == 7 and n.GetDegree() == 1 and n.GetTotalNumHs() == 2]
        if any(mol.GetBondBetweenAtoms(atom, n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in outside if n.GetIdx() not in oxo):
            raise UnsupportedStructure("an exocyclic multiple bond on a steroid is not supported")
        if oxo:
            if len(outside) != 1 or (position not in _RING_POSITIONS and center.GetTotalNumHs() > 0):
                raise UnsupportedStructure("an aldehyde or acyl group on a steroid chain is not supported")
            classes.setdefault("ketone", []).append((position, {oxo[0]}, {}))
            continue
        for neighbor in outside:
            idx = neighbor.GetIdx()
            if idx in hydroxyl:
                classes.setdefault("alcohol", []).append((position, {idx}, {}))
            elif idx in amino:
                classes.setdefault("amine", []).append((position, {idx}, {}))
            else:
                if neighbor.GetAtomicNum() == 6 and position >= 22:
                    raise UnsupportedStructure("a carbon substituent on the side chain needs a larger steroid parent")
                carboxyl = _carboxyl(mol, idx, atom, mapped) if neighbor.GetAtomicNum() == 6 else None
                if carboxyl:
                    kind, alkyl = carboxyl
                    extra = {"kind": "c", "ester": alkyl, "anchor": idx, "root": idx}
                    classes.setdefault(kind, []).append((position, {idx}, extra))
                    continue
                acyl = _acyl_oxygen(mol, idx, atom, mapped) if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 2 else None
                if acyl is not None:
                    extra = {"anchor": idx, "root": idx, "acyl": acyl}
                    classes.setdefault("ester_o", []).append((position, {idx}, extra))
                    continue
                branches.append((position, idx))
    return classes, branches


def _acid_anion(mol, oxygen, acyl, mapped):
    """The anion name of the acid part of an O-acyl group, e.g. 'ethanoate'."""
    from .core import smiles_to_iupac

    keep, stack = {acyl}, [acyl]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() != oxygen and n.GetIdx() not in keep:
                if n.GetIdx() in mapped:
                    raise UnsupportedStructure("a ring-closing O-acyl group on a steroid is not supported")
                keep.add(n.GetIdx())
                stack.append(n.GetIdx())
    editable = Chem.RWMol(mol)
    ester_oxygen = editable.GetAtomWithIdx(oxygen)
    ester_oxygen.SetNumExplicitHs(1)
    ester_oxygen.SetNoImplicit(True)
    for idx in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in keep | {oxygen}), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    name = smiles_to_iupac(Chem.MolToSmiles(acid))
    if not name.endswith("ic acid") or " " in name[: -len(" acid")]:
        raise UnsupportedStructure("the acid part of this steroid ester is not a plain acid")
    return name[: -len("ic acid")] + "ate"


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


def _anchor_face(faces, atoms, extra):
    if extra.get("kind") == "o":
        return ""
    return faces.get(extra.get("anchor", next(iter(sorted(atoms)))), "")


def _assemble(mol, stem, mappings, ring_atoms):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    template = _template(stem)
    terminals = {p for p in template if p not in _RING_POSITIONS and len(template[p]) == 1}
    parent_of = {p: next(iter(template[p])) for p in terminals}
    best = None
    for mapping in mappings:
        classes, branches = _groups(mol, mapping, terminals, parent_of)
        principal = next((c for c in _SENIORITY if c in classes), None)
        key = (
            tuple(sorted(m[0] for m in classes.get(principal, []))),
            tuple(mapping[p] for p in sorted(mapping)),
        )
        if best is None or key < best[0]:
            best = (key, mapping, classes, branches, principal)
    _, mapping, classes, branches, principal = best
    if "ester" in classes and "ester_o" in classes:
        raise UnsupportedStructure("esters of both a steroid acid and a steroid alcohol are not supported")
    for cls in _ACYL_CLASSES:
        if cls in classes and cls != principal:
            for position, _, extra in classes.pop(cls):
                if "root" not in extra:
                    raise UnsupportedStructure("a terminal acyl group below the principal group is not supported")
                branches.append((position, extra["root"]))
    position_of = {atom: p for p, atom in mapping.items()}
    double, triple = _bonds(mol, mapping, position_of)
    if not classes and not branches and not double and not triple:
        raise UnsupportedStructure("a bare steroid parent hydride is named by its retained-name module")
    hydrocarbon = all(a.GetAtomicNum() in (1, 6) for a in mol.GetAtoms())
    if hydrocarbon and any(p not in _RING_POSITIONS for p, _ in branches):
        raise UnsupportedStructure("an alkylated steroid chain may need a nor/homo/seco parent")
    if stem == "gon" and any(p in (10, 13) for p, _ in branches):
        raise UnsupportedStructure("a carbon on C-10 or C-13 of gonane is a nor-steroid")
    faces, parent_descriptors, side = _configuration(mol, stem, mapping)
    chain_descriptors = [text for _, text in sorted(side + _double_bond_descriptors(mol, position_of))]

    prefix_groups = {}
    for cls, members in classes.items():
        if cls == principal:
            continue
        for position, atoms, extra in members:
            prefix_groups.setdefault(_PREFIX[cls], {"locants": [], "compound": False})["locants"].append(
                Loc(position, _anchor_face(faces, atoms, extra))
            )
    ester_roots = [extra["ester"][1] for _, _, extra in classes.get("ester", []) if principal == "ester"]
    with cited_branch_stereo(mol, graph, set(mapping.values()), [root for _, root in branches] + ester_roots):
        named = [
            (position, root, *name_branch(graph, root, mapping[position], halogens, frozenset(), mol=mol, unsaturated=True))
            for position, root in branches
        ]
        ester_alkyls = [
            name_branch(graph, extra["ester"][1], extra["ester"][0], halogens, frozenset(), mol=mol, unsaturated=True)
            for _, _, extra in classes.get("ester", [])
            if principal == "ester"
        ]
    repeats = Counter((position, name) for position, _, name, _ in named)
    for position, root, name, compound in named:
        face = "" if repeats[(position, name)] > 1 else faces.get(root, "")
        prefix_groups.setdefault(name, {"locants": [], "compound": compound})["locants"].append(Loc(position, face))
    prefix = format_substituent_prefixes(prefix_groups) if prefix_groups else ""

    word, locants, alkyl_word, anion = "", [], "", ""
    if principal is not None:
        members = sorted(classes[principal], key=lambda m: m[0])
        if principal in ("acid", "ester", "amide"):
            if len({extra["kind"] for _, _, extra in members}) != 1:
                raise UnsupportedStructure("a mix of ring and chain acyl groups on a steroid is not supported")
            suffix = _ACYL_SUFFIX[(principal, members[0][2]["kind"])]
        else:
            suffix = _SUFFIX[principal]
        word = multiplied_word(len(members), suffix)
        locants = [str(Loc(p, _anchor_face(faces, atoms, extra))) for p, atoms, extra in members]
        if principal == "ester":
            alkyls = {}
            for name, compound in ester_alkyls:
                alkyls[name] = compound
            if len(alkyls) != 1:
                raise UnsupportedStructure("esters with different alkyl groups on a steroid are not supported")
            ((name, compound),) = alkyls.items()
            alkyl_word = name
            if len(members) > 1:
                alkyl_word = multiplying_prefix(len(members), compound=compound) + (f"({name})" if compound else name)
        if principal == "ester_o":
            mapped = set(mapping.values())
            anions = {_acid_anion(mol, extra["anchor"], extra["acyl"], mapped) for _, _, extra in members}
            if len(anions) != 1:
                raise UnsupportedStructure("O-acyl groups from different acids on a steroid are not supported")
            (anion,) = anions
            if len(members) > 1:
                if not anion.isalpha():
                    raise UnsupportedStructure("a substituted acid part repeated on a steroid is not supported")
                anion = multiplying_prefix(len(members), compound=False) + anion
    core = _words(stem, [_ene_locant(b) for b in double], [_ene_locant(b) for b in triple], word, locants)
    descriptor = ",".join(parent_descriptors) + "-" if parent_descriptors else ""
    chain = f"({','.join(chain_descriptors)})-" if chain_descriptors else ""
    if prefix and descriptor:
        body = f"{chain}{prefix}-{descriptor}{core}"
    else:
        body = f"{chain}{prefix}{descriptor}{core}"
    return " ".join(part for part in (alkyl_word, body, anion) if part)


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


def _configuration(mol, stem, mapping):
    """({exocyclic atom: face letter}, parent descriptors like '5α', chain CIP descriptors like '20S')."""
    potential = {e.centeredOn: e for e in Chem.FindPotentialStereo(mol) if e.type == Chem.StereoType.Atom_Tetrahedral}
    specified = {a for a, e in potential.items() if e.specified == Chem.StereoSpecified.Specified}
    if not specified & set(mapping.values()):
        return {}, [], []
    embedded = _plane_faces(mol, mapping)
    if embedded is None:
        raise UnsupportedStructure("the steroid could not be embedded in three dimensions")
    hydrogens, coordinates, normal = embedded
    ring_atoms = {mapping[p] for p in _RING_POSITIONS}
    faces, descriptors, side = {}, [], []
    for position, atom in sorted(mapping.items(), key=lambda item: item[0]):
        if atom not in potential:
            continue
        if position not in _RING_POSITIONS:
            if position != 20 or stem not in ("chol", "cholest"):
                if atom in specified:
                    side.append((position, f"{position}{_cip_label(mol, atom)}"))
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
    return faces, [text for _, text in sorted(descriptors)], side


def _double_bond_descriptors(mol, position_of):
    """[(locant, '23E')] for the specified chain double bonds of the skeleton."""
    result = []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bond.GetBondTypeAsDouble() != 2.0 or bond.IsInRing() or a not in position_of or b not in position_of:
            continue
        if bond.GetStereo() in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
            continue
        rdCIPLabeler.AssignCIPLabels(mol)
        if not bond.HasProp("_CIPCode"):
            raise UnsupportedStructure("a steroid chain double bond has no E/Z label")
        low = min(position_of[a], position_of[b])
        result.append((low, f"{low}{bond.GetProp('_CIPCode')}"))
    return result


def _cip_label(mol, atom):
    rdCIPLabeler.AssignCIPLabels(mol)
    center = mol.GetAtomWithIdx(atom)
    if not center.HasProp("_CIPCode") or center.GetProp("_CIPCode") not in ("R", "S"):
        raise UnsupportedStructure("a steroid chain centre has no CIP R/S label")
    return center.GetProp("_CIPCode")


def _c20_signature(hydrogens, mapping):
    position_of = {atom: p for p, atom in mapping.items()}
    return _chirality_signature(hydrogens, mapping[20], lambda n: position_of.get(n, 0))
