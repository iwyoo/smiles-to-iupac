"""Inositol derivatives (P-104.3.1): O-substituted, esterified and amino-deoxy inositols named on the retained parent.

Numbering keeps the parent's face pattern (P-104.2.1) and is chosen by P-104.2.3 (d)-(f); D/L follows the direction of
numbering with hydroxy 1 above the ring, and is omitted for achiral derivatives.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, multiplied_word, numerical_term

_CW, _CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW


def _substituted_ring(mol, ring):
    members = set(ring)
    for a in ring:
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge() or atom.GetIsotope():
            return False
        if atom.GetTotalNumHs() != 1:
            return False
        outside = [n for n in atom.GetNeighbors() if n.GetIdx() not in members]
        if len(outside) != 1 or outside[0].GetAtomicNum() not in (7, 8):
            return False
    return True


def _ring_cycle(mol):
    """The single six-membered carbon ring whose atoms each carry one hydrogen and one O or N substituent, in ring order."""
    found = [ring for ring in mol.GetRingInfo().AtomRings() if len(ring) == 6 and _substituted_ring(mol, ring)]
    return found[0] if len(found) == 1 else None


def _parity(order, reference):
    """Parity of the permutation taking `reference` to `order` (0 even, 1 odd)."""
    positions = [reference.index(x) for x in order]
    inversions = sum(1 for i in range(4) for j in range(i + 1, 4) if positions[i] > positions[j])
    return inversions % 2


def _face(explicit, prev, atom, nxt, sub):
    """+1 or -1: the side of the ring, seen along prev -> atom -> nxt, on which `sub` lies."""
    center = explicit.GetAtomWithIdx(atom)
    tag = center.GetChiralTag()
    if tag not in (_CW, _CCW):
        return None
    neighbors = [b.GetOtherAtomIdx(atom) for b in center.GetBonds()]
    hydrogen = next(n for n in neighbors if n not in (prev, nxt, sub))
    odd = _parity([prev, nxt, sub, hydrogen], neighbors)
    return 1 if (tag == _CCW) != bool(odd) else -1


def _faces(explicit, ring, substituents, start, step):
    order = [ring[(start + step * i) % 6] for i in range(6)]
    faces = []
    for i, atom in enumerate(order):
        face = _face(explicit, order[i - 1], atom, order[(i + 1) % 6], substituents[atom])
        if face is None:
            return None, None
        faces.append(face)
    return order, faces


_PATTERNS = {
    "cis-inositol": {1, 2, 3, 4, 5, 6},
    "epi-inositol": {1, 2, 3, 4, 5},
    "allo-inositol": {1, 2, 3, 4},
    "myo-inositol": {1, 2, 3, 5},
    "muco-inositol": {1, 2, 4, 5},
    "neo-inositol": {1, 2, 3},
    "scyllo-inositol": {1, 3, 5},
    "chiro-inositol": {1, 2, 4},
}


def _l_face():
    """Face value of hydroxy 1 for clockwise numbering, calibrated on 1L-1,2/3,5-cyclohexanetetrol = (1R,2R,3R,5R) (P-104.2.3)."""
    from rdkit.Chem import rdCIPLabeler
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

    ring_locants = [1, 2, 4, 6, 7, 9]
    for isomer in EnumerateStereoisomers(Chem.MolFromSmiles("OC1C(O)C(O)CC(O)C1")):
        rdCIPLabeler.AssignCIPLabels(isomer)
        if all(isomer.GetAtomWithIdx(i).GetProp("_CIPCode") == "R" for i in (1, 2, 4, 7)):
            explicit = Chem.AddHs(isomer)
            return _face(explicit, ring_locants[-1], 1, ring_locants[1], 0)
    raise RuntimeError("calibration isomer not found")


_L_FACE = _l_face()


def ester_anion(mol, oxygen, acyl):
    """Anion name of the acid whose ester oxygen is `oxygen` and acyl-side atom is `acyl`; None when unsupported."""
    from ._acid_derivatives import anion_name
    from .core import smiles_to_iupac

    atom = mol.GetAtomWithIdx(acyl)
    if atom.GetAtomicNum() == 15:
        neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() != oxygen]
        if atom.GetDegree() == 4 and all(n.GetAtomicNum() == 8 and n.GetDegree() == 1 and not n.GetFormalCharge() for n in neighbors):
            return "dihydrogen phosphate"
        return None
    if atom.GetAtomicNum() == 16:
        neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() != oxygen]
        if atom.GetDegree() == 4 and all(n.GetAtomicNum() == 8 and n.GetDegree() == 1 and not n.GetFormalCharge() for n in neighbors):
            return "hydrogen sulfate"
        return None
    graph = adjacency(mol)
    branch, stack = {oxygen}, [acyl]
    while stack:
        current = stack.pop()
        if current in branch:
            continue
        branch.add(current)
        stack.extend(n for n in graph[current] if n not in branch)
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - branch, reverse=True):
        editable.RemoveAtom(index)
    for atom in editable.GetAtoms():
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            return None
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    try:
        return anion_name(smiles_to_iupac(Chem.MolToSmiles(acid)))
    except (UnsupportedStructure, ValueError):
        return None


def _classify(mol, ring):
    """Per ring carbon: (substituent atom, kind, payload) with kind hydroxy/amino/ether/ester; None if unsupported."""
    from ._substituents import name_branch

    graph = adjacency(mol)
    members = set(ring)
    result, covered = {}, set(ring)
    for a in ring:
        hetero = next(n for n in graph[a] if n not in members)
        atom = mol.GetAtomWithIdx(hetero)
        covered.add(hetero)
        if atom.GetAtomicNum() == 7:
            if atom.GetDegree() != 1 or atom.GetFormalCharge() or atom.GetTotalNumHs() != 2:
                return None
            result[a] = (hetero, "amino", None)
            continue
        if atom.GetFormalCharge() or atom.GetIsotope():
            return None
        if atom.GetDegree() == 1:
            result[a] = (hetero, "hydroxy", None)
            continue
        other = next(n for n in graph[hetero] if n != a)
        beyond = mol.GetAtomWithIdx(other)
        if beyond.GetAtomicNum() in (15, 16) or any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(beyond).GetAtomicNum() == 8 for b in beyond.GetBonds()
        ):
            anion = ester_anion(mol, hetero, other)
            if anion is None:
                return None
            result[a] = (hetero, "ester", anion)
        else:
            if beyond.GetAtomicNum() != 6:
                return None
            try:
                result[a] = (hetero, "ether", name_branch(graph, other, hetero, {}, frozenset(), mol=mol))
            except UnsupportedStructure:
                return None
    return result


def _other_stereo(mol, ring):
    members = set(ring)
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        return True
    return any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and a.GetIdx() not in members for a in mol.GetAtoms())


def _numbering_key(order, kinds):
    modified = {}
    for locant, atom in enumerate(order, 1):
        kind, payload = kinds[atom]
        if kind == "hydroxy":
            continue
        label = (kind, payload[0] if kind == "ether" else payload)
        modified.setdefault(label, []).append(locant)
    flat = sorted(l for locs in modified.values() for l in locs)
    by_name = tuple(
        tuple(modified[label]) for label in sorted(modified, key=lambda lab: (alpha_sort_key(_label_text(lab)), lab[0]))
    )
    return tuple(flat), by_name


def _label_text(label):
    kind, payload = label
    return {"amino": "amino", "ether": payload, "ester": payload}[kind]


def has_inositol_derivative_shape(mol) -> bool:
    return _inositol_parts(mol) is not None


def _inositol_parts(mol):
    ring = _ring_cycle(mol)
    if ring is None or _other_stereo(mol, ring):
        return None
    classified = _classify(mol, ring)
    if classified is None or all(kind == "hydroxy" for _, kind, _ in classified.values()):
        return None
    explicit = Chem.AddHs(mol)
    substituents = {a: info[0] for a, info in classified.items()}
    kinds = {a: (info[1], info[2]) for a, info in classified.items()}
    options = []
    for start in range(6):
        for step in (1, -1):
            order, faces = _faces(explicit, ring, substituents, start, step)
            if order is None:
                return None
            up = frozenset(i + 1 for i, f in enumerate(faces) if f == faces[0])
            options.append((order, up, "L" if faces[0] == _L_FACE else "D"))
    parents = {name for name, pattern in _PATTERNS.items() for _, up, _ in options if up == pattern}
    if len(parents) != 1:
        return None
    parent = parents.pop()
    valid = [(order, descriptor) for order, up, descriptor in options if up == _PATTERNS[parent]]
    keyed = [(_numbering_key(order, kinds), descriptor, order) for order, descriptor in valid]
    best = min(key for key, _, _ in keyed)
    winners = [(descriptor, order) for key, descriptor, order in keyed if key == best]
    descriptors = {descriptor for descriptor, _ in winners}
    chosen = "L" if "L" in descriptors else "D"
    order = next(order for descriptor, order in winners if descriptor == chosen)
    return parent, (None if len(descriptors) == 2 else chosen), order, kinds


_COMPLEX = {2: "bis", 3: "tris", 4: "tetrakis", 5: "pentakis", 6: "hexakis"}


def _multiplier(count, compound):
    if count == 1:
        return ""
    return (_COMPLEX[count] if compound else numerical_term(count)) + "-"


def _prefix_text(label, locants):
    kind, payload = label
    if kind == "amino":
        return [("amino", f"{','.join(map(str, locants))}-{multiplied_word(len(locants), 'amino')}"),
                ("deoxy", f"{','.join(map(str, locants))}-{multiplied_word(len(locants), 'deoxy')}")]
    name, compound = payload
    word = f"({name})" if compound else name
    multiplier = _multiplier(len(locants), compound)
    return [(name, f"{','.join(map(str, locants))}-{multiplier}O-{word}")]


def name_inositol_derivative(mol) -> str:
    parent, descriptor, order, kinds = _inositol_parts(mol)
    groups = {}
    for locant, atom in enumerate(order, 1):
        kind, payload = kinds[atom]
        if kind == "hydroxy":
            continue
        label = (kind, payload[0] if kind == "ether" else payload)
        groups.setdefault(label, []).append(locant)
    prefixes, esters = [], []
    for label, locants in groups.items():
        if label[0] == "ester":
            esters.append((label[1], locants))
        elif label[0] == "ether":
            kind_payload = next(kinds[a][1] for a in order if kinds[a][0] == "ether" and kinds[a][1][0] == label[1])
            prefixes.extend(_prefix_text(("ether", kind_payload), locants))
        else:
            prefixes.extend(_prefix_text(label, locants))
    text = "-".join(piece for _, piece in sorted(prefixes, key=lambda item: alpha_sort_key(item[0])))
    stem = f"{parent}-inositol" if not parent.endswith("inositol") else parent
    head = (f"1{descriptor}-" if descriptor else "") + (f"{text}-" if text else "") + stem
    if not esters:
        return head
    return f"{head} {ester_words(esters)}"


def ester_words(esters):
    """'2,4-diacetate 6-(dihydrogen phosphate)' from [(anion name, [locants])], alphabetized by anion."""
    words = []
    for anion, locants in sorted(esters, key=lambda item: alpha_sort_key(item[0])):
        count = len(locants)
        joined = ",".join(map(str, locants))
        compound = " " in anion
        if count == 1:
            words.append(f"{joined}-({anion})" if compound else f"{joined}-{anion}")
        elif compound:
            words.append(f"{joined}-{_COMPLEX[count]}({anion})")
        else:
            words.append(f"{joined}-{numerical_term(count)}{anion}")
    return " ".join(words)
