"""Deoxy, amino, thio, halo, O- and C-substituted monosaccharides, open-chain or cyclic (P-102.5.3 to P-102.5.6.4).

Every skeleton position is classified (hydroxy, deoxy, amino, thio, halogen, ether, ester, C-substituent) and the
molecule is restored to the plain sugar by putting a hydroxy group back. The configurational prefixes come from the
chirality centres that remain, grouped in fours from C-2 (P-102.5.1.1.2, P-102.5.3.4); a deoxy position that removes a
centre therefore gives a systematic name such as 2-deoxy-D-erythro-pentose, while one that leaves every centre keeps the
retained stem. Esters follow the sugar name as separate words, all other substituents are detachable prefixes.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._carbohydrate import _D_2_KETOSE_PATTERNS, _D_ALDOSE_PATTERNS, _FLIP_CIP
from ._cited_group import cited_group, subtree
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, numerical_term
from ._inositol_derivative import ester_anion, ester_words

_HALOGENS = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALIDE_WORDS = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_COMPLEX = {2: "bis", 3: "tris", 4: "tetrakis", 5: "pentakis", 6: "hexakis"}
_NUMBER_STEM = {4: "tetr", 5: "pent", 6: "hex", 7: "hept", 8: "oct", 9: "non", 10: "dec"}
_GROUP = 4
_ELEMENTS = {6, 7, 8, 9, 15, 16, 17, 34, 35, 52, 53}
_CHALCOGEN_WORDS = {16: "thio", 34: "seleno", 52: "telluro"}
_PREFIX = {pattern: ("glycero" if stem == "glyceraldehyde" else stem[:-2]) for pattern, stem in _D_ALDOSE_PATTERNS.items()}
_SIX_DEOXY = {("D", "glucose"): "quinovose", ("L", "galactose"): "fucose", ("L", "mannose"): "rhamnose"}


class _Skeleton:
    def __init__(self, carbons, ketose, hetero=None, closing=None, size=None):
        self.size = size
        self.carbons = carbons
        self.ketose = ketose
        self.hetero = hetero
        self.closing = closing
        self.n = len(carbons)

    @property
    def ring(self):
        return self.hetero is not None

    @property
    def anomeric(self):
        return 2 if self.ketose else 1

    def position(self, atom):
        return self.carbons.index(atom) + 1


def _walk_ring(ring_set, graph, start, away):
    path, previous, current = [], away, start
    while True:
        path.append(current)
        following = [n for n in graph[current] if n in ring_set and n != previous]
        if len(path) == len(ring_set) - 1:
            return path
        previous, current = current, following[0]


def _carbon_neighbors(mol, atom, excluded):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(atom).GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() not in excluded]


def _ring_skeleton(mol, graph):
    rings = [
        r
        for r in mol.GetRingInfo().AtomRings()
        if len(r) in (5, 6)
        and sum(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in r) == 1
        and not any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)
    ]
    if len(rings) != 1:
        return None
    ring = rings[0]
    ring_set = set(ring)
    hetero = next(a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() != 6)
    if mol.GetAtomWithIdx(hetero).GetAtomicNum() not in (8, *_CHALCOGEN_WORDS) or mol.GetAtomWithIdx(hetero).GetDegree() != 2:
        return None
    ends = [n for n in graph[hetero]]

    def exo(c):
        return [n for n in graph[c] if n not in ring_set]

    anomerics = [c for c in ends if any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for n in exo(c))]
    if len(anomerics) != 1:
        return None
    anomeric = anomerics[0]
    other = next(c for c in ends if c != anomeric)
    chain = _walk_ring(ring_set, graph, anomeric, hetero)
    ketose = any(mol.GetAtomWithIdx(n).GetAtomicNum() == 6 for n in exo(anomeric))
    size = len(ring)
    if ketose:
        carbon_exo = [n for n in exo(anomeric) if mol.GetAtomWithIdx(n).GetAtomicNum() == 6]
        if len(carbon_exo) != 1:
            return None
        carbons = [carbon_exo[0]] + chain
        tail = [n for n in _carbon_neighbors(mol, other, ring_set)]
        if size == 6:
            if tail:
                return None
        else:
            if len(tail) > 1:
                return None
            carbons += tail
    else:
        tail = _carbon_neighbors(mol, other, ring_set)
        if size == 6:
            if len(tail) > 1:
                return None
            carbons = chain + tail
        else:
            if len(tail) != 1:
                return None
            sixth = _carbon_neighbors(mol, tail[0], ring_set | {other})
            if len(sixth) > 1:
                return None
            carbons = chain + tail + sixth
    closing = carbons.index(other) + 1
    return _Skeleton(carbons, ketose, hetero, closing, size)


def _chain_skeleton(mol, graph):
    carbonyls = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
        and any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetAtomicNum() in (8, *_CHALCOGEN_WORDS) and b.GetOtherAtom(a).GetDegree() == 1
            for b in a.GetBonds()
        )
        and not any(n.GetAtomicNum() in (7, 8, 16, 34, 52) and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0 for n in a.GetNeighbors())
    ]
    if len(carbonyls) != 1:
        return None
    carbonyl = carbonyls[0]
    around = _carbon_neighbors(mol, carbonyl, set())

    def walk(start, previous):
        path = [start]
        while True:
            following = _carbon_neighbors(mol, path[-1], {previous} | set(path))
            if len(following) > 1:
                return None
            if not following:
                return path
            previous = path[-1]
            path.append(following[0])

    if len(around) == 1 and mol.GetAtomWithIdx(carbonyl).GetTotalNumHs() == 1:
        rest = walk(around[0], carbonyl)
        return None if rest is None else _Skeleton([carbonyl] + rest, False)
    if len(around) == 2:
        walks = [walk(a, carbonyl) for a in around]
        if None in walks:
            return None
        short, long_ = sorted(walks, key=len)
        return _Skeleton(short + [carbonyl] + long_, True) if len(short) == 1 else None
    return None


def _classify(mol, graph, atom, skeleton_atoms, skeleton):
    """(kind, exocyclic atom, root) of the substituent pattern on skeleton carbon `atom`."""
    exo = [n for n in graph[atom] if n not in skeleton_atoms]
    exo = [n for n in exo if n != skeleton.hetero]
    if not exo:
        return ("H", None, None)
    if len(exo) == 2:
        oxygens = [n for n in exo if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1]
        carbons = [n for n in exo if mol.GetAtomWithIdx(n).GetAtomicNum() == 6]
        if len(oxygens) == 1 and len(carbons) == 1:
            return ("C", oxygens[0], carbons[0])
        raise UnsupportedStructure("two substituents on a sugar carbon")
    (x,) = exo
    a = mol.GetAtomWithIdx(x)
    z = a.GetAtomicNum()
    if z in _CHALCOGEN_WORDS and a.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom, x).GetBondTypeAsDouble() == 2.0:
        return ("thione", x, None)
    if z == 8:
        if mol.GetBondBetweenAtoms(atom, x).GetBondTypeAsDouble() == 2.0:
            return ("oxo", x, None)
        if a.GetDegree() == 1:
            return ("OH", x, None)
        other = next(n for n in graph[x] if n != atom)
        b = mol.GetAtomWithIdx(other)
        if subtree(graph, other, x) & skeleton_atoms:
            raise UnsupportedStructure("a cyclic acetal, carbonate or phosphate bridges two skeleton carbons")
        if b.GetAtomicNum() in (15, 16):
            return ("ester", x, other)
        if b.GetAtomicNum() == 6 and any(bd.GetBondTypeAsDouble() == 2.0 and bd.GetOtherAtom(b).GetAtomicNum() == 8 for bd in b.GetBonds()):
            return ("ester", x, other)
        if b.GetAtomicNum() == 6 and other in skeleton_atoms:
            raise UnsupportedStructure("an ether between two skeleton carbons is an anhydro bridge")
        if b.GetAtomicNum() == 6:
            return ("ether", x, other)
        raise UnsupportedStructure("unsupported oxygen substituent")
    if z == 7 and not a.GetFormalCharge():
        return ("N", x, None)
    if z in _CHALCOGEN_WORDS and a.GetDegree() == 1:
        return ("SH", x, None)
    if z in _CHALCOGEN_WORDS and a.GetDegree() == 2:
        other = next(n for n in graph[x] if n != atom)
        if mol.GetAtomWithIdx(other).GetAtomicNum() == 6 and other not in skeleton_atoms:
            return ("Sether", x, other)
    if z == 6 and mol.GetBondBetweenAtoms(atom, x).GetBondTypeAsDouble() == 1.0:
        return ("Csub", x, x)
    if z in _HALOGENS and a.GetDegree() == 1:
        return ("X", x, None)
    raise UnsupportedStructure("unsupported substituent on a sugar carbon")


def _restore(mol, graph, skeleton, decorations):
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    removed = set()
    for carbon, (kind, x, root) in decorations.items():
        if kind == "H":
            oxygen = editable.AddAtom(Chem.Atom(8))
            editable.GetAtomWithIdx(oxygen).SetIntProp("_orig", -1)
            editable.AddBond(carbon, oxygen, Chem.BondType.SINGLE)
        elif kind in ("N", "SH", "X", "thione", "Csub", "Sether"):
            editable.GetAtomWithIdx(x).SetAtomicNum(8)
            editable.GetAtomWithIdx(x).SetIsAromatic(False)
            for n in graph[x]:
                if n != carbon:
                    removed |= subtree(graph, n, x)
        elif kind in ("ether", "ester"):
            removed |= subtree(graph, root, x)
        elif kind == "C":
            removed |= subtree(graph, root, carbon)
    if skeleton.ring:
        editable.GetAtomWithIdx(skeleton.hetero).SetAtomicNum(8)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    restored = editable.GetMol()
    Chem.SanitizeMol(restored)
    return restored, {a.GetIntProp("_orig"): a.GetIdx() for a in restored.GetAtoms() if a.GetIntProp("_orig") >= 0}


def _labels(restored, skeleton, atom_of):
    rdCIPLabeler.AssignCIPLabels(restored)

    def label(position):
        atom = restored.GetAtomWithIdx(atom_of[skeleton.carbons[position - 1]])
        return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None

    first = 3 if skeleton.ketose else 2
    centres = [(p, label(p)) for p in range(first, skeleton.n)]
    if skeleton.ring and skeleton.closing < skeleton.n:
        flip = skeleton.closing - 1
        centres = [(p, _FLIP_CIP[lab] if p == flip and lab else lab) for p, lab in centres]
    return centres, label(skeleton.anomeric) if skeleton.ring else None


def _mirror(letters):
    return tuple(_FLIP_CIP[c] for c in letters)


def _parent(skeleton, centres, anomeric_label, six_deoxy_only):
    """(core name without substituents, dropped_six_deoxy) or None."""
    specified = [(p, lab) for p, lab in centres if lab]
    if not specified:
        return None
    full = len(specified) == len(centres)
    table = _D_2_KETOSE_PATTERNS if skeleton.ketose else _D_ALDOSE_PATTERNS
    ring_word = ("furanose" if skeleton.size == 5 else "pyranose") if skeleton.ring else None
    groups = [specified[i : i + _GROUP] for i in range(0, len(specified), _GROUP)]
    last_letters = tuple(lab for _, lab in groups[-1])
    series = "D" if last_letters[-1] == "R" else "L"
    anomer = ""
    if skeleton.ring and anomeric_label:
        anomer = "α-" if (series == "D") == (anomeric_label == "S") else "β-"
    if full and skeleton.n in (4, 5, 6):
        letters = tuple(lab for _, lab in specified)
        letters = letters if letters[-1] == "R" else _mirror(letters)
        stem = table.get(letters)
        if stem is None:
            return None
        retained = stem
        dropped = False
        if six_deoxy_only and not skeleton.ketose and skeleton.n == 6 and (series, stem) in _SIX_DEOXY:
            retained, dropped = _SIX_DEOXY[(series, stem)], True
        name = f"{series}-{retained}" if not skeleton.ring else f"{anomer}{series}-{retained[:-2]}{ring_word}"
        return name, dropped
    chunks = []
    for group in groups:
        letters = tuple(lab for _, lab in group)
        direction = "D" if letters[-1] == "R" else "L"
        letters = letters if direction == "D" else _mirror(letters)
        prefix = _PREFIX.get(letters)
        if prefix is None:
            return None
        chunks.append((direction, prefix))
    cited = [f"{d}-{p}" for d, p in reversed(chunks)]
    if skeleton.ring:
        cited[-1] = anomer + cited[-1]
    number = _NUMBER_STEM[skeleton.n]
    if skeleton.ketose:
        stem = f"{number}-2-ulo{ring_word}" if skeleton.ring else f"{number}-2-ulose"
    else:
        stem = f"{number}o{ring_word}" if skeleton.ring else f"{number}ose"
    return "-".join(cited) + "-" + stem, False


def _plain_group(mol, graph, root, origin):
    """True for a group of carbon, halogen and ether or hydroxy oxygen only: nothing that outranks the sugar's hydroxy class."""
    for a in subtree(graph, root, origin):
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() not in (6, 8, 9, 17, 35, 53):
            return False
        if atom.GetAtomicNum() == 8 and (atom.GetDegree() > 2 or any(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds())):
            return False
    return True


def _enclose(name, compound):
    return f"({name})" if compound and not name.startswith("(") else name


def _segment(locants, name, kind, compound):
    count = len(locants)
    joined = ",".join(map(str, locants))
    multiplier = "" if count == 1 else (_COMPLEX[count] if compound else numerical_term(count))
    group = _enclose(name, compound)
    if kind:
        return f"{joined}-{multiplier}-{kind}-{group}" if multiplier else f"{joined}-{kind}-{group}"
    return f"{joined}-{multiplier}{group}"


def _oxo_anion(atom):
    """The O(-) of a phosphate or sulfate ester group, cited in the ester word as 'phosphate' or 'sulfate'."""
    return (
        atom.GetAtomicNum() == 8
        and atom.GetFormalCharge() == -1
        and atom.GetDegree() == 1
        and atom.GetNeighbors()[0].GetAtomicNum() in (15, 16)
    )


def _substituted_sugar_name(mol, bridges=()):
    if any(a.GetAtomicNum() not in _ELEMENTS or (a.GetFormalCharge() and not _oxo_anion(a)) or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    skeleton = _ring_skeleton(mol, graph) or _chain_skeleton(mol, graph)
    if skeleton is None or skeleton.n < 4 or len(set(skeleton.carbons)) != skeleton.n:
        return None
    skeleton_atoms = set(skeleton.carbons) | ({skeleton.hetero} if skeleton.ring else set())
    decorations, reached = {}, set(skeleton_atoms)
    anomeric_kind = "OH"
    for carbon in skeleton.carbons:
        position = skeleton.position(carbon)
        if skeleton.ring and position == skeleton.closing:
            if [n for n in graph[carbon] if n not in skeleton_atoms]:
                return None
            continue
        kind, x, root = _classify(mol, graph, carbon, skeleton_atoms, skeleton)
        if kind == "oxo":
            reached.add(x)
            continue
        if skeleton.ring and position == skeleton.anomeric and kind not in ("OH", "ester", "ether", "X", "N", "SH", "Sether"):
            return None
        at_anomeric = skeleton.ring and position == skeleton.anomeric
        if (kind == "Sether" and not at_anomeric) or (kind == "Csub" and (at_anomeric or position in (1, skeleton.n))):
            return None
        if kind == "C" and position in (1, skeleton.n):
            return None
        decorations[carbon] = (kind, x, root)
        if skeleton.ring and position == skeleton.anomeric:
            anomeric_kind = kind
        if x is not None:
            reached |= subtree(graph, x, carbon) if kind != "C" else {x} | subtree(graph, root, carbon)
    if reached != set(range(mol.GetNumAtoms())):
        return None
    ring_replacement = _CHALCOGEN_WORDS.get(mol.GetAtomWithIdx(skeleton.hetero).GetAtomicNum()) if skeleton.ring else None
    if sum(kind in ("OH", "ether", "ester") for kind, _, _ in decorations.values()) < 2 - (anomeric_kind in ("SH", "Sether")):
        return None
    if all(kind == "OH" for kind, _, _ in decorations.values()) and not ring_replacement and not skeleton.ring and not bridges:
        return None
    try:
        restored, atom_of = _restore(mol, graph, skeleton, decorations)
        centres, anomeric_label = _labels(restored, skeleton, atom_of)
        deoxy_six = skeleton.n == 6 and not skeleton.ketose and decorations[skeleton.carbons[5]][0] == "H"
        found = _parent(skeleton, centres, anomeric_label, deoxy_six)
        if found is None:
            return None
        core, dropped = found
        entries, esters = {}, {}
        aglycone = None
        for carbon, (kind, x, root) in decorations.items():
            position = skeleton.position(carbon)
            if skeleton.ring and position == skeleton.anomeric and kind in ("ether", "X", "N", "Sether"):
                if kind in ("ether", "Sether"):
                    if not _plain_group(mol, graph, root, x):
                        return None
                    aglycone = cited_group(mol, graph, root, x)[0]
                    if kind == "Sether":
                        entries.setdefault((_CHALCOGEN_WORDS[mol.GetAtomWithIdx(x).GetAtomicNum()],) * 2 + ("", False), []).append(position)
                elif kind == "N" and not (mol.GetAtomWithIdx(x).GetDegree() == 1 and mol.GetAtomWithIdx(x).GetTotalNumHs() == 2):
                    return None
                continue
            if kind == "H":
                if not (position == 6 and dropped):
                    entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
            elif kind == "N":
                name, compound = cited_group(mol, graph, x, carbon)
                entries.setdefault((name, name, "", compound), []).append(position)
                entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
            elif kind == "X":
                text = _HALOGENS[mol.GetAtomWithIdx(x).GetAtomicNum()]
                entries.setdefault((text, text, "", False), []).append(position)
                entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
            elif kind in ("SH", "thione"):
                word = _CHALCOGEN_WORDS[mol.GetAtomWithIdx(x).GetAtomicNum()]
                entries.setdefault((word, word, "", False), []).append(position)
            elif kind == "Csub":
                if not _plain_group(mol, graph, x, carbon):
                    return None
                name, compound = cited_group(mol, graph, x, carbon)
                entries.setdefault((name, name, "", compound), []).append(position)
                entries.setdefault(("deoxy", "deoxy", "", False), []).append(position)
            elif kind == "ether":
                if not _plain_group(mol, graph, root, x):
                    return None
                name, compound = cited_group(mol, graph, root, x)
                entries.setdefault((name, name, "O", compound), []).append(position)
            elif kind == "C":
                if not _plain_group(mol, graph, root, carbon):
                    return None
                name, compound = cited_group(mol, graph, root, carbon)
                entries.setdefault((name, name, "C", compound), []).append(position)
            elif kind == "ester" and anomeric_kind in ("X", "N"):
                name, compound = cited_group(mol, graph, root, x)
                entries.setdefault((name, name, "O", compound), []).append(position)
            elif kind == "ester":
                anion = ester_anion(mol, x, root)
                if anion is None:
                    return None
                esters.setdefault(anion, []).append(position)
        if ring_replacement:
            entries.setdefault((ring_replacement, ring_replacement, "", False), []).append(skeleton.closing)
    except UnsupportedStructure:
        return None
    cited = [(key[0], _segment(sorted(locants), key[0], key[2], key[3])) for key, locants in entries.items()]
    for first, second in bridges:
        if first not in skeleton.carbons or second not in skeleton.carbons:
            return None
        cited.append(("anhydro", "{},{}-anhydro".format(*sorted((skeleton.position(first), skeleton.position(second))))))
    segments = [text for _, text in sorted(cited, key=lambda item: (alpha_sort_key(item[0]), item[1]))]
    prefix = "-".join(segments) + "-" if segments else ""
    if anomeric_kind in ("ether", "Sether"):
        return f"{aglycone} {prefix}{core[:-1]}ide" + (f" {ester_words([(a, sorted(l)) for a, l in esters.items()])}" if esters else "")
    if anomeric_kind == "X":
        halide = _HALIDE_WORDS[mol.GetAtomWithIdx(decorations[skeleton.carbons[skeleton.anomeric - 1]][1]).GetAtomicNum()]
        return f"{prefix}{core[:-1]}yl {halide}"
    if anomeric_kind == "N":
        return f"{prefix}{core[:-1]}ylamine"
    name = prefix + core
    return f"{name} {ester_words([(a, sorted(l)) for a, l in esters.items()])}" if esters else name



def _opened_ether(mol, oxygen, carbon):
    """`mol` with the bond between the ether `oxygen` and `carbon` replaced by a hydroxy group on `carbon`, the
    configuration of `carbon` kept."""
    editable = Chem.RWMol(mol)
    atom = editable.GetAtomWithIdx(carbon)
    before = [n.GetIdx() for n in atom.GetNeighbors()]
    tag = atom.GetChiralTag()
    editable.RemoveBond(oxygen, carbon)
    hydroxy = editable.AddAtom(Chem.Atom(8))
    editable.AddBond(carbon, hydroxy, Chem.BondType.SINGLE)
    after = [n.GetIdx() for n in editable.GetAtomWithIdx(carbon).GetNeighbors()]
    order = [hydroxy if n == oxygen else n for n in before]
    swaps = sum(1 for i in range(len(order)) for j in range(i + 1, len(order)) if after.index(order[i]) > after.index(order[j]))
    if tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW) and swaps % 2:
        editable.GetAtomWithIdx(carbon).SetChiralTag(
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW else Chem.ChiralType.CHI_TETRAHEDRAL_CW
        )
    opened = editable.GetMol()
    for a in opened.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(opened)
    return opened


def _anhydro_name(mol):
    """P-102.5.6.7.1: an ether between two skeleton carbons of a monosaccharide (1,6-anhydro-β-D-glucopyranose): the
    ether is opened into two hydroxy groups, the sugar named, and 'anhydro' cited with the two locants."""
    names = set()
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 2 or not oxygen.IsInRing() or oxygen.GetIsAromatic():
            continue
        ends = [n.GetIdx() for n in oxygen.GetNeighbors()]
        if any(mol.GetAtomWithIdx(e).GetAtomicNum() != 6 for e in ends):
            continue
        for kept, opened_end in (ends, ends[::-1]):
            try:
                opened = _opened_ether(mol, oxygen.GetIdx(), opened_end)
                found = _substituted_sugar_name(opened, bridges=[(kept, opened_end)])
            except (UnsupportedStructure, ValueError, RuntimeError, Chem.rdchem.MolSanitizeException):
                continue
            if found is not None:
                names.add(found)
    return next(iter(names)) if len(names) == 1 else None


def substituted_sugar_name(mol):
    errors = (UnsupportedStructure, ValueError, RuntimeError, Chem.rdchem.MolSanitizeException)
    try:
        found = _substituted_sugar_name(mol)
    except errors:
        found = None
    if found is not None:
        return found
    try:
        return _anhydro_name(mol)
    except errors:
        return None


def has_substituted_sugar_shape(mol) -> bool:
    return substituted_sugar_name(mol) is not None


def name_substituted_sugar(mol) -> str:
    return substituted_sugar_name(mol)


def sugar_substituent_group(mol, graph, root, coming_from):
    """(name, True) of a monosaccharide cited as a substituent group through an atom other than its anomeric carbon
    (P-102.6.2): the sugar is named with a hydrogen atom in place of the bond to the parent, its final e is dropped
    and the position of the free valence is cited before 'yl', 'O-yl' or 'C-yl'."""
    from .core import smiles_to_iupac

    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() not in (6, 8) or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        return None
    atoms = subtree(graph, root, coming_from)
    if coming_from in atoms or not 5 <= len(atoms) <= 40 or not any(mol.GetAtomWithIdx(a).IsInRing() for a in atoms):
        return None
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() not in _ELEMENTS or mol.GetAtomWithIdx(a).GetIsAromatic() for a in atoms):
        return None
    editable = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - atoms - {coming_from}, reverse=True):
        editable.RemoveAtom(index)
    for a in editable.GetAtoms():
        a.SetIntProp("_orig", -1)
    original = sorted(atoms | {coming_from})
    for new, old in enumerate(original):
        editable.GetAtomWithIdx(new).SetIntProp("_orig", old)
    editable.GetAtomWithIdx(original.index(coming_from)).SetAtomicNum(1)
    editable.GetAtomWithIdx(original.index(coming_from)).SetIsAromatic(False)
    try:
        sugar = editable.GetMol()
        Chem.SanitizeMol(sugar)
        sugar = Chem.RemoveHs(sugar)
        sugar_graph = adjacency(sugar)
        if sugar.GetRingInfo().NumRings() > 1:
            return None
        skeleton = _ring_skeleton(sugar, sugar_graph) or _chain_skeleton(sugar, sugar_graph)
        if skeleton is None or skeleton.n < 4:
            return None
        index = {a.GetIntProp("_orig"): a.GetIdx() for a in sugar.GetAtoms()}
        carbon = index[root] if atom.GetAtomicNum() == 6 else next(
            n.GetIdx() for n in sugar.GetAtomWithIdx(index[root]).GetNeighbors() if n.GetIdx() in skeleton.carbons
        )
        if carbon not in skeleton.carbons:
            return None
        position = skeleton.position(carbon)
        if skeleton.ring and position == skeleton.anomeric:
            return None
        if atom.GetAtomicNum() == 8:
            kind = "O"
        else:
            skeleton_atoms = set(skeleton.carbons) | ({skeleton.hetero} if skeleton.ring else set())
            exo = [n for n in sugar_graph[carbon] if n not in skeleton_atoms]
            kind = "C" if exo else ""
        name = smiles_to_iupac(Chem.MolToSmiles(sugar))
    except (UnsupportedStructure, ValueError, RuntimeError, KeyError, StopIteration, Chem.rdchem.MolSanitizeException):
        return None
    if not name.endswith("e") or name.count(" ") > 1 or ("(" in name and " " in name) or len(name) < 8:
        return None
    if " " in name and not name.endswith("oside"):
        return None
    marker = f"-{position}-{kind}-yl" if kind else f"-{position}-yl"
    from ._substituents import BRANCH_STEREO

    context = BRANCH_STEREO.get()
    if context:
        context["used"].update(("atom", a) for a in atoms)
    return f"{name[:-1]}{marker}", True
