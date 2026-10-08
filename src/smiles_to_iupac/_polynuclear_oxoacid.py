"""Di- and polynuclear noncarbon oxoacids with their derivatives (P-67.2, P-67.3): a chain of acid centres (P, As, Sb, B,
Si, S, Se, Te) joined directly ('hypo') or through one bridging atom takes the retained polyacid name, functional replacement
is cited by prefixes with chain locants, esters and organic groups by words, and chains that match no retained acid are
named on the senior mononuclear acid (P-67.3.1) or as dithioxane-diones (P-67.3.2)."""

from dataclasses import dataclass, field

from rdkit import Chem

from ._alkoxy import alkoxy_prefix
from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._multiplicative_text import enclose
from ._numerals import numerical_term
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_CENTER_ELEMENTS = {5, 14, 15, 16, 33, 34, 51, 52}
_CHALCOGENS = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_HALOGENS = {9: "F", 17: "Cl", 35: "Br", 53: "I"}
_ELEMENT_ORDER = {"O": 0, "S": 1, "Se": 2, "Te": 3}
_REPLACEMENT_WORD = {"S": "thio", "Se": "seleno", "Te": "telluro", "NH": "imido"}
_PHOSPHORUS_STEMS = {15: ("phosphor", "phosphon"), 33: ("arsor", "arson"), 51: ("stibor", "stibon")}
_SLOTS = {15: 3, 33: 3, 51: 3, 5: 3, 14: 4, 16: 2, 34: 2, 52: 2}
_CHALCOGEN_ACID = {16: "sulfuric", 34: "selenic", 52: "telluric"}
_IRREGULAR_ANIONS = (("phosphoric", "phosphate"), ("phosphorous", "phosphite"), ("sulfuric", "sulfate"), ("sulfurous", "sulfite"))
_PSEUDOHALIDES = {
    "CN": ("cyano", "cyanide", 1),
    "NCO": ("isocyanato", "isocyanate", 2),
    "NCS": ("isothiocyanato", "isothiocyanate", 3),
    "NCSe": ("isoselenocyanato", "isoselenocyanate", 4),
    "NCTe": ("isotellurocyanato", "isotellurocyanate", 5),
}
_HALIDES = {"F": ("fluoro", "fluoride"), "Cl": ("chloro", "chloride"), "Br": ("bromo", "bromide"), "I": ("iodo", "iodide")}


@dataclass
class Ligand:
    atom: int
    element: str
    kind: str = "hydroxy"
    oxo: bool = False
    group: object = None
    extra: tuple = ()


@dataclass
class Center:
    atom: int
    z: int
    hydrogens: int
    oxo: list = field(default_factory=list)
    singles: list = field(default_factory=list)


@dataclass
class Link:
    elements: tuple


def _is_center(atom):
    z = atom.GetAtomicNum()
    if z not in _CENTER_ELEMENTS or atom.GetFormalCharge() or atom.GetIsotope() or atom.IsInRing():
        return False
    mol = atom.GetOwningMol()
    neighbors = atom.GetNeighbors()
    if any(n.GetAtomicNum() == 6 and not _is_cyanide_carbon(mol, n, atom) for n in neighbors):
        return False
    if z in (16, 34, 52):
        return any(
            (n.GetDegree() == 1 or (n.GetAtomicNum() == 7 and n.GetDegree() == 2))
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in neighbors
        )
    return any(n.GetAtomicNum() in _CHALCOGENS or n.GetAtomicNum() in (7, *_HALOGENS) or n.GetAtomicNum() == 6 for n in neighbors)


def _is_cyanide_carbon(mol, carbon, center):
    return (
        carbon.GetDegree() == 2
        and not carbon.GetTotalNumHs()
        and any(
            n.GetAtomicNum() == 7 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 3.0
            for n in carbon.GetNeighbors()
        )
    )


def _scan(mol):
    """(centres, links) of a linear chain of acid centres whose other groups are hydrocarbyl or halogen, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    center_atoms = [a.GetIdx() for a in mol.GetAtoms() if _is_center(a)]
    if len(center_atoms) < 2:
        return None
    center_set = set(center_atoms)
    graph = adjacency(mol)
    centers = {}
    adjacent = {c: [] for c in center_atoms}
    seen_links = set()
    for c in center_atoms:
        atom = mol.GetAtomWithIdx(c)
        center = Center(c, atom.GetAtomicNum(), atom.GetTotalNumHs())
        for n in graph[c]:
            if not _classify(mol, graph, center_set, center, n, adjacent, seen_links):
                return None
        centers[c] = center
    ordered = _order_chain(centers, adjacent)
    if ordered is None:
        return None
    skeleton = set(centers)
    for center in centers.values():
        for lig in center.oxo + center.singles:
            skeleton.add(lig.atom)
            skeleton.update(lig.extra)
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in skeleton and atom.GetAtomicNum() not in (6, *_HALOGENS) and not _bridge_atom(atom, center_set):
            return None
    return ordered


def _bridge_atom(atom, center_set):
    return atom.GetAtomicNum() in (7, *_CHALCOGENS) and any(
        n.GetIdx() in center_set or n.GetAtomicNum() in _CHALCOGENS for n in atom.GetNeighbors()
    )


def _classify(mol, graph, center_set, center, n, adjacent, seen_links):
    """Sort neighbour `n` of `center` into an oxo, single or link entry; False when it is none of them."""
    c = center.atom
    neighbor = mol.GetAtomWithIdx(n)
    z = neighbor.GetAtomicNum()
    order = mol.GetBondBetweenAtoms(c, n).GetBondTypeAsDouble()
    if n in center_set:
        if order != 1.0:
            return False
        _add_link(adjacent, seen_links, c, n, ())
        return True
    if neighbor.GetFormalCharge() or neighbor.GetIsotope():
        return False
    if z in _HALOGENS and neighbor.GetDegree() == 1:
        center.singles.append(Ligand(n, _HALOGENS[z], kind="halide"))
    elif z in _CHALCOGENS and order == 2.0 and neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 0:
        center.oxo.append(Ligand(n, _CHALCOGENS[z], oxo=True))
    elif z == 7 and order == 2.0 and neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 1:
        center.oxo.append(Ligand(n, "NH", oxo=True))
    elif z == 7 and order == 2.0 and neighbor.GetDegree() == 2 and neighbor.GetTotalNumHs() == 0:
        others = [m for m in graph[n] if m != c]
        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(n, others[0]).GetBondTypeAsDouble() != 1.0:
            return False
        center.oxo.append(Ligand(n, "NH", oxo=True, group=others[0]))
    elif z == 6 and order == 1.0 and _is_cyanide_carbon(mol, neighbor, mol.GetAtomWithIdx(c)):
        center.singles.append(Ligand(n, "CN", kind="pseudo", extra=tuple(m for m in graph[n] if m != c)))
    elif order != 1.0:
        return False
    elif z in _CHALCOGENS or z == 7:
        parsed = _single_or_bridge(mol, graph, center_set, c, n)
        if parsed is None:
            return False
        kind, payload = parsed
        if kind == "single":
            center.singles.append(payload)
        else:
            _add_link(adjacent, seen_links, c, payload[0], payload[1])
    else:
        return False
    return True


def _add_link(adjacent, seen_links, a, b, elements):
    key = frozenset((a, b))
    if key in seen_links:
        return
    seen_links.add(key)
    adjacent[a].append((b, elements))
    adjacent[b].append((a, tuple(reversed(elements))))


def _single_or_bridge(mol, graph, center_set, c, n):
    atom = mol.GetAtomWithIdx(n)
    z = atom.GetAtomicNum()
    others = [m for m in graph[n] if m != c]
    if not others:
        if z in _CHALCOGENS and atom.GetTotalNumHs() == 1:
            return "single", Ligand(n, _CHALCOGENS[z], kind="hydroxy")
        if z == 7 and atom.GetTotalNumHs() == 2:
            return "single", Ligand(n, "NH2", kind="amide")
        return None
    if len(others) != 1:
        return None
    m = others[0]
    other = mol.GetAtomWithIdx(m)
    if atom.GetFormalCharge() or other.GetFormalCharge():
        return None
    order = mol.GetBondBetweenAtoms(n, m).GetBondTypeAsDouble()
    if z == 7 and order == 2.0 and other.GetAtomicNum() == 6:
        return _isocyanate(mol, graph, n, m)
    if order != 1.0:
        return None
    if z in _CHALCOGENS and other.GetAtomicNum() == 6:
        return "single", Ligand(n, _CHALCOGENS[z], kind="ester", group=m)
    if z == 7 and other.GetAtomicNum() == 7 and other.GetDegree() == 1 and other.GetTotalNumHs() == 2 and atom.GetTotalNumHs() == 1:
        return "single", Ligand(n, "NHNH2", kind="hydrazide", extra=(m,))
    if m in center_set:
        element = _CHALCOGENS.get(z, "NH" if z == 7 and atom.GetTotalNumHs() == 1 else None)
        return None if element is None else ("link", (m, (element,)))
    if z in _CHALCOGENS and other.GetAtomicNum() in _CHALCOGENS and other.GetDegree() == 2 and not atom.GetTotalNumHs():
        far = [x for x in graph[m] if x != n]
        if len(far) == 1 and far[0] in center_set and not other.GetTotalNumHs():
            return "link", (far[0], (_CHALCOGENS[z], _CHALCOGENS[other.GetAtomicNum()]))
    return None


def _isocyanate(mol, graph, n, carbon):
    ends = [x for x in graph[carbon] if x != n]
    if len(ends) != 1 or mol.GetBondBetweenAtoms(carbon, ends[0]).GetBondTypeAsDouble() != 2.0:
        return None
    end = mol.GetAtomWithIdx(ends[0])
    element = _CHALCOGENS.get(end.GetAtomicNum())
    if element is None or end.GetDegree() != 1 or end.GetTotalNumHs():
        return None
    return "single", Ligand(n, f"NC{element}", kind="pseudo", extra=(carbon, ends[0]))


def _order_chain(centers, adjacent):
    ends = [c for c in adjacent if len(adjacent[c]) == 1]
    if len(ends) != 2 or any(len(v) > 2 or not v for v in adjacent.values()):
        return None
    order, links, previous = [ends[0]], [], None
    while True:
        onward = [(n, el) for n, el in adjacent[order[-1]] if n != previous]
        if not onward:
            break
        n, elements = onward[0]
        links.append(Link(elements))
        previous = order[-1]
        order.append(n)
    if len(order) != len(centers):
        return None
    return [centers[c] for c in order], links


def _signature(center):
    return center.z, len(center.oxo), center.hydrogens


def _retained_base(center, n, direct):
    """The retained name of the polyacid made of `n` such centres, without 'acid' (P-67.2.1)."""
    z, oxo, hydrogens = _signature(center)
    multiplier = numerical_term(n)
    hypo = "hypo" if direct else ""
    if z in _PHOSPHORUS_STEMS and hydrogens in (0, 1) and oxo in (0, 1):
        return f"{hypo}{multiplier}{_PHOSPHORUS_STEMS[z][hydrogens]}{'ic' if oxo else 'ous'}"
    if z in _CHALCOGEN_ACID and oxo == 2:
        if direct and z == 16 and n == 2:
            return "dithionic"
        if direct and n > 2:
            return None
        return f"{hypo}{multiplier}{_CHALCOGEN_ACID[z]}"
    if z == 16 and oxo == 1 and direct and n == 2:
        return "dithionous"
    if z == 5 and oxo == 0 and hydrogens in (0, 1):
        return f"{hypo}{multiplier}{'boric' if hydrogens == 0 else 'boronic'}"
    if z == 14 and oxo == 0 and hydrogens == 0:
        return f"{hypo}{multiplier}silicic"
    return None


def _positions(centers, links, reverse):
    """{centre index: locant}, {link index: locant or None} numbering the chain from one end to the other (P-67.2.2.1)."""
    order = list(range(len(centers)))
    link_order = list(range(len(links)))
    if reverse:
        order.reverse()
        link_order.reverse()
    center_pos, link_pos, pos = {}, {}, 1
    for i, c in enumerate(order):
        center_pos[c] = pos
        pos += 1
        if i < len(link_order):
            li = link_order[i]
            link_pos[li] = pos if links[li].elements else None
            pos += 1 if links[li].elements else 0
    return center_pos, link_pos


def _link_replacement(elements):
    if elements in ((), ("O",)):
        return None
    if elements == ("O", "O"):
        return "peroxy"
    if elements == ("S", "S"):
        return "(dithioperoxy)"
    if len(elements) == 2:
        return "(thioperoxy)" if set(elements) == {"O", "S"} else None
    return _REPLACEMENT_WORD.get(elements[0])


def _class_rank(lig):
    """Seniority of an acid derivative class (P-67.2.3, P-67.2.4): halides, then pseudohalides, amides, hydrazides."""
    if lig.kind == "halide":
        return 0
    if lig.kind == "pseudo":
        return _PSEUDOHALIDES[lig.element][2]
    return {"amide": 10, "hydrazide": 11}.get(lig.kind)


def _class_word(lig):
    if lig.kind == "halide":
        return _HALIDES[lig.element][1]
    if lig.kind == "pseudo":
        return _PSEUDOHALIDES[lig.element][1]
    return lig.kind


def _prefix_word(lig):
    if lig.kind == "halide":
        return _HALIDES[lig.element][0]
    if lig.kind == "pseudo":
        return _PSEUDOHALIDES[lig.element][0]
    return {"amide": "amido", "hydrazide": "hydrazido"}[lig.kind]


def _senior_class(centers):
    ligands = [lig for center in centers for lig in center.singles]
    if any(lig.kind in ("hydroxy", "ester") for lig in ligands):
        return None
    ranks = [_class_rank(lig) for lig in ligands]
    return min(ranks) if ranks else None


def _replacements(centers, links, center_pos, link_pos, senior):
    """[(prefix word, locant)] of the functional replacements; members of the cited class are left out."""
    items = []
    for ci, center in enumerate(centers):
        for lig in center.oxo:
            if lig.element != "O":
                items.append((_REPLACEMENT_WORD[lig.element], center_pos[ci]))
        for lig in center.singles:
            if lig.kind in ("hydroxy", "ester"):
                if lig.element != "O":
                    items.append((_REPLACEMENT_WORD[lig.element], center_pos[ci]))
            elif _class_rank(lig) != senior or senior is None:
                items.append((_prefix_word(lig), center_pos[ci]))
    for li, link in enumerate(links):
        word = _link_replacement(link.elements)
        if word is not None:
            items.append((word, link_pos[li]))
    return items


def _n_substituents(mol, graph, halogens, centers, center_pos):
    """[(name, compound, locant text)] of the groups on imido nitrogens (P-67.2.2.2: N1, N2)."""
    found = []
    for ci, center in enumerate(centers):
        for lig in center.oxo:
            if lig.group is not None:
                name, compound = name_branch(graph, lig.group, lig.atom, halogens, mol=mol)
                found.append((name, compound, f"N{center_pos[ci]}"))
    return found


def _format_n_substituents(found):
    grouped = {}
    for name, compound, locant in found:
        grouped.setdefault((name, compound), []).append(locant)
    parts = []
    for (name, compound), locants in sorted(grouped.items(), key=lambda kv: kv[0][0]):
        count = len(locants)
        multiplied = ((_bis(count) if compound else numerical_term(count)) if count > 1 else "") + (f"({name})" if compound and count > 1 else name)
        parts.append(f"{','.join(sorted(locants))}-{multiplied}")
    return "-".join(parts) + "-" if parts else ""


def _format_replacements(items):
    grouped = {}
    for word, locant in items:
        grouped.setdefault(word, []).append(locant)
    parts = []
    for word in sorted(grouped, key=lambda w: w.strip("()")):
        locants = sorted(grouped[word])
        multiplied = (numerical_term(len(locants)) if len(locants) > 1 else "") + word
        parts.append(f"{','.join(str(x) for x in locants)}-{multiplied}")
    return "-".join(parts)


def _plain_prefixes(items):
    counts = {}
    for word, _ in items:
        counts[word] = counts.get(word, 0) + 1
    return "".join(
        (numerical_term(counts[w]) if counts[w] > 1 else "") + w for w in sorted(counts, key=lambda w: w.strip("()"))
    )


def _descriptors(centers, center_pos):
    """Italic letters for the hydrogen-bearing positions of the centres carrying a replaced S, Se or Te (P-67.1.2.4.1)."""
    cited = []
    for ci, center in enumerate(centers):
        chalcogen_replaced = any(lig.element in ("S", "Se", "Te") for lig in center.oxo) or any(
            lig.element in ("S", "Se", "Te") for lig in center.singles if lig.kind in ("hydroxy", "ester")
        )
        if not chalcogen_replaced:
            continue
        hydroxyls = [lig for lig in center.singles if lig.kind == "hydroxy"]
        for lig in [lig for lig in hydroxyls if lig.element != "O"] or hydroxyls:
            cited.append((_ELEMENT_ORDER[lig.element], lig.element, center_pos[ci]))
    return ",".join(f"{element}{pos}" for _, element, pos in sorted(cited)) + "-" if cited else ""


def _function_locants(centers, center_pos, senior):
    """Locants of the centres that carry the principal function: acid groups, or the cited class of acid derivatives
    (P-67.2.2.1: the principal function receives the lower locant)."""
    def carries(lig):
        return lig.kind in ("hydroxy", "ester") if senior is None else _class_rank(lig) == senior

    return tuple(sorted(center_pos[ci] for ci, c in enumerate(centers) if any(carries(lig) for lig in c.singles)))


def _anion(base):
    for acid, anion in _IRREGULAR_ANIONS:
        if base.endswith(acid):
            return base[: -len(acid)] + anion
    if base.endswith("ic"):
        return base[:-2] + "ate"
    if base.endswith("ous"):
        return base[:-3] + "ite"
    raise UnsupportedStructure("a polyacid name without an 'ic' or 'ous' ending has no anion name")


def _bis(count):
    return {2: "bis", 3: "tris", 4: "tetrakis"}.get(count, "bis")


def _homogeneous(centers, links):
    signatures = {_signature(c) for c in centers}
    direct = {not link.elements for link in links}
    if len(signatures) != 1 or len(direct) != 1:
        return None
    for i, center in enumerate(centers):
        bridges = (1 if i > 0 else 0) + (1 if i < len(centers) - 1 else 0)
        if len(center.singles) + center.hydrogens + bridges != _SLOTS[center.z]:
            return None
    return direct.pop()


def _chain_name(mol, centers, links, priority=False):
    direct = _homogeneous(centers, links)
    if direct is None:
        return None
    base = _retained_base(centers[0], len(centers), direct)
    if base is None:
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    senior = _senior_class(centers)
    candidates = []
    for reverse in (False, True):
        center_pos, link_pos = _positions(centers, links, reverse)
        items = _replacements(centers, links, center_pos, link_pos, senior)
        by_word = {}
        for word, locant in items:
            by_word.setdefault(word, []).append(locant)
        key = (
            _function_locants(centers, center_pos, senior),
            tuple(sorted(loc for _, loc in items)),
            tuple(tuple(sorted(by_word[w])) for w in sorted(by_word, key=lambda w: w.strip("()"))),
        )
        candidates.append((key, center_pos, items))
    _, center_pos, items = min(candidates, key=lambda c: c[0])
    n_found = _n_substituents(mol, graph, halogens, centers, center_pos)
    if priority and not items and senior is None and not n_found and centers[0].z != 5:
        raise UnsupportedStructure("a plain homogeneous polyacid or ester is named by the part-wise modules")
    n_prefix = _format_n_substituents(n_found)
    if senior is not None:
        return n_prefix + _class_name(centers, center_pos, items, base, senior)
    esters, hydrogens = [], 0
    for ci, center in enumerate(centers):
        for lig in center.singles:
            if lig.kind == "ester":
                name, compound = name_branch(graph, lig.group, lig.atom, halogens, mol=mol)
                esters.append((name, compound, lig.element, center_pos[ci]))
            elif lig.kind == "hydroxy":
                hydrogens += 1
    if not esters:
        prefix = n_prefix + _format_replacements(items)
        return f"{prefix}{base} {_descriptors(centers, center_pos)}acid"
    return _ester_name(esters, hydrogens, items, base, n_prefix)


def _class_name(centers, center_pos, items, base, senior):
    members = []
    for ci, center in enumerate(centers):
        for lig in center.singles:
            if _class_rank(lig) == senior:
                members.append((_class_word(lig), center_pos[ci]))
    names = sorted({m for m, _ in members})
    if len(names) == 1:
        count = len(members)
        words = f"{numerical_term(count) if count > 1 else ''}{names[0]}"
    else:
        parts = []
        for name in names:
            locants = sorted(loc for m, loc in members if m == name)
            parts.append(f"{','.join(str(x) for x in locants)}-{numerical_term(len(locants)) if len(locants) > 1 else ''}{name}")
        words = " ".join(parts)
    return f"{_format_replacements(items)}{base} {words}"


def _ester_name(esters, hydrogens, items, base, n_prefix=""):
    replaced_chalcogen = any(word in ("thio", "seleno", "telluro") for word, _ in items)
    elements = [e for _, _, e, _ in esters]
    cite_locants = replaced_chalcogen and len(set(elements)) < len(elements)
    cite_letters = replaced_chalcogen and len(set(elements)) == len(elements) and len(esters) > 1
    grouped = {}
    for name, compound, element, pos in esters:
        grouped.setdefault((name, compound), []).append((element, pos))
    words = []
    for (name, compound), where in sorted(grouped.items(), key=lambda kv: kv[0][0]):
        count = len(where)
        multiplied = ((numerical_term(count) if not compound else _bis(count)) if count > 1 else "")
        shown = f"({name})" if compound and count > 1 else name
        letters = ""
        if cite_locants:
            letters = ",".join(f"{e}{p}" for e, p in sorted(where, key=lambda w: (_ELEMENT_ORDER[w[0]], w[1]))) + "-"
        elif cite_letters:
            letters = ",".join(e for e, _ in sorted(where, key=lambda w: _ELEMENT_ORDER[w[0]])) + "-"
        words.append(f"{letters}{multiplied}{shown}")
    if hydrogens:
        words.append((numerical_term(hydrogens) if hydrogens > 1 else "") + "hydrogen")
    return " ".join(words) + " " + n_prefix + _plain_prefixes(items) + _anion(base)


def _prefix_entry(mol, graph, halogens, lig):
    if lig.kind == "hydroxy" and lig.element == "O":
        return "hydroxy", False
    if lig.kind == "hydroxy":
        return {"S": "sulfanyl", "Se": "selanyl", "Te": "tellanyl"}[lig.element], False
    if lig.kind == "ester" and lig.element == "O":
        name, compound = name_branch(graph, lig.group, lig.atom, halogens, mol=mol)
        return alkoxy_prefix(name, compound)
    if lig.kind == "halide":
        return _HALIDES[lig.element][0], False
    if lig.kind == "amide":
        return "amino", False
    return None


def _dithioxane_name(mol, centers, links):
    """P-67.3.2: sulfur(IV) centres joined through oxygen are an oxane-type chain, not 'disulfurous acid'."""
    if any(c.z != 16 or len(c.oxo) != 1 for c in centers) or any(link.elements != ("O",) for link in links):
        return None
    if len({c.oxo[0].element for c in centers}) != 1:
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    n = len(centers)
    best = None
    for reverse in (False, True):
        order = list(reversed(centers)) if reverse else list(centers)
        entries = []
        for i, center in enumerate(order):
            for lig in center.singles:
                entry = _prefix_entry(mol, graph, halogens, lig)
                if entry is None:
                    return None
                entries.append((entry[0], entry[1], 2 * i + 1))
        key = tuple(loc for _, _, loc in sorted(entries, key=lambda e: e[0].strip("([{").lower()))
        if best is None or key < best[0]:
            best = (key, entries)
    grouped = {}
    for name, compound, locant in best[1]:
        info = grouped.setdefault(name, {"locants": [], "compound": compound})
        info["locants"].append(locant)
    locants = ",".join(f"{2 * i + 1}λ4" for i in range(n))
    parent = "dithioxane" if n == 2 else f"{numerical_term(n)}sulfoxane"
    suffix_word = {"O": "one", "S": "thione", "Se": "selone", "Te": "tellone", "NH": "imine"}[centers[0].oxo[0].element]
    suffix_locants = ",".join(str(2 * i + 1) for i in range(n))
    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{'-' if prefix else ''}{locants}-{parent}-{suffix_locants}-{numerical_term(n)}{suffix_word}"


_ELEMENT_RANK = {15: 0, 33: 1, 51: 2, 5: 3, 14: 4}
_ACID_ORDER = {(1, 2): 0, (0, 2): 1, (1, 1): 2, (0, 1): 3}
_STEM_WORD = {
    15: ("phosphin", "phosphon", "phosphor"), 33: ("arsin", "arson", "arsor"), 51: ("stibin", "stibon", "stibor"),
    5: ("borin", "boron", "bor"),
}
_BRIDGE_PREFIX = {"O": "oxy", "S": "sulfanyl", "Se": "selanyl", "Te": "tellanyl"}
_GROUP_STEM = {15: "phosph", 33: "ars", 51: "stib"}


def _parent_acid(center, link_elements):
    """(rank, acid stem, single positions) of `center` as the mononuclear parent of P-67.3.1, or None."""
    if center.z not in _STEM_WORD or len(center.oxo) > 1:
        return None
    oxo = len(center.oxo)
    effective = len(center.singles) + (1 if link_elements == ("O",) and center.hydrogens >= 1 else 0)
    if effective not in (1, 2, 3) or (center.z == 5 and oxo):
        return None
    stem = _STEM_WORD[center.z][effective - 1]
    if center.z == 5:
        return (_ELEMENT_RANK[5], 0, 0, effective), stem + "ic", effective
    kind = 0 if effective in (1, 2) else 1
    order = _ACID_ORDER.get((oxo, effective), effective if kind else 0)
    return (_ELEMENT_RANK[center.z], 0 if oxo else 1, kind, order), stem + ("ic" if oxo else "ous"), effective


def _group_name(mol, graph, halogens, center):
    """Substituent group of an acid centre attached to the parent (P-67.1.4.1, P-67.2.6), or None."""
    z, oxo = center.z, len(center.oxo)
    if any(lig.element != "O" for lig in center.oxo + center.singles) or center.hydrogens > 1:
        return None
    entries = []
    for lig in center.singles:
        entry = _prefix_entry(mol, graph, halogens, lig)
        if entry is None:
            return None
        entries.append(entry)
    two_hydroxyls = len(entries) == 2 and all(e == ("hydroxy", False) for e in entries)
    if z in _GROUP_STEM:
        stem = _GROUP_STEM[z]
        if oxo and center.hydrogens == 0:
            if two_hydroxyls:
                return stem + "ono", False
            return format_mononuclear_prefixes(entries) + stem + "oryl", bool(entries)
        if oxo:
            return format_mononuclear_prefixes(entries) + stem + "onoyl", bool(entries)
        return format_mononuclear_prefixes(entries) + stem + "anyl", bool(entries)
    if z == 5 and not oxo:
        return format_mononuclear_prefixes(entries) + "boranyl", bool(entries)
    if z == 14 and not oxo:
        return format_mononuclear_prefixes(entries) + "silyl", bool(entries)
    if z in _CHALCOGEN_ACID and oxo in (1, 2):
        word = {16: "sulf", 34: "selen", 52: "tellur"}[z]
        sole_hydroxyl = len(entries) == 1 and entries[0] == ("hydroxy", False)
        if oxo == 2:
            return (word + "o", False) if sole_hydroxyl else (format_mononuclear_prefixes(entries) + word + "onyl", bool(entries))
        return (word + "ino", False) if sole_hydroxyl else (format_mononuclear_prefixes(entries) + word + "inyl", bool(entries))
    return None


def _substitutive_name(mol, centers, links):
    """P-67.3.1: the senior mononuclear acid is the parent and the other centre a substituent group."""
    if len(centers) != 2 or links[0].elements not in ((), ("O",), ("S",), ("Se",), ("Te",)):
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    elements = links[0].elements
    choices = []
    for parent, other in ((0, 1), (1, 0)):
        acid = _parent_acid(centers[parent], elements)
        if acid is None or any(lig.element != "O" for lig in centers[parent].oxo + centers[parent].singles):
            continue
        group = _group_name(mol, graph, halogens, centers[other])
        if group is not None:
            choices.append((acid[0], parent, acid, group))
    if not choices:
        return None
    _, parent, acid, group = min(choices, key=lambda c: c[0])
    center = centers[parent]
    bridge = _BRIDGE_PREFIX[elements[0]] if elements else ""
    inner = enclose(group[0]) if group[1] else group[0]
    substituent = enclose(inner + bridge) if bridge else inner
    esters, hydrogens = [], 0
    for lig in center.singles:
        if lig.kind == "ester":
            name, compound = name_branch(graph, lig.group, lig.atom, halogens, mol=mol)
            esters.append((name, compound))
        elif lig.kind == "hydroxy":
            hydrogens += 1
        else:
            return None
    if not esters:
        return f"{substituent}{acid[1]} acid"
    grouped = {}
    for key in esters:
        grouped[key] = grouped.get(key, 0) + 1
    words = []
    for (name, compound), count in sorted(grouped.items(), key=lambda kv: kv[0][0]):
        multiplied = ((_bis(count) if compound else numerical_term(count)) if count > 1 else "")
        words.append(multiplied + (f"({name})" if compound and count > 1 else name))
    if hydrogens:
        words.append((numerical_term(hydrogens) if hydrogens > 1 else "") + "hydrogen")
    return " ".join(words) + " " + substituent + _anion(acid[1])


def name_polynuclear_oxoacid(mol, priority=False):
    """The name of a chain of acid centres. With `priority`, plain homogeneous acids and esters (no replacement, no class
    name) are left to the part-wise ester and acid modules, which name them already."""
    scanned = _scan(mol)
    if scanned is None:
        raise UnsupportedStructure("this is not a chain of acid centres")
    centers, links = scanned
    name = _chain_name(mol, centers, links, priority)
    if name is None:
        name = _dithioxane_name(mol, centers, links) or _substitutive_name(mol, centers, links)
    if name is None:
        raise UnsupportedStructure("this chain of acid centres matches no retained polyacid name yet")
    return name
