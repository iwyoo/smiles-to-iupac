"""Mononuclear oxoacids of P, As, Sb, S, Se and Te modified by infixes (P-67.1.2.3), and the halides, pseudohalides,
amides and hydrazides they form (P-67.1.2.5, P-67.1.2.6): the replacing groups are cited as infixes while an acidic
hydrogen remains, else the senior replaced class becomes the class name, N- and P-substituents take letter locants."""

from itertools import permutations

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, group_substituents, halogen_substituents
from ._numerals import numerical_term
from ._multiplicative_text import enclose
from ._phosphate import format_ester_words
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_CHALCOGEN = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_CHALCOGEN_PREFIX = {16: "thio", 34: "seleno", 52: "telluro"}
_ORDER = {8: 0, 16: 1, 34: 2, 52: 3}
_PENTAVALENT_STEMS = {15: ("phosphin", "phosphon", "phosphor"), 33: ("arsin", "arson", "arsor"), 51: ("stibin", "stibon", "stibor")}
_CHALCOGEN_STEMS = {16: "sulfur", 34: "selen", 52: "tellur"}
_LOCANT = {15: "P", 33: "As", 51: "Sb"}
_HALIDE = {9: ("fluorid", "fluoride"), 17: ("chlorid", "chloride"), 35: ("bromid", "bromide"), 53: ("iodid", "iodide")}
_PSEUDOHALIDES = (
    ("[O;D2]-C#N", "cyanatid", "cyanate"),
    ("[N;D2]=C=O", "isocyanatid", "isocyanate"),
    ("[S;D2]-C#N", "thiocyanatid", "thiocyanate"),
    ("[N;D2]=C=S", "isothiocyanatid", "isothiocyanate"),
    ("[Se;D2]-C#N", "selenocyanatid", "selenocyanate"),
    ("[Te;D2]-C#N", "tellurocyanatid", "tellurocyanate"),
    ("[N;D2]=C=[Se]", "isoselenocyanatid", "isoselenocyanate"),
    ("[C;D2]#N", "cyanid", "cyanide"),
    ("[N;D2]#[C;D1]", "isocyanid", "isocyanide"),
    ("[N;D2]=[N+]=[N-]", "azid", "azide"),
)
_PSEUDOHALIDE_PATTERNS = [(Chem.MolFromSmarts(s), infix, term) for s, infix, term in _PSEUDOHALIDES]
_CLASS_ORDER = ["bromide", "chloride", "fluoride", "iodide", "azide", "cyanide", "isocyanide", "isocyanate", "cyanate",
                "thiocyanate", "isothiocyanate", "selenocyanate", "isoselenocyanate", "tellurocyanate"]
_PARENTHESIZED = ("thiocyanatid", "selenocyanatid", "tellurocyanatid", "peroxo")
_ACYL_NAMES = {15: "phosphoryl", 16: "sulfuryl", 34: "selenonyl", 52: "telluronyl"}
_VOWELS = "aeiou"
_HALOGEN_ACIDS = {
    9: ("hypofluorous", "fluorous", "fluoric", "perfluoric"),
    17: ("hypochlorous", "chlorous", "chloric", "perchloric"),
    35: ("hypobromous", "bromous", "bromic", "perbromic"),
    53: ("hypoiodous", "iodous", "iodic", "periodic"),
}
_HALOGEN_PREFIX = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_PSEUDOHALIDE_PREFIX = {
    "cyanatid": "cyanato", "isocyanatid": "isocyanato", "thiocyanatid": "thiocyanato", "isothiocyanatid": "isothiocyanato",
    "selenocyanatid": "selenocyanato", "isoselenocyanatid": "isoselenocyanato", "tellurocyanatid": "tellurocyanato",
    "cyanid": "cyano", "isocyanid": "isocyano", "azid": "azido",
}


def _plain_tree(mol, graph, root, behind):
    """Whether the group starting at carbon `root` holds only carbon, halogen and sp3 oxygen atoms, so that nothing in
    it outranks the acid derivative being named."""
    seen, stack = {behind}, [root]
    while stack:
        index = stack.pop()
        if index in seen:
            continue
        seen.add(index)
        atom = mol.GetAtomWithIdx(index)
        z = atom.GetAtomicNum()
        if atom.GetIsotope() or atom.GetFormalCharge():
            return False
        if z == 8:
            if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
                return False
        elif z != 6 and z not in _HALIDE:
            return False
        if z == 6 and any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8 for b in atom.GetBonds()
        ):
            return False
        stack.extend(graph[index])
    return True


def _carbon_substituents(mol, graph, halogens, atom, skip):
    """[(name, compound)] of the carbon groups on `atom` outside `skip`, or None when another atom is bonded."""
    found = []
    for n in graph[atom.GetIdx()]:
        if n in skip:
            continue
        neighbor = mol.GetAtomWithIdx(n)
        if neighbor.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(atom.GetIdx(), n).GetBondTypeAsDouble() != 1.0:
            return None
        if not _plain_tree(mol, graph, n, atom.GetIdx()):
            return None
        found.append(name_branch(graph, n, atom.GetIdx(), halogens, mol=mol))
    return found


def _peroxo_name(first, last):
    elements = sorted({first, last} - {8}, key=lambda z: _CHALCOGEN_PREFIX[z])
    if not elements:
        return "peroxo"
    if len(elements) == 1 and first == last:
        return "di" + _CHALCOGEN_PREFIX[elements[0]] + "peroxo"
    return "".join(_CHALCOGEN_PREFIX[z] for z in elements) + "peroxo"


def _classify(mol, graph, halogens, centre, neighbor):
    """The ligand description of `neighbor` bonded to the acid centre, or None when it is not an acid ligand."""
    order = mol.GetBondBetweenAtoms(centre.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble()
    z = neighbor.GetAtomicNum()
    if neighbor.GetIsotope():
        return None
    if order == 1.0:
        for pattern, infix, term in _PSEUDOHALIDE_PATTERNS:
            for match in mol.GetSubstructMatches(pattern):
                if match[0] == neighbor.GetIdx():
                    return {"role": "single", "kind": "pseudohalide", "infix": infix, "term": term, "H": False}
    if z in _CHALCOGEN and not neighbor.GetFormalCharge():
        if neighbor.GetDegree() == 1:
            if order == 2.0 and not neighbor.GetTotalNumHs():
                return {"role": "ylidene", "kind": "chalcogen", "z": z, "H": False}
            if order == 1.0 and neighbor.GetTotalNumHs() == 1:
                return {"role": "single", "kind": "chalcogen", "z": z, "H": True}
        elif neighbor.GetDegree() == 2 and order == 1.0 and not neighbor.GetTotalNumHs():
            other = next(n for n in neighbor.GetNeighbors() if n.GetIdx() != centre.GetIdx())
            if other.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(neighbor.GetIdx(), other.GetIdx()).GetBondTypeAsDouble() == 1.0:
                if not _plain_tree(mol, graph, other.GetIdx(), neighbor.GetIdx()):
                    return None
                return {"role": "single", "kind": "ester", "z": z, "H": False, "R": other.GetIdx()}
            if (
                other.GetAtomicNum() in _CHALCOGEN
                and other.GetDegree() == 1
                and other.GetTotalNumHs() == 1
                and not other.GetFormalCharge()
                and mol.GetBondBetweenAtoms(neighbor.GetIdx(), other.GetIdx()).GetBondTypeAsDouble() == 1.0
            ):
                return {"role": "single", "kind": "peroxo", "infix": _peroxo_name(z, other.GetAtomicNum()), "H": True, "z": z, "z2": other.GetAtomicNum()}
        return None
    if neighbor.GetFormalCharge():
        return None
    if z in _HALIDE and order == 1.0 and neighbor.GetDegree() == 1:
        infix, term = _HALIDE[z]
        return {"role": "single", "kind": "halide", "infix": infix, "term": term, "H": False}
    if z != 7 or neighbor.IsInRing():
        return None
    if order == 3.0 and neighbor.GetDegree() == 1 and not neighbor.GetTotalNumHs():
        return {"role": "ylidene", "kind": "nitride", "H": False}
    others = [n for n in neighbor.GetNeighbors() if n.GetIdx() != centre.GetIdx()]
    nitrogens = [n for n in others if n.GetAtomicNum() == 7]
    if len(nitrogens) > 1:
        return None
    if nitrogens:
        far = nitrogens[0]
        if far.GetFormalCharge() or mol.GetBondBetweenAtoms(neighbor.GetIdx(), far.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
        near = _carbon_substituents(mol, graph, halogens, neighbor, {centre.GetIdx(), far.GetIdx()})
        beyond = _carbon_substituents(mol, graph, halogens, far, {neighbor.GetIdx()})
        if near is None or beyond is None or len(beyond) + far.GetTotalNumHs() != 2:
            return None
        if len(near) + neighbor.GetTotalNumHs() != (0 if order == 2.0 else 1):
            return None
        hydrazone = order == 2.0
        return {"role": "ylidene" if hydrazone else "single", "kind": "hydrazone" if hydrazone else "hydrazide",
                "infix": "hydrazon" if hydrazone else "hydrazid", "term": "hydrazide", "H": False, "subs": [near, beyond]}
    subs = _carbon_substituents(mol, graph, halogens, neighbor, {centre.GetIdx()})
    if subs is None:
        return None
    if order == 2.0 and len(subs) + neighbor.GetTotalNumHs() == 1:
        return {"role": "ylidene", "kind": "imide", "infix": "imid", "H": False, "subs": [subs]}
    if order == 1.0 and len(subs) + neighbor.GetTotalNumHs() == 2:
        return {"role": "single", "kind": "amide", "infix": "amid", "term": "amide", "H": False, "subs": [subs]}
    return None


def _parse(mol, centre):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    if centre.GetFormalCharge() or centre.GetIsotope() or centre.IsInRing() or centre.GetNumRadicalElectrons():
        return None
    z = centre.GetAtomicNum()
    ylidenes, singles, organyl = [], [], []
    hydrogens = centre.GetTotalNumHs()
    for neighbor in centre.GetNeighbors():
        order = mol.GetBondBetweenAtoms(centre.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble()
        if neighbor.GetAtomicNum() == 6 and order == 1.0 and not _plain_tree(mol, graph, neighbor.GetIdx(), centre.GetIdx()):
            return None
        if neighbor.GetAtomicNum() == 6 and order == 1.0 and not any(
            p for p, _, _ in _PSEUDOHALIDE_PATTERNS if any(m[0] == neighbor.GetIdx() for m in mol.GetSubstructMatches(p))
        ):
            organyl.append(neighbor)
            continue
        found = _classify(mol, graph, halogens, centre, neighbor)
        if found is None:
            return None
        found["atom"] = neighbor.GetIdx()
        (ylidenes if found["role"] == "ylidene" else singles).append(found)
    k = len(organyl) + hydrogens
    nitrido = sum(1 for y in ylidenes if y["kind"] == "nitride")
    if z in _PENTAVALENT_STEMS:
        valence_v = len(ylidenes) == 1
        if len(ylidenes) > 1 or len(singles) != 3 - k - nitrido:
            return None
    else:
        if k or nitrido or len(singles) != 2 or len(ylidenes) not in (1, 2):
            return None
        valence_v = len(ylidenes) == 2
    modified = any(
        (s["kind"] not in ("chalcogen", "ester")) or s["z"] != 8 for s in singles if s["kind"] != "ester" or s["z"] != 8
    ) or any(y["kind"] != "chalcogen" or y["z"] != 8 for y in ylidenes)
    if not modified or not singles:
        return None
    esters = [x for x in singles if x["kind"] == "ester"]
    if esters and any(x["kind"] in ("hydrazide",) for x in singles):
        return None
    if z in _CHALCOGEN_STEMS and any(p["kind"] in ("hydrazide", "hydrazone") for p in (*ylidenes, *singles)):
        return None
    if z in _PENTAVALENT_STEMS and not organyl and not ylidenes and not any(p["kind"] in ("chalcogen", "peroxo", "ester") for p in singles):
        return None
    return {"centre": centre, "z": z, "ylidenes": ylidenes, "singles": singles, "organyl": organyl, "k": k, "pentavalent": valence_v, "graph": graph, "halogens": halogens}


def _parse_prefixed(mol, centre):
    """Silicic, nitrous and halogen oxoacids, whose replacements are cited as prefixes (P-67.1.2.2)."""
    if centre.GetFormalCharge() or centre.GetIsotope() or centre.IsInRing() or centre.GetTotalNumHs() or centre.GetNumRadicalElectrons():
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    z = centre.GetAtomicNum()
    ylidenes, singles = [], []
    for neighbor in centre.GetNeighbors():
        found = _classify(mol, graph, halogens, centre, neighbor)
        if found is None or found["kind"] not in ("chalcogen", "halide", "pseudohalide", "ester"):
            return None
        if found["kind"] == "ester" and not _plain_tree(mol, graph, found["R"], neighbor.GetIdx()):
            return None
        found["atom"] = neighbor.GetIdx()
        (ylidenes if found["role"] == "ylidene" else singles).append(found)
    if z == 14:
        ok = not ylidenes and len(singles) == 4
    elif z == 7:
        ok = len(ylidenes) == 1 and len(singles) == 1 and singles[0]["kind"] in ("chalcogen", "ester")
    elif z in _HALOGEN_ACIDS:
        ok = len(ylidenes) <= 3 and len(singles) == 1 and singles[0]["kind"] in ("chalcogen", "ester")
    else:
        return None
    if not ok or not any(s["H"] or s["kind"] == "ester" for s in singles):
        return None
    if not any(p["kind"] not in ("chalcogen", "ester") or p["z"] != 8 for p in (*ylidenes, *singles)):
        return None
    return {"centre": centre, "z": z, "ylidenes": ylidenes, "singles": singles, "family": "prefix", "graph": graph, "halogens": halogens, "pentavalent": True}


def _acid_parts(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    parts = []
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z in _LOCANT or z in _CHALCOGEN_STEMS:
            found = _parse(mol, atom)
        elif z == 14 or z == 7 or z in _HALOGEN_ACIDS:
            found = _parse_prefixed(mol, atom)
        else:
            continue
        if found is not None:
            parts.append(found)
    return parts[0] if len(parts) == 1 else None


def has_noncarbon_oxoacid_shape(mol) -> bool:
    try:
        return _acid_parts(mol) is not None
    except UnsupportedStructure:
        return False


def _assign_letters(groups):
    """Locants N, N', N'' ... for the amide and imide nitrogens (P-67.1.2.6.1): lowest set first, then the
    alphabetically first substituent."""
    by_kind = {}
    for group in groups:
        by_kind.setdefault(group["kind"], []).append(group)
    kinds = sorted(by_kind, key=lambda kind: {"amide": "amid", "imide": "imid"}[kind])
    best = None
    ordered_kinds = [by_kind[kind] for kind in kinds]
    for choice in _product_permutations(ordered_kinds):
        flat = [g for part in choice for g in part]
        levels = []
        for level, group in enumerate(flat):
            for name, _ in group["subs"][0]:
                levels.append((name, level))
        key = (sorted(level for _, level in levels), [level for _, level in sorted(levels)])
        if best is None or key < best[0]:
            best = (key, flat)
    return best[1]


def _product_permutations(parts):
    if not parts:
        yield []
        return
    for first in permutations(parts[0]):
        for rest in _product_permutations(parts[1:]):
            yield [list(first)] + rest


def _complete_substitution_prefix(groups):
    """'hexamethyl' when every substitutable nitrogen position carries the same group, so locants are unnecessary."""
    capacity = {"amide": 2, "imide": 1, "hydrazide": 3, "hydrazone": 2}
    items = [item for g in groups for part in g["subs"] for item in part]
    if len({name for name, _ in items}) != 1 or len(items) != sum(capacity[g["kind"]] for g in groups):
        return None
    return format_mononuclear_prefixes(items)


def _prime(count):
    return "'" * count


def _substituent_prefix(parts, groups):
    """Prefix text for the carbon groups on the acid centre and on the nitrogens of amide, imide and hydrazide groups."""
    centre_entries = parts["centre_entries"]
    if not groups:
        return format_mononuclear_prefixes(centre_entries)
    entries = {}
    letter_groups = [g for g in groups if g["kind"] in ("amide", "imide")]
    for level, group in enumerate(_assign_letters(letter_groups) if letter_groups else []):
        for item in group["subs"][0]:
            entries.setdefault("N" + _prime(level), []).append(item)
    hydrazines = [g for g in groups if g["kind"] in ("hydrazide", "hydrazone")]
    for level, group in enumerate(hydrazines):
        near, far = group["subs"]
        for item in near:
            entries.setdefault("1" + _prime(level), []).append(item)
        for item in far:
            entries.setdefault("2" + _prime(level), []).append(item)
    for item in centre_entries:
        entries.setdefault(_LOCANT.get(parts["z"], "S"), []).append(item)
    grouped = group_substituents(entries)
    return format_substituent_prefixes(grouped) if grouped else ""


def _compose(stem, tokens, ending):
    """Stem, infixes and ending with the euphonic 'o' of P-67.1.2.3.5; `tokens` is [(text, parenthesized)]."""
    parts = [{"text": stem, "paren": False, "trail": ""}]
    for text, paren in tokens:
        previous = parts[-1]
        if text[0] not in _VOWELS and not previous["text"].endswith("o"):
            previous["trail"] += "o"
        parts.append({"text": text, "paren": paren, "trail": ""})
    last = parts[-1]
    if ending == "ous" and last["text"].endswith("o") and len(parts) > 1:
        last["text"] = last["text"][:-1]
    last["trail"] += ending
    return "".join(("(" if p["paren"] else "") + p["text"] + p["trail"] + (")" if p["paren"] else "") for p in parts)


def _infix_tokens(positions):
    counts = {}
    for text in positions:
        counts[text] = counts.get(text, 0) + 1
    tokens = []
    for text in sorted(counts, key=lambda t: t.lstrip("di") if t.startswith("dithio") or t.startswith("diseleno") else t):
        paren = any(marker in text for marker in _PARENTHESIZED)
        count = counts[text]
        if count > 1:
            tokens.append((numerical_term(count) + text, paren))
        else:
            tokens.append((text, paren))
    return tokens


def _position_infixes(parts, exclude_class):
    texts = []
    for p in (*parts["ylidenes"], *parts["singles"]):
        if p["kind"] in ("chalcogen", "ester"):
            if p["z"] != 8:
                texts.append(_CHALCOGEN_PREFIX[p["z"]])
        elif p["kind"] == "nitride":
            texts.append("nitrid")
        elif p in exclude_class:
            continue
        else:
            texts.append(p["infix"])
    return texts


def _acid_locants(parts):
    """The element symbols of the hydrogen-bearing atoms that tell tautomers apart (P-67.1.2.4.1), '' when all the
    chalcogen positions are of one element."""
    positions = [p for p in (*parts["ylidenes"], *parts["singles"]) if p["kind"] in ("chalcogen", "peroxo", "ester")]
    elements = {z for p in positions for z in (p["z"], p.get("z2", p["z"]))}
    if not any(y["kind"] == "chalcogen" for y in parts["ylidenes"]) or len(elements) < 2:
        return ""
    replaced = [p for p in positions if p["kind"] == "peroxo" or p["z"] != 8]
    with_hydrogen = [p for p in replaced if p["H"]]
    if with_hydrogen:
        cited = [_CHALCOGEN[p["z"]] + (_CHALCOGEN[p["z2"]] if p["kind"] == "peroxo" else "") for p in with_hydrogen]
    else:
        cited = [_CHALCOGEN[p["z"]] for p in positions if p["H"]]
    return ",".join(sorted(cited, key=lambda c: _ORDER[[z for z, e in _CHALCOGEN.items() if e == c[0]][0]]))


def name_noncarbon_oxoacid(mol) -> str:
    parts = _acid_parts(mol)
    if parts is None:
        raise UnsupportedStructure("this is not a mononuclear noncarbon oxoacid modified by functional replacement")
    if parts.get("family") == "prefix":
        return _name_prefixed(mol, parts)
    graph, halogens, centre = parts["graph"], parts["halogens"], parts["centre"]
    parts["centre_entries"] = [name_branch(graph, c.GetIdx(), centre.GetIdx(), halogens, mol=mol) for c in parts["organyl"]]
    z = parts["z"]
    singles, ylidenes = parts["singles"], parts["ylidenes"]
    acidic = any(s["H"] for s in singles)
    esters = [x for x in singles if x["kind"] == "ester"]

    class_members, class_terms = [], []
    if not acidic and not esters:
        for group in (("halide",), ("pseudohalide",), ("amide",), ("hydrazide",)):
            members = [s for s in singles if s["kind"] in group]
            if members:
                class_members = members
                break
        if not class_members:
            raise UnsupportedStructure("an acid centre with no replaceable class is not named here")
        class_terms = sorted((s["term"] for s in class_members), key=lambda t: (_CLASS_ORDER.index(t) if t in _CLASS_ORDER else 99, t))

    groups = [p for p in (*ylidenes, *singles) if p["kind"] in ("amide", "imide", "hydrazide", "hydrazone")]
    n_substituted = any(any(sub for sub in p["subs"]) for p in groups)
    if class_members:
        acyl = _acyl_class_name(parts, class_members, class_terms)
        if acyl is not None:
            amide = [g for g in groups if g["kind"] == "amide"]
            prefix = _substituent_prefix(dict(parts, centre_entries=[]), amide) if n_substituted and amide else ""
            return prefix + acyl
    prefix = _substituent_prefix(parts, groups if n_substituted else [])
    only_amide_modifies = acidic and z in _PENTAVALENT_STEMS and all(p["kind"] == "amide" or (p["kind"] == "chalcogen" and p["z"] == 8) for p in (*ylidenes, *singles))
    if n_substituted and not parts["centre_entries"] and (len(groups) > 1 or only_amide_modifies):
        prefix = _complete_substitution_prefix(groups) or prefix
    infixes = _infix_tokens(_position_infixes(parts, class_members))
    pentavalent = parts["pentavalent"]
    if z in _PENTAVALENT_STEMS:
        stem = _PENTAVALENT_STEMS[z][len(singles) + (1 if any(y["kind"] == "nitride" for y in ylidenes) else 0) - 1]
    else:
        stem = _CHALCOGEN_STEMS[z]
        if z == 16 and pentavalent and infixes and infixes[0][0] == "amid" and all(t in ("thio", "seleno", "telluro") for t, _ in infixes[1:]):
            stem, infixes = "sulfam", infixes[1:]
    if esters:
        return _ester_name(mol, parts, esters, prefix, stem, infixes, acidic)
    ending = "ic" if pentavalent else "ous"
    body = _compose(stem, infixes, ending)
    if class_terms:
        return prefix + body + " " + _class_phrase(class_terms)
    locants = _acid_locants(parts)
    return prefix + (body + " " + locants + "-acid" if locants else body + " acid")


def _class_phrase(terms):
    counts = {}
    for term in terms:
        counts[term] = counts.get(term, 0) + 1
    return " ".join((numerical_term(c) if c > 1 else "") + t for t, c in sorted(counts.items(), key=lambda it: (_CLASS_ORDER.index(it[0]) if it[0] in _CLASS_ORDER else 99, it[0])))


def _acyl_class_name(parts, members, terms):
    """'phosphoryl trichloride' and its kin: identical halides or pseudohalides of phosphoric, sulfuric, selenic and
    telluric acid (P-67.1.2.5.1)."""
    z = parts["z"]
    if z not in _ACYL_NAMES or parts["k"] or not parts["pentavalent"]:
        return None
    if any(y["kind"] != "chalcogen" or y["z"] != 8 for y in parts["ylidenes"]):
        return None
    if members[0]["kind"] not in ("halide", "pseudohalide") or len({m.get("infix") for m in members}) != 1:
        return None
    others = [s for s in parts["singles"] if s not in members]
    if not others:
        return _ACYL_NAMES[z] + " " + _class_phrase(terms)
    if z == 16 and len(others) == 1 and others[0]["kind"] == "amide" and len(members) == 1:
        return "sulfamoyl " + terms[0]
    return None


def _name_prefixed(mol, parts):
    z, ylidenes, singles = parts["z"], parts["ylidenes"], parts["singles"]
    counts = {}
    for p in (*ylidenes, *singles):
        if p["kind"] in ("chalcogen", "ester"):
            text = _CHALCOGEN_PREFIX.get(p["z"])
        elif p["kind"] == "halide":
            text = next(v for k, v in _HALOGEN_PREFIX.items() if _HALIDE[k][0] == p["infix"])
        else:
            text = _PSEUDOHALIDE_PREFIX[p["infix"]]
        if text:
            counts[text] = counts.get(text, 0) + 1
    prefix = "".join((numerical_term(c) if c > 1 else "") + t for t, c in sorted(counts.items()))
    if z == 14:
        base = "silicic"
    elif z == 7:
        base = "nitrous"
    else:
        base = _HALOGEN_ACIDS[z][len(ylidenes)]
    esters = [x for x in singles if x["kind"] == "ester"]
    if esters:
        return _prefixed_ester_name(mol, parts, esters, prefix, base)
    locants = _acid_locants(parts) if z != 14 else ""
    return prefix + base + (" " + locants + "-acid" if locants else " acid")


def _prefixed_ester_name(mol, parts, esters, prefix, base):
    words = _ester_words(mol, parts, esters, _acid_locants_present(parts))
    anion = base[:-2] + "ate" if base.endswith("ic") else base[:-3] + "ite"
    return " ".join(words) + " " + prefix + anion


def _ester_words(mol, parts, esters, cited):
    graph, halogens = parts["graph"], parts["halogens"]
    by_name = {}
    for ester in esters:
        name, _ = name_branch(graph, ester["R"], ester["atom"], halogens, mol=mol)
        by_name.setdefault(name, []).append(_CHALCOGEN[ester["z"]])
    words = []
    for name in sorted(by_name, key=alpha_sort_key):
        letters = sorted(by_name[name], key=lambda e: _ORDER[[z for z, sym in _CHALCOGEN.items() if sym == e][0]])
        if cited and len(letters) == 1 and (name[0].isdigit() or name[0] in "([{"):
            word = enclose(name)
        else:
            word = format_ester_words([name] * len(letters))
        words.append((",".join(letters) + "-" if cited else "") + word)
    hydrogens = sum(1 for s in parts["singles"] if s["H"])
    if hydrogens:
        words.append((numerical_term(hydrogens) if hydrogens > 1 else "") + "hydrogen")
    return words


def _ester_name(mol, parts, esters, prefix, stem, infixes, acidic):
    """Ester words, 'hydrogen' words for the acidic positions left, then the anion name (P-67.1.3.2)."""
    words = _ester_words(mol, parts, esters, _acid_locants_present(parts))
    anion = _compose(stem, infixes, "ate" if parts["pentavalent"] else "ite")
    return " ".join(words) + " " + prefix + anion


def _acid_locants_present(parts):
    positions = [p for p in (*parts["ylidenes"], *parts["singles"]) if p["kind"] in ("chalcogen", "peroxo", "ester")]
    elements = {z for p in positions for z in (p["z"], p.get("z2", p["z"]))}
    return len(elements) >= 2
