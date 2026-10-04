"""Assembly of a complete natural-product name from a stereoparent candidate (P-101.7)."""

from collections import Counter

from rdkit import Chem
from dataclasses import dataclass

from ._common import UnsupportedStructure, adjacency, halogen_substituents, multiplied_word
from ._np_config import configuration
from ._np_core import FACES, loc_key
from ._np_name import (
    _ACYL_CLASSES,
    _ACYL_SUFFIX,
    _PREFIX,
    _SENIORITY,
    _SUFFIX,
    Loc,
    _acid_anion,
    alkyl_count,
    classify,
    final_labels,
)
from ._np_rings import bridge_prefixes, components
from ._np_text import core_text, locant_pair, op_prefixes, unsaturation
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch


_MAX_UNSATURATION_CHANGES = 6
_MAX_HYDRO_PAIRS = 3
_MAX_TOTAL_COST = 3
_MAX_SKELETAL_MODIFICATIONS = 2


@dataclass
class Built:
    name: str
    cost: int
    key: tuple


def _components(view, mapped):
    seen, comps = set(mapped), []
    for start in view.adj:
        if start in seen:
            continue
        comp, stack = {start}, [start]
        seen.add(start)
        while stack:
            for n in view.adj[stack.pop()]:
                if n not in seen:
                    seen.add(n)
                    comp.add(n)
                    stack.append(n)
        comps.append(comp)
    return comps


def _check_single_attachment(view, mapped):
    for comp in _components(view, mapped):
        links = [(a, n) for n in comp for a in view.adj[n] if a in mapped]
        if len(links) > 1:
            raise UnsupportedStructure("a ring or bridge added to a natural-product parent is not supported")


def nondetachable_cost(cand):
    return len(cand.skel.ops) + len(cand.cyclo) + len(cand.replaced)


def build(cand, view):
    skel, mapping, parent = cand.skel, cand.mapping, cand.parent
    mapped = set(mapping.values())
    if any(op[0] == "seco" and frozenset(op[1:]) in {frozenset(c) for c in cand.cyclo} for op in skel.ops):
        raise UnsupportedStructure("a cleaved bond that is formed again")
    ring_comps = components(view, mapped)
    bridge_atoms = set().union(*(c.atoms for c in ring_comps)) if ring_comps else set()
    groups = classify(cand, view)
    groups.branches = [(loc, root) for loc, root in groups.branches if root not in bridge_atoms]
    alkyls = alkyl_count(cand, view, groups)
    if alkyls is None:
        raise UnsupportedStructure("an alkyl group on an acyclic part of the parent")
    cost = nondetachable_cost(cand) + alkyls + len(ring_comps)

    classes = groups.classes
    principal = next((c for c in _SENIORITY if c in classes), None)
    if "ester" in classes and "ester_o" in classes:
        raise UnsupportedStructure("esters of both an acid and an alcohol of the parent are not supported")
    branches = list(groups.branches)
    for cls in _ACYL_CLASSES:
        if cls in classes and cls != principal:
            for loc, _, extra in classes.pop(cls):
                if "root" not in extra:
                    raise UnsupportedStructure("a terminal acyl group below the principal group is not supported")
                branches.append((loc, extra["root"]))

    enes_bonds, hydro, dehydro = unsaturation(cand, view)
    if len(hydro) % 2 or len(dehydro) % 2:
        raise UnsupportedStructure("indicated hydrogen would be needed")
    if len(enes_bonds) + len(hydro) // 2 + len(dehydro) // 2 > _MAX_UNSATURATION_CHANGES or (len(hydro) + len(dehydro)) // 2 > _MAX_HYDRO_PAIRS and not parent.name.endswith("carotene"):
        raise UnsupportedStructure("too many changes of the degree of hydrogenation for this parent")
    if _pair_count(hydro, dehydro) + cost > _MAX_TOTAL_COST and not parent.name.endswith("carotene"):
        raise UnsupportedStructure("too many modifications for this parent")
    cost += _pair_count(hydro, dehydro)
    final = _final(skel)
    hydro = sorted((final(x) for x in hydro), key=loc_key)
    dehydro = sorted((final(x) for x in dehydro), key=loc_key)
    enes = [locant_pair(skel, final, a, b) for a, b, o in enes_bonds if o == 2]
    ynes = [locant_pair(skel, final, a, b) for a, b, o in enes_bonds if o == 3]
    enes.sort(key=lambda t: loc_key(t.split("(")[0]))
    ynes.sort(key=lambda t: loc_key(t.split("(")[0]))

    config = configuration(cand, view)
    graph = adjacency(view.mol)
    halogens = halogen_substituents(view.mol)

    prefix_groups = {}
    for cls, members in classes.items():
        if cls == principal:
            continue
        for loc, atoms, extra in members:
            atom = next(iter(sorted(atoms)))
            prefix_groups.setdefault(_PREFIX[cls], {"locants": [], "compound": False})["locants"].append(
                Loc(final(loc), config.faces.get(atom, ""))
            )
    named = [
        (loc, root, *name_branch(graph, root, mapping[loc], halogens, frozenset(), mol=view.mol, unsaturated=True))
        for loc, root in branches
    ]
    repeats = Counter((loc, name) for loc, _, name, _ in named)
    for loc, root, name, compound in named:
        face = "" if repeats[(loc, name)] > 1 else config.faces.get(root, "")
        prefix_groups.setdefault(name, {"locants": [], "compound": compound})["locants"].append(Loc(final(loc), face))
    for entry in prefix_groups.values():
        entry["locants"].sort(key=lambda t: loc_key(t.base))
    prefix = format_substituent_prefixes(prefix_groups) if prefix_groups else ""

    suffix, locants, alkyl_word, anion = "", [], "", ""
    if principal is not None:
        members = sorted(classes[principal], key=lambda m: loc_key(m[0]))
        if principal in ("acid", "ester", "amide"):
            if len({extra["kind"] for _, _, extra in members}) != 1:
                raise UnsupportedStructure("a mix of ring and chain acyl groups is not supported")
            word = _ACYL_SUFFIX[(principal, members[0][2]["kind"])]
        else:
            word = _SUFFIX[principal]
        suffix = multiplied_word(len(members), word)
        locants = [str(Loc(final(loc), "" if extra.get("kind") == "o" else config.faces.get(extra.get("anchor", next(iter(sorted(atoms)))), ""))) for loc, atoms, extra in members]
        if principal == "ester":
            alkyls_named = {}
            for _, _, extra in members:
                oxygen, alkyl = extra["ester"]
                name, compound = name_branch(graph, alkyl, oxygen, halogens, frozenset(), mol=view.mol, unsaturated=True)
                alkyls_named[name] = compound
            if len(alkyls_named) != 1:
                raise UnsupportedStructure("esters with different alkyl groups are not supported")
            ((name, compound),) = alkyls_named.items()
            alkyl_word = name
            if len(members) > 1:
                alkyl_word = multiplying_prefix(len(members), compound=compound) + (f"({name})" if compound else name)
        if principal == "ester_o":
            anions = {_acid_anion(view, extra["anchor"], extra["acyl"], mapped) for _, _, extra in members}
            if len(anions) != 1:
                raise UnsupportedStructure("O-acyl groups from different acids are not supported")
            (anion,) = anions
            if len(members) > 1:
                if not anion.isalpha():
                    raise UnsupportedStructure("a substituted acid part repeated on a natural product is not supported")
                anion = multiplying_prefix(len(members), compound=False) + anion

    stem_core = core_text(parent.name, enes, ynes, suffix, locants)
    descriptor = ",".join(text for _, text in sorted(config.parent)) + "-" if config.parent else ""
    hydro_text = _hydro_text(dehydro, "dehydro") + _hydro_text(hydro, "hydro")
    cyclo_text = _cyclo_text(cand, final, config, view)
    ops = op_prefixes(cand, final, cyclo_text)
    replacement = _replacement_text(cand, final, view)
    if replacement:
        ops = [replacement] + ops
    bridges = bridge_prefixes(ring_comps, cand, view, config, final)
    if nondetachable_cost(cand) + len(ring_comps) > _MAX_SKELETAL_MODIFICATIONS:
        raise UnsupportedStructure("too many skeletal modifications")
    ops = bridges + ops
    nondetachable = "-".join(ops) + ("-" if ops and descriptor else "")
    side = f"({','.join(text for _, text in sorted(config.side))})-" if config.side else ""
    detachable = f"{prefix}-" if prefix and (hydro_text or nondetachable or descriptor) else prefix
    if hydro_text and prefix:
        detachable = f"{prefix}-"
    body = f"{side}{detachable}{hydro_text}{nondetachable}{descriptor}{stem_core}"
    name = " ".join(part for part in (alkyl_word, body, anion) if part)
    nondet = nondetachable_cost(cand)
    rearranged = any(op[0] == "seco" for op in cand.skel.ops) or bool(cand.cyclo)
    removed = sorted((loc_key(op[1])[1] for op in cand.skel.ops if op[0] == "nor"), reverse=True)
    return Built(name, cost, (cost, nondet > 0, rearranged, -len(mapping), tuple(-n for n in removed)))


def _final(skel):
    labels = final_labels(skel)
    return lambda loc: labels.get(loc, loc)


def _hydro_text(locants, word, reverse=False):
    if not locants:
        return ""
    prefix = multiplying_prefix(len(locants), compound=False) or ""
    return f"{','.join(locants)}-{prefix}{word}" + "-"


_A_PREFIX = {
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora", "Al": "alumina",
    "Ga": "galla", "In": "inda", "Tl": "thalla", "C": "carba",
}
_A_ORDER = list(_A_PREFIX)


def _replacement_text(cand, final, view):
    groups = {}
    for loc in cand.replaced:
        groups.setdefault(view.elem[cand.mapping[loc]], []).append(final(loc))
    parts = []
    for element in sorted(groups, key=_A_ORDER.index):
        locs = sorted(groups[element], key=loc_key)
        word = _A_PREFIX[element]
        parts.append(f"{','.join(locs)}-{multiplying_prefix(len(locs), compound=False) if len(locs) > 1 else ''}{word}")
    return "-".join(parts)


def _cyclo_text(cand, final, config, view):
    texts = []
    parent = cand.parent
    for a, b in sorted(cand.cyclo, key=lambda ab: tuple(sorted((loc_key(final(ab[0])), loc_key(final(ab[1])))))):
        pair = sorted((a, b), key=lambda x: loc_key(final(x)))
        cited = []
        for x, y in ((pair[0], pair[1]), (pair[1], pair[0])):
            stereo = view.mol.GetAtomWithIdx(cand.mapping[x]).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            face = config.faces.get(cand.mapping[y], "") if x not in parent.centers and stereo else ""
            cited.append(f"{final(x)}{FACES.get(face, '')}")
        texts.append(",".join(cited))
    return texts


def _pair_count(hydro, dehydro):
    return (len(hydro) + len(dehydro)) // 2
