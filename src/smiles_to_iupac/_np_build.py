"""Assembly of a complete natural-product name from a stereoparent candidate (P-101.7)."""

import re
from collections import Counter
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, halogen_substituents, multiplied_word
from ._np_chain import chain_template
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
    principal_groups_left_outside,
)
from ._np_fusion import assign_primes, fusion_text, hydro_texts, indicated_texts, name_fused
from ._np_rings import bridge_prefixes, components, name_spiro, split_components
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
    implied_total: int = 0
    implied_cited: int = 0
    first_face: str = ""


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
    bridge_comps, fused_comps, spiro_comps = split_components(ring_comps, cand, view)
    if len(spiro_comps) > 1 or (spiro_comps and (bridge_comps or fused_comps)):
        raise UnsupportedStructure("a spiro ring together with other added rings is not supported")
    spiro = name_spiro(spiro_comps[0], cand, view) if spiro_comps else None
    parent_prime = ring_prime = ""
    if spiro:
        ring_key = spiro.ring_name.lstrip("0123456789,[]-")
        if ring_key < parent.name:
            parent_prime = "′"
        else:
            ring_prime = "′"
    groups = classify(cand, view, bridge_atoms)
    groups.branches = [(loc, root) for loc, root in groups.branches if root not in bridge_atoms]
    alkyls = alkyl_count(cand, view, groups)
    if alkyls is None:
        raise UnsupportedStructure("an alkyl group on an acyclic part of the parent")
    cost = nondetachable_cost(cand) + alkyls + len(ring_comps)

    classes = groups.classes
    principal = next((c for c in _SENIORITY if c in classes), None)
    chain = None
    if principal is None:
        chain = chain_template(cand, view, groups.branches)
        if chain is not None:
            groups.branches = [b for b in groups.branches if b[1] != chain[1]]
            classes["yl"] = [(chain[0], {chain[1]}, {"anchor": chain[1], "root": chain[1]})]
            principal = "yl"
    if principal not in ("yl", "diyl", "ester_o") and principal_groups_left_outside(view, classes, principal, set(cand.mapping.values())):
        raise UnsupportedStructure("a principal characteristic group outside the parent needs a multiplicative or other parent")
    if "ester" in classes and "ester_o" in classes:
        raise UnsupportedStructure("esters of both an acid and an alcohol of the parent are not supported")
    branches = list(groups.branches)
    for cls in _ACYL_CLASSES:
        if cls in classes and cls != principal:
            for loc, _, extra in classes.pop(cls):
                if "root" not in extra:
                    raise UnsupportedStructure("a terminal acyl group below the principal group is not supported")
                branches.append((loc, extra["root"]))

    final = _primed(_final(skel), parent_prime)
    fused = [name_fused(c, cand, view, final, parent.centers, {}) for c in fused_comps]
    fusion_atoms = frozenset(loc for f in fused for loc in f.fusion_locs)
    enes_bonds, hydro, dehydro, retro = unsaturation(cand, view, fusion_atoms)
    cost += 1 if retro else 0
    shared = set(hydro) & set(dehydro)
    if shared:
        left_hydro = [x for x in hydro if x not in shared]
        left_dehydro = [x for x in dehydro if x not in shared]
        if not len(left_hydro) % 2 and not len(left_dehydro) % 2:
            hydro, dehydro = left_hydro, left_dehydro
    mancude_indicated = []
    if len(hydro) % 2 and (cand.skel.ops or cand.replaced):
        mancude_indicated = [min(hydro, key=loc_key)]
        hydro = list(hydro)
        hydro.remove(mancude_indicated[0])
    if len(hydro) % 2 or len(dehydro) % 2:
        raise UnsupportedStructure("indicated hydrogen would be needed")
    if len(enes_bonds) + len(hydro) // 2 + len(dehydro) // 2 > _MAX_UNSATURATION_CHANGES or (len(hydro) + len(dehydro)) // 2 > _MAX_HYDRO_PAIRS and not parent.name.endswith("carotene"):
        raise UnsupportedStructure("too many changes of the degree of hydrogenation for this parent")
    if _pair_count(hydro, dehydro) + cost > _MAX_TOTAL_COST and not parent.name.endswith("carotene"):
        raise UnsupportedStructure("too many modifications for this parent")
    cost += _pair_count(hydro, dehydro)
    modifications = len(skel.ops) + len(cand.cyclo) + len(cand.replaced) + len(bridge_comps) + len(fused_comps) + len(spiro_comps)
    if not view.has_stereo and (modifications > 1 or cost > 2):
        raise UnsupportedStructure("a heavily modified parent needs the configuration to be a natural product")
    hydro = sorted((final(x) for x in hydro), key=loc_key)
    dehydro = sorted((final(x) for x in dehydro), key=loc_key)
    enes = [locant_pair(skel, final, a, b) for a, b, o in enes_bonds if o == 2]
    ynes = [locant_pair(skel, final, a, b) for a, b, o in enes_bonds if o == 3]
    enes.sort(key=lambda t: loc_key(t.split("(")[0]))
    ynes.sort(key=lambda t: loc_key(t.split("(")[0]))

    config = configuration(cand, view)
    config.parent = [(key, _locant_through(text, final)) for key, text in config.parent]
    config.side = [(key, _locant_through(text, final)) for key, text in config.side]
    if spiro:
        _add_spiro_center(config, cand, view, spiro, spiro_comps[0], final, parent_prime)
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
    if spiro:
        if enes or ynes or hydro or dehydro:
            raise UnsupportedStructure("unsaturation of a spiro parent is not supported")
        for locant, name, compound in spiro.substituents:
            prefix_groups.setdefault(name, {"locants": [], "compound": compound})["locants"].append(Loc(f"{locant}{ring_prime}", ""))
        for entry in prefix_groups.values():
            entry["locants"].sort(key=lambda t: loc_key(t.base))
    prefix = format_substituent_prefixes(prefix_groups) if prefix_groups else ""

    suffix, locants, alkyl_word, anion = "", [], "", ""
    if principal is not None:
        members = sorted(classes[principal], key=lambda m: loc_key(m[0]))
        if principal == "yl":
            word = "yl"
        elif principal in ("acid", "ester", "amide"):
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
        if principal == "diyl":
            if len(members) != 1:
                raise UnsupportedStructure("several cyclic acetals on one parent are not supported")
            loc, _, extra = members[0]
            locants = sorted([str(Loc(final(loc), config.faces.get(extra["anchors"][0], ""))),
                              str(Loc(final(extra["second"]), config.faces.get(extra["anchors"][1], "")))], key=lambda t: loc_key(t.rstrip("αβξ")))
            suffix = "diyl"
            anion, alkyl_word = _acetal_words(view, extra)
        if principal == "ester_o":
            anions = {_acid_anion(view, extra["anchor"], extra["acyl"], mapped) for _, _, extra in members}
            if len(anions) != 1:
                raise UnsupportedStructure("O-acyl groups from different acids are not supported")
            (anion,) = anions
            if len(members) > 1:
                if not anion.isalpha():
                    raise UnsupportedStructure("a substituted acid part repeated on a natural product is not supported")
                anion = multiplying_prefix(len(members), compound=False) + anion

    stem_core = core_text(_parent_text_name(parent.name, cand.skel.ops), enes, ynes, suffix, locants)
    if spiro:
        center = final(spiro_center(cand, spiro_comps[0]))
        ring_loc = f"{spiro.spiro_locant}{ring_prime}"
        if parent_prime:
            joined = f"spiro[{spiro.ring_name}-{spiro.spiro_locant},{center}-{parent.name}]"
        else:
            joined = f"spiro[{parent.name}-{center},{ring_loc}-{spiro.ring_name}]"
        stem_core = joined + (f"-{','.join(locants)}-{suffix}" if suffix else "")
    descriptor = ",".join(text for _, text in sorted(config.parent)) + "-" if config.parent else ""
    groups = assign_primes(fused)
    fused_hydro = hydro_texts(groups)
    if fused_hydro:
        hydro = sorted(hydro + fused_hydro, key=lambda t: (int("".join(c for c in t if c.isdigit())), loc_key(t)))
    hydro_text = _hydro_text(dehydro, "dehydro") + _hydro_text(hydro, "hydro")
    cyclo_text = _cyclo_text(cand, final, config, view)
    ops = op_prefixes(cand, final, cyclo_text, retro)
    replacement = _replacement_text(cand, final, view)
    if replacement:
        ops = [replacement] + ops
    bridges = bridge_prefixes(bridge_comps, cand, view, config, final)
    if nondetachable_cost(cand) + len(ring_comps) > _MAX_SKELETAL_MODIFICATIONS:
        raise UnsupportedStructure("too many skeletal modifications")
    fused_prefix = fusion_text(groups) if fused else ""
    indicated = indicated_texts(groups)
    front_stereo, front_plain = [], []
    for f in fused:
        for locant, atom, role in f.fusion_h:
            face = config.hfaces.get(atom, "")
            if face:
                front_stereo.append((loc_key(locant), f"{locant}{FACES[face]}H"))
            elif role == "indicated":
                front_plain.append((loc_key(locant), f"{locant}H"))
    for loc in cand.replaced:
        atom = cand.mapping[loc]
        if view.elem[atom] == "C" and view.mol.GetAtomWithIdx(atom).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            face = config.hfaces.get(atom, "")
            if face:
                front_stereo.append((loc_key(final(loc)), f"{final(loc)}{FACES[face]}H"))
    indicated = indicated + [f"{final(loc)}H" for loc in mancude_indicated] + [t for _, t in sorted(front_plain)]
    indicated = sorted(indicated, key=lambda t: loc_key(t[:-1]))
    bridge_head = ""
    if parent.name.endswith("carotene") and bridges:
        bridge_head = "-".join(b.rstrip("-") for b in bridges) + "-"
        bridges = []
    ops = ([fused_prefix] if fused_prefix else []) + bridges + ops
    nondetachable = "-".join(p.rstrip("-") for p in ops) + ("-" if ops and (descriptor or ops[-1].endswith("-")) else "")
    if ops and not nondetachable.endswith("-") and not descriptor and stem_core and not stem_core[0].isascii():
        nondetachable += "-"
    side_items = sorted(config.side + front_stereo)
    side = f"({','.join(text for _, text in side_items)})-" if side_items else ""
    indicated_text = f"{','.join(indicated)}-" if indicated else ""
    detachable = f"{prefix}-" if prefix and (hydro_text or indicated_text or nondetachable or descriptor) else prefix
    if hydro_text and prefix:
        detachable = f"{prefix}-"
    if hydro_text.endswith("-") and not indicated_text and nondetachable and not nondetachable[0].isdigit():
        hydro_text = hydro_text[:-1]
    body = f"{side}{detachable}{bridge_head}{hydro_text}{indicated_text}{nondetachable}{descriptor}{stem_core}"
    name = " ".join(part for part in (alkyl_word, body, anion) if part)
    if chain is not None:
        name = chain[2].replace(chain[3], _enclose_group(body), 1)
    nondet = nondetachable_cost(cand)
    rearranged = any(op[0] == "seco" for op in cand.skel.ops) or bool(cand.cyclo)
    removed = sorted((loc_key(op[1])[1] for op in cand.skel.ops if op[0] == "nor"), reverse=True)
    inserted = sorted((_homo_position(op) for op in cand.skel.ops if op[0] == "homo"), reverse=True)
    first = next((text for _, text in sorted(config.parent) if text[-1] in "αβ"), "")
    return Built(
        name, cost, (cost, nondet > 0 or bool(fused_comps) or bool(bridge_comps) or bool(spiro_comps), -len(mapping), rearranged, tuple(-n for n in removed), tuple(-n for n in inserted), len(config.parent) + len(config.side)),
        config.implied_total, config.implied_cited, first[-1] if first else "",
    )


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


def _parent_text_name(name, ops):
    """Carotene names drop the designation of an end group removed by 'apo' beyond its ring (P-101.3.4.2)."""
    if name.endswith("-carotene") and any(op[0] == "apo" for op in ops):
        ends = name[: -len("-carotene")].split(",")
        return f"{ends[0]}-carotene"
    return name


def _acetal_words(view, extra):
    """(functional class word, leading carbonyl-compound name) of a cyclic acetal on two skeleton oxygens."""
    from .core import smiles_to_iupac

    carbon, others = extra["carbon"], extra["others"]
    oxo = [n for n in others if view.elem[n] == "O" and view.order[frozenset((carbon, n))] == 2]
    if oxo and len(others) == 1:
        return "carbonate", ""
    editable = Chem.RWMol(view.mol)
    keep = {carbon, *others}
    stack = list(others)
    while stack:
        a = stack.pop()
        for n in view.adj[a]:
            if n != carbon and n not in keep and n not in extra["anchors"]:
                keep.add(n)
                stack.append(n)
    for idx in sorted((i for i in view.adj if i not in keep), reverse=True):
        editable.RemoveAtom(idx)
    carbonyl = editable.GetMol()
    for atom in carbonyl.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    rw = Chem.RWMol(carbonyl)
    new_o = rw.AddAtom(Chem.Atom(8))
    c_index = sorted(keep).index(carbon)
    rw.AddBond(c_index, new_o, Chem.BondType.DOUBLE)
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    name = smiles_to_iupac(Chem.MolToSmiles(result))
    has_h = view.mol.GetAtomWithIdx(carbon).GetTotalNumHs() > 0
    return ("acetal" if has_h else "ketal"), name


def _primed(final, prime):
    if not prime:
        return final
    return lambda loc: f"{final(loc)}{prime}"


def spiro_center(cand, comp):
    image = {a: loc for loc, a in cand.mapping.items()}
    return image[comp.links[0][0]]


def _enclose_group(text):
    from ._substituents import wrap_marks

    return wrap_marks(text)


_LOCANT_TEXT = re.compile(r"^(\d+[a-c]?[¹²³]*)(.*)$")


def _locant_through(text, final):
    match = _LOCANT_TEXT.match(text)
    return f"{final(match.group(1))}{match.group(2)}" if match else text


def _add_spiro_center(config, cand, view, spiro, comp, final, parent_prime):
    center = spiro_center(cand, comp)
    atom = cand.mapping[center]
    if not any(e.centeredOn == atom and e.type == Chem.StereoType.Atom_Tetrahedral for e in Chem.FindPotentialStereo(view.mol)):
        return
    cited = spiro.spiro_locant if parent_prime else final(center)
    if view.mol.GetAtomWithIdx(atom).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
        label = "ξ"
    else:
        rdCIPLabeler.AssignCIPLabels(view.mol)
        label = view.mol.GetAtomWithIdx(atom).GetProp("_CIPCode")
    config.side.append((loc_key(cited), f"{cited}{label}"))


def _homo_position(op):
    """Locant number of the atom after which a methylene is inserted: equivalent connectors take the highest."""
    kind, data = op[1], op[2]
    if kind == "terminal":
        return loc_key(data)[1]
    if kind == "atomic":
        return loc_key(data[0])[1]
    return max(loc_key(data[0])[1], loc_key(data[1])[1])
