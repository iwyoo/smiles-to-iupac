"""Divalent and polyvalent linking-group components of multiplicative names
(P-15.3.1.2, P-29): carbon chains (principal chain by P-44.3, free valences
first), single and homonuclear heteroatom groups, and ring components.
"""

from dataclasses import dataclass

from ._common import ENE_BOND_ORDER, UnsupportedStructure, alpha_sort_key, multiplied_word, suffix_body
from ._multiplicative_prefix import SIMPLE_PREFIXES, prefix_name, subtree
from ._multiplicative_ring import name_ring_component
from ._numerals import alkane_name
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes

_SINGLE_ATOM_WORDS = {8: "oxy", 16: "sulfanediyl", 34: "selanediyl", 52: "tellanediyl", 7: "azanediyl"}
_SUBSTITUTABLE_WORDS = {
    7: ("azanediyl", "nitrilo"),
    5: ("boranediyl", "boranetriyl"),
    14: ("silanediyl", "silanetriyl"),
    15: ("phosphanediyl", "phosphanetriyl"),
    32: ("germanediyl", "germanetriyl"),
    33: ("arsanediyl", "arsanetriyl"),
}
_CHAIN_STEMS = {14: "silane"}
_HOMO_RUN_WORDS = {(8, 2): "peroxy", (16, 2): "disulfanediyl", (34, 2): "diselanediyl", (52, 2): "ditellanediyl"}


class DecompositionRejected(Exception):
    pass


@dataclass
class Part:
    text: str
    has_prefix: bool
    has_locants: bool


def _pendants(mol, atoms, attachment_pairs):
    atom_set = set(atoms)
    result = []
    for a in atoms:
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            idx = n.GetIdx()
            if idx not in atom_set and (a, idx) not in attachment_pairs:
                result.append((a, idx))
    return result


def _entry(mol, owner, root, ctx):
    if ctx.entry is not None:
        return ctx.entry(mol, owner, root, ctx)
    return _substituent_entry(mol, owner, root, ctx)


def _substituent_entry(mol, owner, root, ctx):
    atom = mol.GetAtomWithIdx(root)
    bond = mol.GetBondBetweenAtoms(owner, root)
    if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and bond.GetBondTypeAsDouble() == 2:
        return "oxo", False
    atoms = subtree(mol, root, owner)
    for g in ctx.groups:
        if g.anchor == root and g.atoms - {owner} == atoms and g.name in SIMPLE_PREFIXES and g.name != "ketone":
            if g.name in ("amide", "sulfonamide", "amine") and not _primary(mol, g):
                break
            return SIMPLE_PREFIXES[g.name], False
    return prefix_name(mol, root, owner, ctx.suffix_group, ctx.name_function, ctx.groups)


def _primary(mol, group):
    hetero = [a for a in group.atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    return len(hetero) == 1 and mol.GetAtomWithIdx(hetero[0]).GetTotalNumHs() == 2


def _hydride_chain_part(mol, atoms, attachments, pend, ctx, stem):
    ends = [a for a, _, _ in attachments]
    walk, seen = [ends[0]], {ends[0]}
    while len(walk) < len(atoms):
        nxt = [n.GetIdx() for n in mol.GetAtomWithIdx(walk[-1]).GetNeighbors() if n.GetIdx() in atoms and n.GetIdx() not in seen]
        if len(nxt) != 1:
            raise UnsupportedStructure("a branched heteroatom chain is not supported as a multiplicative linker")
        walk.append(nxt[0])
        seen.add(nxt[0])
    best = None
    for order in (walk, walk[::-1]):
        position = {a: i + 1 for i, a in enumerate(order)}
        cited = sorted((position[a], *_entry(mol, a, root, ctx)) for a, root in pend)
        key = (tuple(p for p, _, _ in cited), tuple(n for _, n, _ in cited))
        if best is None or key < best[0]:
            best = (key, cited)
    grouped = {}
    for p, name, compound in best[1]:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(p)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    n = len(atoms)
    body = f"{multiplied_word(n, stem)}-1,{n}-diyl"
    return Part(prefix + body, bool(prefix), True)


def _hetero_part(mol, atoms, attachments, ctx):
    pairs = {(a, b) for a, b, _ in attachments}
    pend = _pendants(mol, atoms, pairs)
    z = mol.GetAtomWithIdx(atoms[0]).GetAtomicNum()
    if len(atoms) >= 2:
        if z in _CHAIN_STEMS and len(attachments) == 2 and all(mol.GetAtomWithIdx(a).GetAtomicNum() == z for a in atoms):
            return _hydride_chain_part(mol, atoms, attachments, pend, ctx, _CHAIN_STEMS[z])
        word = _HOMO_RUN_WORDS.get((z, len(atoms)))
        if len(atoms) != 2 or word is None or pend or len(attachments) != 2:
            raise UnsupportedStructure("this heteroatom chain is not supported as a multiplicative linker")
        return Part(word, False, False)
    count = len(attachments)
    if z in _SUBSTITUTABLE_WORDS:
        divalent, trivalent = _SUBSTITUTABLE_WORDS[z]
        if count == 3 and not pend:
            return Part(trivalent, False, False)
        if count != 2:
            raise UnsupportedStructure("this heteroatom linker is not supported")
        entries = [_entry(mol, atoms[0], root, ctx) for _, root in pend]
        prefix = format_mononuclear_prefixes(entries) if entries else ""
        return Part(prefix + divalent, bool(entries), False)
    if count != 2:
        raise UnsupportedStructure("this heteroatom linker is not supported")
    oxo = [r for _, r in pend if mol.GetAtomWithIdx(r).GetAtomicNum() == 8 and mol.GetAtomWithIdx(r).GetDegree() == 1]
    if len(oxo) != len(pend):
        raise UnsupportedStructure("a substituted chalcogen linker is not supported")
    if z == 16 and len(oxo) == 1:
        return Part("sulfinyl", False, False)
    if z == 16 and len(oxo) == 2:
        return Part("sulfonyl", False, False)
    if oxo or z not in _SINGLE_ATOM_WORDS:
        raise UnsupportedStructure("this heteroatom linker is not supported")
    return Part(_SINGLE_ATOM_WORDS[z], False, False)


def _carbon_closure(mol, atoms):
    seen = set(atoms)
    stack = list(atoms)
    while stack:
        a = stack.pop()
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            idx = n.GetIdx()
            if idx in seen or n.GetAtomicNum() != 6 or n.IsInRing():
                continue
            seen.add(idx)
            stack.append(idx)
    return seen


def _tree_path(adj, s, t):
    parent = {s: None}
    stack = [s]
    while stack:
        a = stack.pop()
        for n in adj[a]:
            if n not in parent:
                parent[n] = a
                stack.append(n)
    path = [t]
    while path[-1] != s:
        path.append(parent[path[-1]])
    return list(reversed(path))


def _chain_descriptors(mol, ctx, position):
    result = []
    for kind, idx, code in ctx.stereo:
        if kind == "atom" and idx in position:
            result.append((position[idx], code, kind, idx))
        elif kind == "bond":
            bond = mol.GetBondWithIdx(idx)
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in position and b in position:
                result.append((min(position[a], position[b]), code, kind, idx))
    return result


def _carbon_part(mol, atoms, attachments, directed, ctx):
    pairs = {(a, b) for a, b, _ in attachments}
    if len(atoms) == 1 and len(attachments) == 2:
        pend = _pendants(mol, atoms, pairs)
        if len(pend) == 1:
            o = mol.GetAtomWithIdx(pend[0][1])
            if o.GetAtomicNum() == 8 and o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atoms[0], o.GetIdx()).GetBondTypeAsDouble() == 2:
                return Part("carbonyl", False, False)
    universe = _carbon_closure(mol, atoms)
    adj = {a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in universe] for a in universe}
    free_atoms = [a for a, _, _ in attachments]
    required = set(free_atoms)
    ordered = sorted(universe)
    candidates = []
    for i, s in enumerate(ordered):
        for t in ordered[i:]:
            path = _tree_path(adj, s, t)
            if required <= set(path):
                candidates.append(path)
    if not candidates:
        raise DecompositionRejected("free valences of a carbon linker do not lie on one chain")
    best = None
    for chain in candidates:
        bond_orders = [mol.GetBondBetweenAtoms(chain[i], chain[i + 1]).GetBondTypeAsDouble() for i in range(len(chain) - 1)]
        multiple = sum(1 for o in bond_orders if o > 1)
        for direction in (chain, list(reversed(chain))):
            position = {a: i + 1 for i, a in enumerate(direction)}
            if directed is not None:
                free_key = (position[directed[0]], position[directed[1]])
            else:
                free_key = tuple(sorted(position[a] for a in free_atoms))
            ene, yne = [], []
            for i in range(len(direction) - 1):
                order = mol.GetBondBetweenAtoms(direction[i], direction[i + 1]).GetBondTypeAsDouble()
                if order == 2:
                    ene.append(i + 1)
                elif order == 3:
                    yne.append(i + 1)
            entries = []
            on_chain = set(direction)
            for a in direction:
                for n in mol.GetAtomWithIdx(a).GetNeighbors():
                    idx = n.GetIdx()
                    if idx in on_chain or (a, idx) in pairs:
                        continue
                    if (a, idx) not in ctx.entry_cache:
                        ctx.entry_cache[(a, idx)] = _entry(mol, a, idx, ctx)
                    name, compound = ctx.entry_cache[(a, idx)]
                    entries.append((position[a], name, compound))
            key = (
                -len(direction),
                -multiple,
                free_key,
                tuple(sorted(ene + yne)),
                tuple(sorted(p for p, _, _ in entries)),
                tuple(p for p, _ in sorted(((p, n) for p, n, _ in entries), key=lambda e: (alpha_sort_key(e[1]), e[0]))),
                tuple(sorted(loc for loc, code, _, _ in _chain_descriptors(mol, ctx, position) if code in "ZR")),
            )
            if best is None or key < best[0]:
                best = (key, direction, position, ene, yne, entries)
    _, chain, position, ene, yne, entries = best
    length = len(chain)
    descriptors = _chain_descriptors(mol, ctx, position)
    for _, _, kind, idx in descriptors:
        ctx.used.add((kind, idx))
    stereo_prefix = f"({','.join(f'{loc}{code}' for loc, code, _, _ in sorted(descriptors))})-" if descriptors else ""
    if directed is not None:
        cited = [position[directed[1]], position[directed[0]]]
    else:
        cited = sorted(position[a] for a in free_atoms)
    count = len(cited)
    grouped = {}
    for p, name, compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": compound})
        info["locants"].append(p)
    if length == 1 and not ene and not yne:
        prefix = format_mononuclear_prefixes([(n, c) for _, n, c in entries]) if entries else ""
        base = {2: "methylene", 3: "methanetriyl", 4: "methanetetrayl"}.get(count)
        if base is None:
            raise UnsupportedStructure("unsupported carbon linker valence")
        return Part(prefix + base, bool(entries), False)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    word = multiplied_word(count, "yl")
    sorted_join = ",".join(str(x) for x in sorted(cited))
    cited_join = ",".join(str(x) for x in cited)
    if not ene and not yne:
        name = f"{alkane_name(length)}-{cited_join}-{word}"
    elif length == 2 and len(ene + yne) == 1:
        name = f"{'ethene' if ene else 'ethyne'}-{cited_join}-{word}"
    else:
        body, _ = suffix_body(ene, yne, word, sorted(cited))
        tail = f"{sorted_join}-{word}"
        body = body[: -len(tail)] + f"{cited_join}-{word}"
        stem = alkane_name(length)[:-3]
        needs_a = (len(ene) >= 2) if ene else (len(yne) >= 2)
        name = f"{stem}{'a' if needs_a else ''}-{body}"
    return Part(stereo_prefix + prefix + name, bool(prefix or stereo_prefix), True)


def name_component(mol, kind, atoms, attachments, ctx, directed=None):
    """`attachments`: [(atom_in_component, external_atom, bond_order)];
    `directed`: (unit_side_atom, center_side_atom) for an arm component."""
    if kind == "ring":
        ring_attachments = [(a, b) for a, b, _ in attachments]
        result = name_ring_component(mol, atoms, ring_attachments, ctx.groups, ctx.suffix_group, ctx.name_function, directed)
        if result is None:
            raise UnsupportedStructure("this ring is not supported as a multiplicative linker component")
        text, has_prefix = result
        return Part(text, has_prefix, True)
    if any(order != 1 for _, _, order in attachments):
        raise UnsupportedStructure("a multiple bond to the multiplied units is not supported yet")
    if kind == "carbon":
        return _carbon_part(mol, list(atoms), attachments, directed, ctx)
    return _hetero_part(mol, list(atoms), attachments, ctx)
