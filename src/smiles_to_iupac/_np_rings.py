"""Atoms and rings added to a natural-product parent: atomic bridges (P-101.5.2)."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import UnsupportedStructure
from ._np_core import FACES, loc_key
from ._numerals import multiplying_prefix

_BRIDGES = {
    ("O",): "epoxy",
    ("O", "O"): "epidioxy",
    ("S",): "epithio",
    ("S", "S"): "epidithio",
    ("N",): "epimino",
    ("C",): "methano",
    ("C", "C"): "ethano",
    ("C", "C", "C"): "propano",
    ("C", "C", "C", "C"): "butano",
}


@dataclass
class Component:
    atoms: set
    links: list


def components(view, mapped):
    """Connected groups of non-skeleton atoms attached to the skeleton at two or more bonds."""
    seen, found = set(mapped), []
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
        links = [(a, n) for n in comp for a in view.adj[n] if a in mapped]
        if len(links) > 1 and not _closes_acetal(comp, links, view):
            found.append(Component(comp, links))
    return found


def _closes_acetal(comp, links, view):
    """Two oxygens on skeleton atoms bonded to one carbon: a cyclic acetal, ketal or carbonate (P-101.7.4)."""
    if len(links) != 2 or len({a for a, _ in links}) != 2:
        return False
    oxygens = [n for _, n in links]
    if len(set(oxygens)) != 2 or any(view.elem[o] != "O" or len(view.adj[o]) != 2 for o in oxygens):
        return False
    shared = [set(view.adj[o]) - {a for a, n in links if n == o} for o in oxygens]
    return shared[0] == shared[1] and len(shared[0]) == 1 and view.elem[next(iter(shared[0]))] == "C"


def _path(comp, view):
    """The atoms of an unbranched acyclic component from one attached end to the other."""
    atoms = comp.atoms
    ends = {n for _, n in comp.links}
    if len(comp.links) != 2:
        raise UnsupportedStructure("a ring added to a natural-product parent at more than two bonds is not supported")
    if len(atoms) == 1:
        return list(atoms)
    inner = {a: [n for n in view.adj[a] if n in atoms] for a in atoms}
    if any(len(v) > 2 for v in inner.values()) or sum(len(v) for v in inner.values()) // 2 != len(atoms) - 1:
        raise UnsupportedStructure("a branched or cyclic added component is not supported")
    start = sorted(ends)[0]
    order, previous = [start], None
    while True:
        nxt = [n for n in inner[order[-1]] if n != previous]
        if not nxt:
            break
        previous = order[-1]
        order.append(nxt[0])
    if len(order) != len(atoms):
        raise UnsupportedStructure("a branched added component is not supported")
    return order


def bridge_prefixes(comps, cand, view, config, final):
    """['3α,8-epidioxy', ...] for the atomic bridges of a candidate, grouped like P-101.5.2 and P-14.2."""
    parent = cand.parent
    image = {a: loc for loc, a in cand.mapping.items()}
    grouped = {}
    for comp in comps:
        order = _path(comp, view)
        elements = tuple(view.elem[a] for a in order)
        word = _BRIDGES.get(elements)
        if word is None:
            raise UnsupportedStructure("this added bridge is not named")
        if elements[0] == "C" and not _fused_ok(parent):
            raise UnsupportedStructure("a carbon bridge needs a polycyclic parent")
        for a in order:
            extra = [n for n in view.adj[a] if n not in comp.atoms and n not in image]
            if extra or (view.elem[a] == "N" and view.mol.GetAtomWithIdx(a).GetTotalNumHs() != 1):
                raise UnsupportedStructure("a substituted bridge is not named")
        texts = []
        for end_atom, skeleton_atom in sorted(((n, a) for a, n in comp.links), key=lambda t: loc_key(final(image[t[1]]))):
            loc = image[skeleton_atom]
            stereo = view.mol.GetAtomWithIdx(skeleton_atom).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            face = config.faces.get(end_atom, "") if loc not in parent.centers and stereo else ""
            texts.append(f"{final(loc)}{FACES.get(face, '')}")
        grouped.setdefault(word, []).append(",".join(texts))
    result = []
    for word in sorted(grouped):
        pairs = sorted(grouped[word], key=lambda t: loc_key(t.split(",")[0].rstrip("αβξ")))
        mult = "" if len(pairs) == 1 else multiplying_prefix(len(pairs), compound=False)
        result.append(f"{':'.join(pairs)}-{mult}{word}")
    return result


def _fused_ok(parent):
    from ._np import _parent_facts

    return _parent_facts(parent.name)[4] >= 3


def split_components(comps, cand, view):
    """(bridges, fused rings): a short unbranched chain is a bridge unless it closes a ring on adjacent atoms."""
    image = {a: loc for loc, a in cand.mapping.items()}
    bridges, fused = [], []
    spiro = []
    for comp in comps:
        skeleton = {a for a, _ in comp.links}
        if len(skeleton) == 1 and len(comp.links) == 2:
            spiro.append(comp)
            continue
        if len(skeleton) != 2 or len(comp.links) != 2:
            raise UnsupportedStructure("an added ring attached at other than two skeleton atoms")
        sa, sb = sorted(skeleton)
        adjacent = image[sb] in cand.skel.adj[image[sa]]
        try:
            order = _path(comp, view)
            elements = tuple(view.elem[a] for a in order)
        except UnsupportedStructure:
            elements = None
        if elements in _BRIDGES and (not adjacent or elements[0] != "C"):
            bridges.append(comp)
        else:
            fused.append(comp)
    return bridges, fused, spiro


@dataclass
class Spiro:
    ring_name: str
    spiro_locant: str
    substituents: list


def name_spiro(comp, cand, view):
    """The saturated ring spiro-joined to one skeleton atom (P-101.5.3, P-24.5)."""
    from ._np_fusion import _mono_numberings, ring_cycle
    from ._substituents import name_branch
    from ._common import adjacency, halogen_substituents

    spiro_atom = comp.links[0][0]
    l1, l2 = comp.links[0][1], comp.links[1][1]
    inside = comp.atoms
    previous = {l1: None}
    frontier = [l1]
    while frontier and l2 not in previous:
        nxt = []
        for a in frontier:
            for n in view.adj[a]:
                if n in inside and n not in previous:
                    previous[n] = a
                    nxt.append(n)
        frontier = nxt
    path = [l2]
    while previous[path[-1]] is not None:
        path.append(previous[path[-1]])
    ring = set(path) | {spiro_atom}
    cycle = ring_cycle(ring, view.adj)
    if cycle is None or len(ring) > 8:
        raise UnsupportedStructure("a spiro component that is not a simple ring")
    if any(view.order[frozenset((a, b))] != 1 for a in ring for b in view.adj[a] if b in ring):
        raise UnsupportedStructure("an unsaturated spiro component is not supported")
    elem = {a: view.elem[a] for a in ring}
    subs = [(a, n) for a in path for n in view.adj[a] if n not in ring and n in inside]
    if any(view.elem[a] == "C" for a in ()):
        pass
    graph = adjacency(view.mol)
    halogens = halogen_substituents(view.mol)
    best = None
    for numbering in _mono_numberings(cycle, elem):
        spiro_number = numbering[spiro_atom]
        cites = sorted(numbering[a] for a, _ in subs)
        key = (spiro_number, cites)
        if best is None or key < best[0]:
            best = (key, numbering)
    numbering = best[1]
    entries = []
    for a, n in subs:
        name, compound = name_branch(graph, n, a, halogens, frozenset(), mol=view.mol, unsaturated=True)
        entries.append((str(numbering[a]), name, compound))
    return Spiro(_ring_name(ring, view, elem), str(numbering[spiro_atom]), entries)


def _ring_name(ring, view, elem):
    from .core import smiles_to_iupac

    editable = Chem.RWMol(view.mol)
    for idx in sorted((i for i in view.adj if i not in ring), reverse=True):
        editable.RemoveAtom(idx)
    free = editable.GetMol()
    for atom in free.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    Chem.SanitizeMol(free)
    return smiles_to_iupac(Chem.MolToSmiles(free))
