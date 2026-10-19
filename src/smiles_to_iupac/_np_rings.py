"""Atoms and rings added to a natural-product parent: atomic bridges (P-101.5.2)."""

from dataclasses import dataclass

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_bridged import GENERAL_BRIDGE_WORDS, Part, _kekule, _prefix_text, _units_of, bridge_assemblies
from ._np_core import FACES, loc_key


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


def _backbone(comp, view):
    """The atoms of the bridge itself: the shortest chain between the atoms attached to the skeleton and every ring that
    chain runs through; the other atoms of the component are substituents of the bridge."""
    graph = nx.Graph([(a, b) for a in comp.atoms for b in view.adj[a] if b in comp.atoms])
    graph.add_nodes_from(comp.atoms)
    ends = sorted({n for _, n in comp.links})
    through = set(nx.shortest_path(graph, ends[0], ends[-1]))
    rings = []
    for cycle in nx.cycle_basis(graph):
        merged = set(cycle)
        for ring in rings[:]:
            if len(ring & merged) >= 2:
                merged |= ring
                rings.remove(ring)
        rings.append(merged)
    for ring in rings:
        if len(ring & through) >= 2:
            through |= ring
    return through


@dataclass
class Bridges:
    mol: Chem.Mol
    atoms: set
    entries: list
    numbering: dict
    backbone: set


def bridge_layout(comps, cand, view, final):
    """The words, attachment locants and atom numbers of the atomic bridges of a candidate, the bridge atoms numbered on
    from the highest locant of the parent (P-101.5.2, P-25.4)."""
    if not comps:
        return Bridges(view.mol, set(), [], {}, set())
    image = {a: loc for loc, a in cand.mapping.items()}
    mapped = set(image)
    if any(view.elem[a] == "C" for comp in comps for a in _backbone(comp, view)) and not _fused_ok(cand.parent):
        raise UnsupportedStructure("a carbon bridge across an acyclic part of a parent that is not polycyclic")
    kekule = _kekule(view.mol)
    if kekule is None:
        raise UnsupportedStructure("an added bridge whose bonds cannot be written with alternating double bonds")
    graph = nx.Graph([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in kekule.GetBonds()])
    backbone = {}
    parts = []
    for comp in comps:
        backbone[id(comp)] = _backbone(comp, view)
        parts.append(Part(frozenset(backbone[id(comp)]), _units_of(kekule, graph, backbone[id(comp)]), True, 2))
    atoms = mapped | set().union(*backbone.values())
    position = {a: str(final(loc)) for a, loc in image.items()}
    token = GENERAL_BRIDGE_WORDS.set(True)
    try:
        assemblies = list(bridge_assemblies(kekule, atoms, mapped, position, parts))
    finally:
        GENERAL_BRIDGE_WORDS.reset(token)
    _, entries, numbering = min(assemblies, key=lambda found: found[0])
    return Bridges(kekule, atoms, entries, numbering, set().union(*backbone.values()))


def bridge_prefixes(layout, cand, view, config, final):
    """['3α,8-epidioxy', ...] for the atomic bridges of a candidate, grouped like P-101.5.2 and P-14.2."""
    parent = cand.parent
    image = {a: loc for loc, a in cand.mapping.items()}
    entries = []
    for text, locants, option, part in layout.entries:
        cited = []
        for token in option["cited"]:
            loc = image[token.outer]
            stereo = view.mol.GetAtomWithIdx(token.outer).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            face = config.faces.get(token.inner, "") if loc not in parent.centers and stereo else ""
            cited.append(f"{final(loc)}{FACES.get(face, '')}")
        entries.append((text, cited, option, part))
    words = sorted({(not part.independent, text) for text, _, _, part in entries}, key=lambda k: (not k[0], k[1].lstrip("0123456789,[]")))
    return [_prefix_text([e for e in entries if (not e[3].independent, e[0]) == word]) for word in words]


def _fused_ok(parent):
    from ._np import _parent_facts

    return _parent_facts(parent.name)[4] >= 3


def _peri_fused(comp, image, sa, sb, cand, view):
    """Whether an unsaturated component closes a ring on the two bonds at a ring-junction atom of the skeleton, which is
    named by peri-fusion (P-25.3.1.3) rather than by a bridge."""
    middle = cand.skel.adj[image[sa]] & cand.skel.adj[image[sb]]
    if len(middle) != 1:
        return False
    junction = sum(1 for ring in nx.cycle_basis(nx.Graph([(a, b) for a, n in cand.skel.adj.items() for b in n])) if next(iter(middle)) in ring)
    return junction >= 2 and any(
        view.order[frozenset((a, b))] != 1 for a in comp.atoms | {sa, sb} for b in view.adj[a] if b in comp.atoms | {sa, sb} and a < b
    )


def split_components(comps, cand, view):
    """(bridges, fused rings): a component joining two nonadjacent skeleton atoms is a bridge unless it closes a ring on the
    two bonds at a ring junction, as is a single oxygen or sulfur atom joining adjacent ones; any other component closes a
    ring on the bond between adjacent atoms (P-101.5.2)."""
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
        backbone = _backbone(comp, view)
        single_atom = len(backbone) == 1 and view.elem[next(iter(backbone))] in ("O", "S")
        if (not adjacent and not _peri_fused(comp, image, sa, sb, cand, view)) or (adjacent and single_atom):
            bridges.append(comp)
        else:
            ring_atoms = cand.skel.cycle_atoms()
            if not (image[sa] in ring_atoms and image[sb] in ring_atoms):
                raise UnsupportedStructure("a ring cannot be fused across an acyclic part of the parent")
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
