"""Fusion name of a carbocyclic ring system made of one parent component and ortho-fused benzo/cyclopenta/... rings
(P-25.3.1.2, P-25.3.1.3): attached components are cited as prefixes in alphabetical order, each followed by the letters
of the parent's peripheral bonds it occupies, as low as possible; the parent is phenanthrene or the monocycle
[7]annulene (the parents of systems no retained-name module covers). Indicated hydrogen comes from the numbering of the
whole system (P-25.3.3).
"""

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit import Chem

from ._common import UnsupportedStructure

_PREFIX_OF_SIZE = {3: "cyclopropa", 4: "cyclobuta", 5: "cyclopenta", 6: "benzo", 7: "cyclohepta", 8: "cycloocta"}
_MULTIPLIER = {1: "", 2: "di", 3: "tri", 4: "tetra"}
_LETTERS = "abcdefghijklmnopqrstuvwxyz"


def _phenanthrene_parent():
    from ._ring_diyl_numbering import _RETAINED_NUMBERINGS

    locants, _, interior = _RETAINED_NUMBERINGS["phenanthrene"][:3]
    index = {loc: i for i, loc in enumerate(locants)}
    graph = nx.Graph()
    graph.add_nodes_from(range(len(locants)))
    graph.add_edges_from((i, (i + 1) % len(locants)) for i in range(len(locants)))
    graph.add_edges_from((index[a], index[b]) for a, b in interior)
    return "phenanthrene", locants, graph


def _annulene_parent(size):
    locants = tuple(str(i) for i in range(1, size + 1))
    graph = nx.cycle_graph(size)
    return f"[{size}]annulene", locants, graph


def _ring_graph(mol):
    graph = nx.Graph()
    graph.add_nodes_from(a.GetIdx() for a in mol.GetAtoms())
    graph.add_edges_from((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds())
    return graph


def _attached_rings(system, base_atoms, peripheral_pairs):
    """[(ring size, frozenset(u, v))] for the components outside the parent, or None when any is not an ortho-fused
    path closing on a peripheral bond between two atoms that lie in a single parent ring."""
    rest = set(system) - set(base_atoms)
    attached = []
    seen = set()
    for start in rest:
        if start in seen:
            continue
        component, stack = {start}, [start]
        while stack:
            for n in system[stack.pop()]:
                if n in rest and n not in component:
                    component.add(n)
                    stack.append(n)
        seen |= component
        sub = system.subgraph(component)
        if any(d > 2 for _, d in sub.degree()) or not nx.is_connected(sub):
            return None
        ends = [a for a in component if sub.degree(a) <= 1] if len(component) > 1 else [start]
        if len(component) > 1 and len(ends) != 2:
            return None
        contacts = []
        for atom in component:
            base_neighbors = [n for n in system[atom] if n in base_atoms]
            if atom in ends:
                if len(component) == 1:
                    if len(base_neighbors) != 2:
                        return None
                elif len(base_neighbors) != 1:
                    return None
            elif base_neighbors:
                return None
            contacts.extend(base_neighbors)
        pair = frozenset(contacts)
        if len(pair) != 2 or pair not in peripheral_pairs:
            return None
        attached.append((len(component) + 2, pair))
    return attached


def _embeddings(system, parent):
    name, locants, graph = parent
    nonfusion = {i for i in graph if graph.degree(i) == 2}
    periphery = [frozenset((i, (i + 1) % len(locants))) for i in range(len(locants))]
    letter_of = {pair: _LETTERS[i] for i, pair in enumerate(periphery)}
    usable = {pair for pair in periphery if all(a in nonfusion for a in pair)}
    for mapping in GraphMatcher(system, graph).subgraph_isomorphisms_iter():
        to_system = {b: s for s, b in mapping.items()}
        base_atoms = set(mapping)
        pairs = {frozenset(to_system[a] for a in pair): pair for pair in usable}
        attached = _attached_rings(system, base_atoms, set(pairs))
        if attached is None:
            continue
        found = []
        for size, pair in attached:
            if size not in _PREFIX_OF_SIZE:
                break
            found.append((_PREFIX_OF_SIZE[size], letter_of[pairs[pair]]))
        else:
            yield name, found, base_atoms


_ANTHRACENE = nx.Graph(
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (4, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 0)]
)


def _name(mol, system):
    parents = [_phenanthrene_parent()]
    if any(len(r) == 7 for r in mol.GetRingInfo().AtomRings()):
        parents.append(_annulene_parent(7))
    for parent in parents:
        options = [(sorted(letter for _, letter in found), found, name) for name, found, _ in _embeddings(system, parent)]
        if options:
            letters, found, name = min(options, key=lambda option: option[0])
            return name, found, letters
    return None


def _prefix_text(name, found):
    by_prefix = {}
    for prefix, letter in found:
        by_prefix.setdefault(prefix, []).append(letter)
    parts = []
    for prefix in sorted(by_prefix):
        letters = sorted(by_prefix[prefix])
        word = _MULTIPLIER[len(letters)] + prefix
        parts.append(word + f"[{','.join(letters)}]")
    if name.startswith("[") and len(found) == 1:
        parts = [parts[0].split("[")[0]]
    return "".join(parts)


def _candidate(mol):
    if mol.GetNumAtoms() < 10 or any(a.GetAtomicNum() != 6 or a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    ring_info = mol.GetRingInfo()
    if mol.GetNumAtoms() != len({a for r in ring_info.AtomRings() for a in r}) or len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(len(r) not in (5, 6, 7) for r in ring_info.AtomRings()):
        return None
    system = _ring_graph(mol)
    found = _name(mol, system)
    if found is None:
        return None
    if found[0] == "phenanthrene" and GraphMatcher(system, _ANTHRACENE).subgraph_is_isomorphic():
        return None
    return found


def name_fusion_prefix_system(mol) -> str:
    from ._fusion_numbering_oriented import oriented_peripheral_numberings
    from ._ring_diyl_numbering import _retained_numberings

    candidate = _candidate(mol)
    if candidate is None:
        raise UnsupportedStructure("no parent component with attached benzo or cyclopenta rings fits this ring system")
    base, found, _ = candidate
    prefixes = _prefix_text(base, found)
    sp3 = [a.GetIdx() for a in mol.GetAtoms() if a.GetTotalNumHs() == 2]
    indicated = ""
    if sp3:
        if base == "phenanthrene" and [p for p, _ in found] == ["cyclopenta"] and found[0][1] == "a":
            numberings = _retained_numberings(mol, "cyclopenta[a]phenanthrene")
        else:
            numberings = oriented_peripheral_numberings(mol, ignore_indicated=False)
        if not numberings:
            raise UnsupportedStructure("this ring system has no verified numbering for its indicated hydrogen")
        best = min(sorted(int(n[a].rstrip("abcdefgh")) for a in sp3) for n in numberings)
        indicated = ",".join(f"{loc}H" for loc in best) + "-"
    return f"{indicated}{prefixes}{base}"
