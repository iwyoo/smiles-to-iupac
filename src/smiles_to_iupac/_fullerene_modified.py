"""Fullerenes of the numbered C60-Ih and C70-D5h(6) cages with structure-modifying and replacement prefixes (P-27.4,
P-27.5): 'nor' deletes cage atoms, 'seco' cleaves cage bonds, 'homo' inserts a methylene group into a bond, 'cyclo' forms
a bond and 'aza', 'sila', ... replace cage atoms. The molecule is embedded into the cage graph with as few operations as
possible; the numbering is the one of the cage that gives the lowest locants (P-27.4.5, P-27.5.2)."""

from dataclasses import dataclass, field
from functools import lru_cache
from itertools import combinations

import networkx as nx
from rdkit import Chem

from ._fullerene import _FULLERENE_C60_SMILES, _FULLERENE_C70_SMILES, _cage_graph
from ._fullerene_numbering import _Cage, fullerene_numberings
from ._numerals import multiplying_prefix

_CAGES = (
    (_FULLERENE_C60_SMILES, "(C60-Ih)[5,6]fullerene", 31, 11),
    (_FULLERENE_C70_SMILES, "(C70-D5h(6))[5,6]fullerene", 36, 13),
)
_MAX_HOMO = 2
_MAX_CYCLO = 2
_MAX_SOLUTIONS = 4000
_REPLACEMENT_PREFIXES = {
    8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma",
    14: "sila", 32: "germa", 50: "stanna", 82: "plumba", 5: "bora", 13: "aluma", 31: "galla", 49: "indiga", 81: "thalla",
}
_SENIORITY = {z: rank for rank, z in enumerate(_REPLACEMENT_PREFIXES)}


class _Reference:
    def __init__(self, smiles, stem, atoms, rings):
        mol = Chem.MolFromSmiles(smiles)
        self.stem = stem
        self.least_atoms = atoms
        self.least_rings = rings
        self.graph = _cage_graph(mol, range(mol.GetNumAtoms()))
        self.adjacent = {a: frozenset(self.graph[a]) for a in self.graph}
        self.faces = [tuple(face) for face in _Cage(self.graph).faces]
        self.numberings = fullerene_numberings(self.graph)
        base = self.numberings[0]
        inverse = [{locant: atom for atom, locant in numbering.items()} for numbering in self.numberings]
        orbits = {frozenset(inv[base[a]] for inv in inverse) for a in self.graph}
        self.representatives = sorted(min(orbit) for orbit in orbits)
        self.orbit_count = len(orbits)


@lru_cache(maxsize=None)
def _references():
    return tuple(_Reference(*cage) for cage in _CAGES)


@dataclass
class _Edit:
    """The operations that turn a cage into the molecule, as cage atoms and bonds."""

    deleted: list = field(default_factory=list)
    seco: list = field(default_factory=list)
    cyclo: list = field(default_factory=list)
    homo: list = field(default_factory=list)
    replaced: dict = field(default_factory=dict)
    saturated: list = field(default_factory=list)

    @property
    def cost(self):
        return len(self.deleted) + len(self.seco) + len(self.cyclo) + len(self.homo)


def _skeleton(mol):
    """(graph, hydrogens) of the heavy atoms when they can form a modified cage, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or not 31 <= mol.GetNumAtoms() <= 90:
        return None
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 10 or sum(1 for a in mol.GetAtoms() if not ring_info.NumAtomRings(a.GetIdx())) > _MAX_HOMO:
        return None
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() != 6 and atom.GetAtomicNum() not in _REPLACEMENT_PREFIXES
            or atom.GetFormalCharge()
            or atom.GetNumRadicalElectrons()
            or atom.GetIsotope()
            or atom.GetDegree() > 3
            or any(b.GetBondTypeAsDouble() > 2.0 for b in atom.GetBonds())
        ):
            return None
    graph = nx.Graph((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds())
    return graph, {a.GetIdx(): a.GetTotalNumHs() for a in mol.GetAtoms()}


def _embeddings(graph, reference):
    """Injective maps of `graph` into the cage in which every bond of `graph` is a cage bond."""
    if graph.number_of_nodes() > reference.graph.number_of_nodes() or not nx.is_connected(graph):
        return []
    seed = max(graph, key=lambda a: (graph.degree(a), -a))
    order, parent = [seed], {seed: None}
    for u, v in nx.bfs_edges(graph, seed):
        order.append(v)
        parent[v] = u
    cage = reference.adjacent
    neighbours = {a: frozenset(graph[a]) for a in graph}
    found, placed, taken = [], {}, set()

    def extend(i):
        if len(found) >= _MAX_SOLUTIONS:
            return
        if i == len(order):
            found.append(dict(placed))
            return
        v = order[i]
        for g in cage[placed[parent[v]]]:
            if g not in taken and all(placed[w] in cage[g] for w in neighbours[v] if w in placed):
                placed[v] = g
                taken.add(g)
                extend(i + 1)
                del placed[v]
                taken.discard(g)

    for representative in reference.representatives:
        placed[seed] = representative
        taken.add(representative)
        extend(1)
        del placed[seed]
        taken.discard(representative)
    return found


def _cyclo_choices(graph, forced, count):
    """The connected graphs left when `count` bonds that are not `forced` are taken out of `graph`."""
    for removed in combinations([e for e in graph.edges if frozenset(e) not in forced], count):
        rest = graph.copy()
        rest.remove_edges_from(removed)
        if nx.is_connected(rest):
            yield rest


def name_modified_fullerene(mol):
    """Name of `mol` as a nor-, seco-, homo-, cyclo- or heteroatom-replaced C60-Ih or C70-D5h(6) fullerene, else None. A
    fragment keeps the fullerene name only with more than half of the atoms and a third of the rings of the cage
    (P-52.2.6.1)."""
    skeleton = _skeleton(mol)
    if skeleton is None:
        return None
    graph, hydrogens = skeleton
    element = {a: mol.GetAtomWithIdx(a).GetAtomicNum() for a in graph}
    if any(element[a] == 6 and graph.degree(a) + hydrogens[a] not in (3, 4) for a in graph):
        return None
    if sum(1 for a in graph if hydrogens[a] and (graph.degree(a) + hydrogens[a] == 4 or element[a] != 6)) > 1:
        return None
    homo_candidates = [
        x for x in graph
        if graph.degree(x) == 2
        and not graph.has_edge(*graph[x])
        and all(graph.degree(n) == 3 for n in graph[x])
        and (hydrogens[x] == 2 or element[x] != 6)
    ]
    attempts = []
    for homo_count in range(min(_MAX_HOMO, len(homo_candidates)) + 1):
        for homo in combinations(homo_candidates, homo_count):
            reduced, forced = graph.copy(), set()
            for x in homo:
                a, b = graph[x]
                reduced.remove_node(x)
                reduced.add_edge(a, b)
                forced.add(frozenset((a, b)))
            for reference in _references():
                deleted = reference.graph.number_of_nodes() - reduced.number_of_nodes()
                if deleted >= 0 and reduced.number_of_nodes() >= reference.least_atoms:
                    attempts.append((deleted + homo_count, len(attempts), reference, homo, reduced, forced))
    best = None
    for lower, _, reference, homo, reduced, forced in sorted(attempts, key=lambda a: a[:2]):
        if best is not None and lower > best[0]:
            break
        deficiency = sum(3 - d for _, d in reduced.degree)
        deleted = reference.graph.number_of_nodes() - reduced.number_of_nodes()
        for cyclo_count in range(_MAX_CYCLO + 1):
            cleaved = max(0, -(-(deficiency + 2 * cyclo_count - 3 * deleted) // 2))
            if best is not None and lower + cyclo_count + cleaved > best[0]:
                continue
            for rest in _cyclo_choices(reduced, forced, cyclo_count):
                for mapping in _embeddings(rest, reference):
                    found = _edit(graph, hydrogens, element, reference, mapping, homo)
                    if found is None:
                        continue
                    candidate = _name(reference, *found)
                    if candidate is not None and (best is None or candidate[:2] < best[:2]):
                        best = candidate
    return None if best is None else best[2]


def _edit(graph, hydrogens, element, reference, mapping, homo):
    """(edit, number of intact cage rings) of an embedding, None when the hydrogens do not fit a name."""
    cage = reference.graph
    image = set(mapping.values())
    inverse = {g: a for a, g in mapping.items()}
    edit = _Edit(deleted=[g for g in cage if g not in image])
    homo_bonds = {x: tuple(mapping[n] for n in graph[x]) for x in homo}
    bridged = {frozenset(bond) for bond in homo_bonds.values()}
    for g, h in cage.edges:
        if g in image and h in image and not graph.has_edge(inverse[g], inverse[h]) and frozenset((g, h)) not in bridged:
            edit.seco.append((g, h))
    edit.cyclo = [(mapping[a], mapping[b]) for a, b in graph.edges if a in mapping and b in mapping and not cage.has_edge(mapping[a], mapping[b])]
    for a, g in mapping.items():
        if element[a] != 6:
            if hydrogens[a]:
                return None
            edit.replaced[g] = element[a]
        elif graph.degree(a) + hydrogens[a] == 4:
            edit.saturated.append(g)
    for x, bond in homo_bonds.items():
        edit.homo.append((bond, element[x]))
        if hydrogens[x]:
            edit.saturated.append(bond)
    if len(edit.saturated) > 1 or not (edit.deleted or edit.seco or edit.cyclo or edit.homo or edit.replaced):
        return None
    cut = {frozenset(bond) for bond in edit.seco} | bridged
    rings = sum(
        1 for face in reference.faces
        if image.issuperset(face) and not any(frozenset(pair) in cut for pair in zip(face, face[1:] + face[:1]))
    )
    return edit, rings


def _name(reference, edit, rings):
    """(cost, key, name) under the numbering of the cage that gives the lowest locants: 'nor' first, then
    'homo', then the 'seco' and 'cyclo' locants together and in the order of citation, then the replaced atoms and the
    indicated hydrogen; None for a fragment that has too few atoms or rings to keep the fullerene name."""
    atoms = reference.graph.number_of_nodes() - len(edit.deleted)
    if atoms < reference.least_atoms or rings < reference.least_rings:
        return None
    best = None
    for numbering in reference.numberings:
        def pair(bond):
            return tuple(sorted(numbering[g] for g in bond))

        nor = sorted(numbering[g] for g in edit.deleted)
        homo = sorted(pair(bond) for bond, _ in edit.homo)
        cyclo = sorted(pair(bond) for bond in edit.cyclo)
        seco = sorted(pair(bond) for bond in edit.seco)
        joined = sorted(n for p in cyclo + seco for n in p)
        cited = [n for p in cyclo + seco for n in p]
        replaced = sorted(numbering[g] for g in edit.replaced)
        indicated = sorted(pair(s) + (1,) if isinstance(s, tuple) else (numbering[s], 0, 0) for s in edit.saturated)
        key = (nor, homo, joined, cited, replaced, indicated)
        if best is None or key < best[0]:
            best = (key, numbering)
    key, numbering = best
    return edit.cost, key, _text(reference, edit, numbering)


def _plural(count, word):
    return (multiplying_prefix(count) if count > 1 else "") + word


def _text(reference, edit, numbering):
    def pair(bond):
        return tuple(sorted(numbering[g] for g in bond))

    replaced_homo = {pair(bond)[0]: z for bond, z in edit.homo if z != 6}
    replaced = {}
    for g, z in edit.replaced.items():
        replaced.setdefault(z, []).append(str(numbering[g]))
    for low, z in replaced_homo.items():
        replaced.setdefault(z, []).append(f"{low}a")
    nor = sorted(numbering[g] for g in edit.deleted)
    homo = sorted(pair(bond) for bond, _ in edit.homo)
    cyclo = sorted(pair(bond) for bond in edit.cyclo)
    seco = sorted(pair(bond) for bond in edit.seco)
    if sum(len(v) for v in replaced.values()) == 1 and reference.orbit_count == 1 and not (nor or homo or cyclo or seco or edit.saturated):
        return _REPLACEMENT_PREFIXES[next(iter(replaced))] + reference.stem
    parts = []
    if edit.saturated:
        texts = []
        for s in edit.saturated:
            if isinstance(s, tuple):
                low, high = pair(s)
                texts.append(f"{low}aH" if low in replaced_homo else f"{low}({high})aH")
            else:
                texts.append(f"{numbering[s]}H")
        parts.append(",".join(texts))
    for z in sorted(replaced, key=_SENIORITY.get):
        locants = sorted(replaced[z], key=lambda v: (int(v.rstrip("abc")), v))
        parts.append(f"{','.join(locants)}-{_plural(len(locants), _REPLACEMENT_PREFIXES[z])}")
    for pairs, word in ((cyclo, "cyclo"), (seco, "seco")):
        if pairs:
            parts.append(":".join(f"{a},{b}" for a, b in pairs) + "-" + _plural(len(pairs), word))
    if homo:
        parts.append(",".join(f"{low}({high})a" for low, high in homo) + "-" + _plural(len(homo), "homo"))
    if nor:
        parts.append(",".join(str(n) for n in nor) + "-" + _plural(len(nor), "nor"))
    return "-".join(parts) + reference.stem
