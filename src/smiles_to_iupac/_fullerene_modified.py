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
UNMODIFIED = object()


class _Reference:
    def __init__(self, smiles, stem, atoms, rings):
        mol = Chem.MolFromSmiles(smiles)
        self.stem = stem
        self.least_atoms = atoms
        self.least_rings = rings
        self.graph = _cage_graph(mol, range(mol.GetNumAtoms()))
        self.adjacent = {a: frozenset(self.graph[a]) for a in self.graph}
        self.faces = [tuple(face) for face in _Cage(self.graph).faces]
        self.signatures = {
            (sum(len(f) == 5 for f in self.faces if a in f), sum(len(f) == 6 for f in self.faces if a in f)) for a in self.graph
        }
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
    fusions: list = field(default_factory=list)
    mapping: dict = field(default_factory=dict)
    homo_atoms: list = field(default_factory=list)

    @property
    def cost(self):
        return len(self.deleted) + len(self.seco) + len(self.cyclo) + len(self.homo)


@dataclass
class Candidate:
    """A numbered modified cage: `fusions` are (cage atom, cage atom, tag) of the rings fused to it."""

    cost: int
    key: tuple
    reference: _Reference
    edit: _Edit
    numbering: dict
    numberings: list = field(default_factory=list)
    ties: list = field(default_factory=list)

    def text(self, fused=None):
        return _text(self.reference, self.edit, self.numbering, fused)

    def atom_locants(self, numbering):
        """The locant of every atom of the cage under `numbering`; a homo atom takes the lower locant of its bond and 'a'."""
        locant_of = {a: str(numbering[g]) for a, g in self.edit.mapping.items()}
        for x, (g, h) in self.edit.homo_atoms:
            locant_of[x] = f"{min(numbering[g], numbering[h])}a"
        return locant_of


def _skeleton(mol, atoms=None, fusion_atoms=frozenset(), bare=frozenset()):
    """(graph, hydrogens) of the heavy atoms `atoms` when they can form a modified cage, else None. A fusion atom takes
    the place of a hydrogen atom with the bond to its fused ring (P-27.6.1); a `bare` atom is judged as if the bonds
    left out of `atoms` were hydrogen atoms (a spiro atom, P-27.6.3)."""
    atoms = set(range(mol.GetNumAtoms())) if atoms is None else set(atoms)
    graph = nx.Graph()
    graph.add_nodes_from(atoms)
    graph.add_edges_from(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
    )
    if not 31 <= len(atoms) <= 90 or not nx.is_connected(graph):
        return None
    bridges = {frozenset(bridge) for bridge in nx.bridges(graph)}
    if graph.number_of_edges() - len(atoms) + 1 < 10:
        return None
    if sum(all(frozenset((a, n)) in bridges for n in graph[a]) for a in atoms) > _MAX_HOMO:
        return None
    for a in atoms:
        atom = mol.GetAtomWithIdx(a)
        if (
            atom.GetAtomicNum() != 6 and atom.GetAtomicNum() not in _REPLACEMENT_PREFIXES
            or atom.GetFormalCharge()
            or atom.GetNumRadicalElectrons()
            or atom.GetIsotope()
            or graph.degree(a) > 3
            or any(b.GetBondTypeAsDouble() > 2.0 for b in atom.GetBonds() if b.GetOtherAtomIdx(a) in atoms)
        ):
            return None
    hydrogens = {
        a: 3 - graph.degree(a) if a in fusion_atoms else 4 - graph.degree(a) if a in bare else mol.GetAtomWithIdx(a).GetTotalNumHs()
        for a in atoms
    }
    return graph, hydrogens


def _embeddings(graph, reference, limit=_MAX_SOLUTIONS):
    """Up to `limit` injective maps of `graph` into the cage in which every bond of `graph` is a cage bond."""
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
        if len(found) >= limit:
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


def _unembeddable_regions(graph, reference):
    """The bond sets around an atom whose small rings no atom of the cage has (a ring of three or four atoms, two
    five-membered rings or three six-membered rings at an atom of C60). A bond formed by 'cyclo' lies in each of them."""
    faces = [cycle for cycle in nx.minimum_cycle_basis(graph) if len(cycle) <= 6]
    regions = [_ring_bonds(cycle) for cycle in faces if len(cycle) < 5]
    for atom in graph:
        at = [cycle for cycle in faces if atom in cycle and len(cycle) >= 5]
        counts = (sum(len(c) == 5 for c in at), sum(len(c) == 6 for c in at))
        if not any(counts[0] <= pentagons and counts[1] <= hexagons for pentagons, hexagons in reference.signatures):
            regions.append(frozenset().union(*(_ring_bonds(cycle) for cycle in at)))
    return regions


def _ring_bonds(cycle):
    return frozenset(frozenset(pair) for pair in zip(cycle, cycle[1:] + cycle[:1]))


def _cyclo_choices(graph, forced, count, regions):
    """The connected graphs left when `count` bonds that are not `forced` and lie in every unembeddable region are taken
    out of `graph`."""
    bonds = [frozenset(bond) for bond in graph.edges if frozenset(bond) not in forced]
    for removed in combinations(bonds, count):
        if any(region.isdisjoint(removed) for region in regions):
            continue
        rest = graph.copy()
        rest.remove_edges_from(tuple(bond) for bond in removed)
        if nx.is_connected(rest):
            yield rest


def name_modified_fullerene(mol):
    """Name of `mol` as a nor-, seco-, homo-, cyclo- or heteroatom-replaced C60-Ih or C70-D5h(6) fullerene, else None. A
    fragment keeps the fullerene name only with more than half of the atoms and a third of the rings of the cage
    (P-52.2.6.1)."""
    found = modified_candidate(mol)
    return None if found is None else found.text()


def modified_candidate(mol, atoms=None, fusions=(), fusion_atoms=frozenset(), extra_key=None, bare=frozenset()):
    """The numbered cage of the lowest cost that the heavy atoms `atoms` form, with `fusions` (atom, atom, (kind, label,
    index)) as the atom pairs that rings are fused to (kind 0) or bridges are attached to (kind 1), else None.
    `extra_key(numbering, edit)` ranks the numberings that tie on all other locants."""
    skeleton = _skeleton(mol, atoms, fusion_atoms, bare)
    if skeleton is None:
        return None
    graph, hydrogens = skeleton
    element = {a: mol.GetAtomWithIdx(a).GetAtomicNum() for a in graph}
    if any(element[a] == 6 and graph.degree(a) + hydrogens[a] not in (3, 4) for a in graph):
        return None
    homo_candidates = [
        x for x in graph
        if graph.degree(x) == 2
        and not graph.has_edge(*graph[x])
        and all(graph.degree(n) == 3 for n in graph[x])
        and (hydrogens[x] == 2 or element[x] != 6)
    ]
    saturated = [a for a in graph if hydrogens[a] and (graph.degree(a) + hydrogens[a] == 4 or element[a] != 6)]
    if len(saturated) > 1 + len(set(saturated) & set(homo_candidates)):
        return None
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
    best, ties = None, []
    for lower, _, reference, homo, reduced, forced in sorted(attempts, key=lambda a: a[:2]):
        if best is not None and lower > best.cost:
            break
        deficiency = sum(3 - d for _, d in reduced.degree)
        deleted = reference.graph.number_of_nodes() - reduced.number_of_nodes()
        regions = None
        for cyclo_count in range(_MAX_CYCLO + 1):
            cleaved = max(0, -(-(deficiency + 2 * cyclo_count - 3 * deleted) // 2))
            if best is not None and lower + cyclo_count + cleaved > best.cost:
                continue
            if cyclo_count and regions is None:
                regions = _unembeddable_regions(reduced, reference)
            if cyclo_count > 1 and not regions:
                continue
            for rest in _cyclo_choices(reduced, forced, cyclo_count, regions or []):
                for mapping in _embeddings(rest, reference):
                    found = _edit(graph, hydrogens, element, reference, mapping, homo, fusions, bare)
                    if found is UNMODIFIED:
                        return None
                    if found is None:
                        continue
                    candidate = _name(reference, *found, extra_key)
                    if candidate is None:
                        continue
                    if best is None or (candidate.cost, candidate.key) < (best.cost, best.key):
                        best, ties = candidate, [candidate]
                    elif (candidate.cost, candidate.key) == (best.cost, best.key):
                        ties.append(candidate)
    if best is not None:
        best.ties = ties
    return best


def _edit(graph, hydrogens, element, reference, mapping, homo, fusions=(), bare=frozenset()):
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
        edit.homo_atoms.append((x, bond))
        if hydrogens[x] and x not in bare:
            edit.saturated.append(bond)
    edit.fusions = [(mapping[a], mapping[b], tag) for a, b, tag in fusions]
    edit.mapping = dict(mapping)
    if not (edit.deleted or edit.seco or edit.cyclo or edit.homo or edit.replaced or edit.fusions):
        return UNMODIFIED
    if sum(1 for s in edit.saturated if not isinstance(s, tuple)) > 1:
        return None
    cut = {frozenset(bond) for bond in edit.seco} | bridged
    rings = sum(
        1 for face in reference.faces
        if image.issuperset(face) and not any(frozenset(pair) in cut for pair in zip(face, face[1:] + face[:1]))
    )
    return edit, rings


def _name(reference, edit, rings, extra_key=None):
    """The candidate under the numbering of the cage that gives the lowest locants: 'nor' first, then 'homo', then
    the 'seco' and 'cyclo' locants together and in the order of citation, then the fusion locants, then the replaced
    atoms and the indicated hydrogen; None for a fragment that has too few atoms or rings to keep the fullerene
    name."""
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
        fused = sorted((tag[1], pair((first, second))) for first, second, tag in edit.fusions if tag[0] == 0)
        fused_joined = sorted(n for _, p in fused for n in p)
        fused_cited = [n for _, p in fused for n in p]
        replaced = sorted(numbering[g] for g in edit.replaced)
        indicated = sorted(pair(s) + (1,) if isinstance(s, tuple) else (numbering[s], 0, 0) for s in edit.saturated)
        attached = extra_key(numbering, edit) if extra_key else ()
        key = (nor, homo, joined, cited, fused_joined, fused_cited, attached, replaced, indicated)
        if best is None or key < best[0]:
            best, tied = (key, numbering), [numbering]
        elif key == best[0]:
            tied.append(numbering)
    key, numbering = best
    return Candidate(edit.cost, key, reference, edit, numbering, tied)


def _plural(count, word):
    return (multiplying_prefix(count) if count > 1 else "") + word


def _text(reference, edit, numbering, fused=None):
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
    stem = reference.stem
    parts, indicated = [], []
    if sum(len(v) for v in replaced.values()) == 1 and reference.orbit_count == 1 and not (nor or homo or cyclo or seco or edit.saturated):
        stem = _REPLACEMENT_PREFIXES[next(iter(replaced))] + stem
    else:
        for s in edit.saturated:
            if isinstance(s, tuple):
                low, high = pair(s)
                indicated.append(((low, 1, high), f"{low}aH" if low in replaced_homo else f"{low}({high})aH"))
            else:
                indicated.append(((numbering[s], 0, 0), f"{numbering[s]}H"))
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
    rest = "-".join(parts)
    if fused is None:
        return "-".join(([",".join(text for _, text in sorted(indicated))] if indicated else []) + parts) + stem
    hydro, ring_indicated, prefix, substituents = fused(numbering, edit)
    indicated = sorted(indicated + ring_indicated)
    front = hydro + ("-" if hydro and indicated else "") + (",".join(text for _, text in indicated) + "-" if indicated else "")
    after = front + prefix
    joiner = "-" if substituents and after[:1].isdigit() else ""
    return substituents + joiner + after + ("-" + rest if rest else "") + stem
