"""Ring assemblies of identical von Baeyer components joined by single bonds (P-28.2, P-28.3, P-28.4, P-31.1.7).

Two components take primed locants (`2,2'-bi(bicyclo[2.2.2]octane)`), three to six the composite locants of the
unbranched chain; unsaturation is cited after the assembly bracket (P-31.1.7.1), skeletal replacement prefixes in
front of it (P-28.4.1, P-31.1.7.3). Locants go in order to junctions, heteroatoms, then multiple bonds."""

from itertools import product

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    multiplied_word,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    superscript_locant,
    von_baeyer_bond_citation,
)
from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._numerals import numerical_term
from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
from ._substituents import format_substituent_prefixes, name_branch

_MIN_COMPONENTS, _MAX_COMPONENTS = 2, 6
_MULTIPLIER = {3: "ter", 4: "quater", 5: "quinque", 6: "sexi"}
_MAX_COMBINATIONS = 200_000

# P-28.4.2: seniority order of the skeletal replacement prefixes
_REPLACEMENT = (
    (8, "O", "oxa", 2), (16, "S", "thia", 2), (34, "Se", "selena", 2), (52, "Te", "tellura", 2), (7, "N", "aza", 3),
    (15, "P", "phospha", 3), (33, "As", "arsa", 3), (51, "Sb", "stiba", 3), (83, "Bi", "bisma", 3),
    (14, "Si", "sila", 4), (32, "Ge", "germa", 4), (50, "Sn", "stanna", 4), (82, "Pb", "plumba", 4),
    (5, "B", "bora", 3), (13, "Al", "alumina", 3), (31, "Ga", "galla", 3), (49, "In", "inda", 3), (81, "Tl", "thalla", 3),
)
_PREFIX = {z: prefix for z, _, prefix, _ in _REPLACEMENT}
_SENIORITY = {z: rank for rank, (z, _, _, _) in enumerate(_REPLACEMENT)}
_VALENCE = {z: valence for z, _, _, valence in _REPLACEMENT}


def _components(mol):
    """Atom sets of the ring systems (ring atoms joined through ring bonds) with their ring bonds."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    ring_bonds = [b for b in mol.GetBonds() if b.IsInRing()]
    for bond in ring_bonds:
        parent[find(bond.GetBeginAtomIdx())] = find(bond.GetEndAtomIdx())
    groups = {}
    for bond in ring_bonds:
        root = find(bond.GetBeginAtomIdx())
        atoms, bonds = groups.setdefault(root, (set(), []))
        atoms.update((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        bonds.append(bond)
    return list(groups.values())


def _chain_order(mol, comps):
    """Component indices in one path order and the (atom, atom) junction of each consecutive pair, else None."""
    of_atom = {a: i for i, (atoms, _) in enumerate(comps) for a in atoms}
    links = {i: [] for i in range(len(comps))}
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in of_atom and b in of_atom and of_atom[a] != of_atom[b]:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return None
            links[of_atom[a]].append((of_atom[b], a, b))
            links[of_atom[b]].append((of_atom[a], b, a))
    n = len(comps)
    if sum(len(v) for v in links.values()) != 2 * (n - 1) or any(len({x[0] for x in v}) != len(v) for v in links.values()):
        return None
    ends = [i for i in links if len(links[i]) == 1]
    if len(ends) != 2:
        return None
    order, previous = [ends[0]], None
    while len(order) < n:
        step = [x for x in links[order[-1]] if x[0] != previous]
        if len(step) != 1:
            return None
        previous = order[-1]
        order.append(step[0][0])
    if order[-1] != ends[1]:
        return None
    junctions = []
    for i in range(n - 1):
        entry = next(x for x in links[order[i]] if x[0] == order[i + 1])
        junctions.append((entry[1], entry[2]))
    return order, junctions


def _skeleton(atoms, bonds):
    """Saturated all-carbon copy of one ring system and the map from its atom indices back."""
    ordered = sorted(atoms)
    index = {a: i for i, a in enumerate(ordered)}
    skeleton = Chem.RWMol()
    for _ in ordered:
        skeleton.AddAtom(Chem.Atom(6))
    for bond in bonds:
        skeleton.AddBond(index[bond.GetBeginAtomIdx()], index[bond.GetEndAtomIdx()], Chem.BondType.SINGLE)
    result = skeleton.GetMol()
    result.UpdatePropertyCache(strict=False)
    return result, ordered


def _von_baeyer_numberings(skeleton, ring_count):
    """(parent name, [atom order]) of every numbering of the lowest descriptor class, else None."""
    if ring_count == 2:
        core = find_bicyclic_core(skeleton)
        if core is None or any(not bridge for bridge in core[2]):
            return None
        return bicyclic_parent_name(core), [order for order in iter_bicyclic_numberings(core)]
    core = find_polycyclic_core(skeleton, ring_count)
    if core is None:
        return None
    found = list(iter_polycyclic_candidates(core, ring_count))
    if not found:
        return None
    best = min(key for _, _, key in found)
    names = {parent for _, parent, key in found if key == best}
    if len(names) != 1:
        return None
    return names.pop(), [order for order, _, key in found if key == best]


def _hetero_ok(mol, atom):
    z = atom.GetAtomicNum()
    return (
        z in _PREFIX
        and not atom.GetFormalCharge()
        and not atom.GetNumRadicalElectrons()
        and not atom.GetIsotope()
        and atom.GetTotalValence() == _VALENCE[z]
    )


def find_vb_ring_assembly_core(mol):
    """(path, junctions, parent, per-component numberings) for an unbranched chain of two to six identical von Baeyer
    carbocyclic or skeletally replaced ring systems whose only substituents are halogens, else None."""
    if mol.GetRingInfo().NumRings() < 4 or specified_stereocenters(mol) is not None:
        return None
    comps = _components(mol)
    if not _MIN_COMPONENTS <= len(comps) <= _MAX_COMPONENTS:
        return None
    ring_atoms = set().union(*(atoms for atoms, _ in comps))
    halogens = halogen_substituents(mol)
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() or (atom.GetIdx() not in ring_atoms and atom.GetIdx() not in halogens):
            return None
        if atom.GetIdx() in ring_atoms and atom.GetAtomicNum() != 6 and not _hetero_ok(mol, atom):
            return None
    chain = _chain_order(mol, comps)
    if chain is None:
        return None
    order, junctions = chain
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for pair in junctions for a in pair):
        return None

    parents, numberings = set(), []
    for i in order:
        atoms, bonds = comps[i]
        ring_count = len(bonds) - len(atoms) + 1
        if ring_count < 2:
            return None
        skeleton, ordered = _skeleton(atoms, bonds)
        found = _von_baeyer_numberings(skeleton, ring_count)
        if found is None:
            return None
        parent, orders = found
        parents.add((parent, Chem.MolToSmiles(skeleton)))
        numberings.append([[ordered[j] for j in o] for o in orders])
    if len(parents) != 1:
        return None
    parent, skeleton_smiles = next(iter(parents))
    from .core import smiles_to_iupac

    try:
        if smiles_to_iupac(skeleton_smiles) != parent:
            return None
    except Exception:
        return None
    return [comps[i][0] for i in order], junctions, parent, numberings


def _locant_text(count, ring, position):
    return f"{position}{chr(39) * ring}" if count == 2 else superscript_locant(ring + 1, position)


def _unsaturation_ending(enes, ynes):
    if not enes and not ynes:
        return ""
    ene_word, yne_word = multiplied_word(len(enes), "ene"), multiplied_word(len(ynes), "yne")
    cited = lambda found: ",".join(found)
    if enes and ynes:
        return f"{cited(enes)}-{ene_word[:-1]}-{cited(ynes)}-{yne_word}"
    return f"{cited(enes)}-{ene_word}" if enes else f"{cited(ynes)}-{yne_word}"


def _hetero_text(mol, locations, count):
    by_element = {}
    for atom, (ring, position) in locations.items():
        z = mol.GetAtomWithIdx(atom).GetAtomicNum()
        if z in _PREFIX:
            by_element.setdefault(z, []).append((ring, position))
    pieces = []
    for z in sorted(by_element, key=_SENIORITY.get):
        cited = sorted(by_element[z], key=lambda rp: _sort_key(count, *rp))
        multiplier = numerical_term(len(cited)) if len(cited) > 1 else ""
        pieces.append(f"{','.join(_locant_text(count, *rp) for rp in cited)}-{multiplier}{_PREFIX[z]}")
    return "-".join(pieces)


def _sort_key(count, ring, position):
    return (position, ring) if count == 2 else (ring, position)


def name_vb_ring_assembly(mol, core) -> str:
    path, junctions, parent, numberings = core
    count = len(path)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    sort_key = lambda ring, position: _sort_key(count, ring, position)
    ring_word = f"bi({parent})" if count == 2 else f"{_MULTIPLIER[count]}{parent}"

    best_key, best_name = None, None
    for reverse in (False, True):
        rings = list(reversed(path)) if reverse else list(path)
        links = [(b, a) for a, b in reversed(junctions)] if reverse else list(junctions)
        orders = list(reversed(numberings)) if reverse else numberings
        candidates = []
        for ring, (atoms, ring_numberings) in enumerate(zip(rings, orders)):
            toward = [None, None]
            if ring > 0:
                toward[0] = links[ring - 1][1]
            if ring < count - 1:
                toward[1] = links[ring][0]
            seen, options = set(), []
            for numbering in ring_numberings:
                position = {atom: i + 1 for i, atom in enumerate(numbering)}
                bonds = []
                for bond in mol.GetBonds():
                    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    if a in atoms and b in atoms and bond.GetBondTypeAsDouble() != 1.0:
                        primary, display, compound = von_baeyer_bond_citation(position, a, b)
                        high = max(position[a], position[b])
                        bonds.append((bond.GetBondTypeAsDouble(), primary, high if compound else None))
                hetero = tuple(sorted((mol.GetAtomWithIdx(a).GetAtomicNum(), position[a]) for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6))
                substituents = []
                for atom in atoms:
                    roots = [nb for nb in graph[atom] if nb not in atoms and nb in halogens]
                    for root in roots:
                        substituents.append((position[atom], name_branch(graph, root, atom, halogens, mol=mol)))
                attach = tuple(position[a] if a is not None else None for a in toward)
                signature = (attach, hetero, tuple(sorted(bonds, key=str)), tuple(sorted(substituents)))
                if signature in seen:
                    continue
                seen.add(signature)
                options.append((position, bonds, substituents))
            candidates.append(options)
        size = 1
        for options in candidates:
            size *= len(options)
        if size > _MAX_COMBINATIONS:
            raise UnsupportedStructure("too many equivalent numberings of the components of this ring assembly")
        for combo in product(*candidates):
            locations = {}
            for ring, (position, _, _) in enumerate(combo):
                for atom, p in position.items():
                    locations[atom] = (ring, p)
            pairs = [(locations[a], locations[b]) for a, b in links]
            if count == 2:
                pairs = [tuple(sorted(pairs[0], key=lambda rp: sort_key(*rp)))]
            junction_set = tuple(sorted((loc for pair in pairs for loc in pair), key=lambda rp: sort_key(*rp)))
            junction_citation = tuple(sort_key(*loc) for pair in pairs for loc in pair)
            hetero_items = [(loc, mol.GetAtomWithIdx(a).GetAtomicNum()) for a, loc in locations.items() if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
            hetero_set = tuple(sorted(sort_key(*loc) for loc, _ in hetero_items))
            hetero_by_element = tuple(
                tuple(sorted(sort_key(*loc) for loc, z in hetero_items if z == element)) for element, *_ in _REPLACEMENT
            )
            enes, ynes, compound_count, primaries, fulls = [], [], 0, [], []
            for ring, (_, bonds, _) in enumerate(combo):
                for order, primary, high in bonds:
                    text = _locant_text(count, ring, primary)
                    primaries.append(sort_key(ring, primary))
                    fulls.append(sort_key(ring, primary))
                    if high is not None:
                        text += f"({_locant_text(count, ring, high)})"
                        compound_count += 1
                        fulls.append(sort_key(ring, high))
                    (enes if order == 2.0 else ynes).append((sort_key(ring, primary), text))
            ene_texts = [t for _, t in sorted(enes)]
            yne_texts = [t for _, t in sorted(ynes)]
            unsaturation = (compound_count, tuple(sorted(primaries)), tuple(sorted(fulls)))
            grouped_input = {}
            for ring, (_, _, substituents) in enumerate(combo):
                for p, text in substituents:
                    grouped_input.setdefault(_locant_text(count, ring, p), []).append(text)
            grouped = group_substituents(grouped_input)
            sub_set, _, sub_citation = substituent_locant_set_and_citation(grouped)
            prefix = format_substituent_prefixes(grouped) if grouped else ""

            junction_text = ":".join(f"{_locant_text(count, *a)},{_locant_text(count, *b)}" for a, b in pairs)
            base = f"{junction_text}-{ring_word}"
            ending = _unsaturation_ending(ene_texts, yne_texts)
            hetero_prefix = _hetero_text(mol, locations, count)
            if ending:
                body = f"{hetero_prefix}[{base}]-{ending}"
            else:
                body = f"{hetero_prefix}-{base}" if hetero_prefix else base
            name = f"{prefix}-{body}" if prefix else body
            key = (junction_set, junction_citation, hetero_set, hetero_by_element, unsaturation, sub_set, sub_citation, name)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
