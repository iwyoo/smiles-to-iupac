"""Polyesters of one polyol (P-65.6.3.3.3): the polyol's best skeleton -- a ring system or an acyclic chain
carrying the most ester oxygens -- is cited as a multivalent group ('cyclohexane-1,2-diyl', 'naphthalene-2,3-diyl',
'butane-1,4-diyl') before the anions; other esters become acyloxy prefixes (P-65.6.3.3.4.2). Numbering follows
P-31.1.4: heteroatoms, indicated hydrogen, free valences, hydro/ene, all prefixes, citation order, anion locants,
CIP descriptors. Two esters on a symmetric group cite no acid locants (P-65.6.3.3.3.2).
"""

import contextvars
import re

from rdkit import Chem

from . import _aromatic
from ._free_valence import SUFFIX_OF_ORDER, attach, suffix_of
from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    ring_cycle,
    substituent_locant_set_and_citation,
)
from ._fullerene import is_fullerene_cage
from ._fullerene_numbering import cage_numberings
from ._diester_anions import acid_anions, anion_locant_key, cip_labels, cite_anions
from ._functional_prefixes import functional_names, nitro_atoms
from ._ring_diyl_numbering import ANION_SUFFIX, SUFFIX_ATOMS, _locs, _yl, chain_numberings, monocycle_numberings, system_numberings
from ._substituents import format_substituent_prefixes, name_branch

_DESCRIPTOR_ORDER = {"R": 0, "S": 1, "r": 2, "s": 3}
PARENT_START = contextvars.ContextVar("parent_start", default=0)


def _isotope_key(position_of, skeleton):
    """Locants of the modified atoms in increasing order, then the higher nuclide first (P-82.5.2)."""
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if not context:
        return ()
    table = Chem.GetPeriodicTable()
    found = []
    for atom, entry in context["labels"].items():
        if atom not in skeleton or atom not in position_of:
            continue
        nuclides = [(-table.GetAtomicNumber("".join(c for c in entry["skeleton"] if c.isalpha())), -int("".join(c for c in entry["skeleton"] if c.isdigit())))] if entry["skeleton"] else []
        nuclides += [(-1, -{"1H": 1, "2H": 2, "3H": 3}[n]) for n in entry["H"]]
        found.append((position_of[atom], tuple(sorted(nuclides))))
    found.sort()
    return (tuple(loc for loc, _ in found), tuple(n for _, n in found))


def _component(graph, start, blocked):
    seen = set()
    stack = [start]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(n for n in graph[node] if n not in blocked and n not in seen)
    return seen


def _diyl_is_symmetric(mol, keep, ester_oxygens):
    rw = Chem.RWMol(mol)
    kept = sorted(keep | set(ester_oxygens))
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(kept), reverse=True):
        rw.RemoveAtom(idx)
    sub = rw.GetMol()
    Chem.SanitizeMol(sub, catchErrors=True)
    ranks = list(Chem.CanonicalRankAtoms(sub, breakTies=False, includeChirality=False))
    position = {idx: i for i, idx in enumerate(kept)}
    return len({ranks[position[o]] for o in ester_oxygens}) == 1


def _system_of(mol, atom):
    ring_info = mol.GetRingInfo()
    rings = [set(r) for r in ring_info.AtomRings()]
    start = next((r for r in rings if atom in r), None)
    if start is None:
        return None
    members = [start]
    atoms = set(start)
    grew = True
    while grew:
        grew = False
        for r in rings:
            if r not in members and r & atoms:
                members.append(r)
                atoms |= r
                grew = True
    return [tuple(r) for r in ring_info.AtomRings() if set(r) in members], atoms


def _chain_paths(graph, component, alcohol_atoms):
    sub = {a: [n for n in graph[a] if n in component] for a in component}
    best = (0, 0)
    paths = []
    for u in component:
        parent = {u: None}
        queue = [u]
        while queue:
            node = queue.pop(0)
            for n in sub[node]:
                if n not in parent:
                    parent[n] = node
                    queue.append(n)
        for v in component:
            if v not in parent or (v < u):
                continue
            path = [v]
            while path[-1] != u:
                path.append(parent[path[-1]])
            score = (sum(1 for a in path if a in alcohol_atoms), len(path))
            if score > best:
                best, paths = score, []
            if score == best:
                paths.append(path)
    return best, paths


def select_skeleton(mol, graph, matches):
    """(kind, rings_or_paths, atoms) of the single best ring system or chain, or None when no skeleton carries
    an ester oxygen; raises when several equally ranked units tie (multiplicative/ring-assembly names)."""
    alcohol = {m[3].GetIdx() for m in matches}
    ring_info = mol.GetRingInfo()
    units = []
    seen = set()
    for idx in sorted(alcohol):
        if ring_info.NumAtomRings(idx):
            rings, atoms = _system_of(mol, idx)
            if frozenset(atoms) not in seen:
                seen.add(frozenset(atoms))
                units.append(("ring", len([a for a in alcohol if a in atoms]), rings, atoms))
    chain_atoms = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not ring_info.NumAtomRings(a.GetIdx())
    }
    done = set()
    for idx in sorted(alcohol):
        if idx not in chain_atoms or idx in done:
            continue
        component = set()
        stack = [idx]
        while stack:
            node = stack.pop()
            if node not in component:
                component.add(node)
                stack.extend(n for n in graph[node] if n in chain_atoms and n not in component)
        done |= component
        (valence, _), paths = _chain_paths(graph, component, alcohol)
        units.append(("chain", valence, paths, set().union(*(set(p) for p in paths))))
    if not units:
        return None
    top = max(u[1] for u in units)
    best = [u for u in units if u[1] == top]
    if any(u[0] == "ring" for u in best):
        best = [u for u in best if u[0] == "ring"]
    if len(best) > 1:
        raise UnsupportedStructure(
            "several equally ranked ring/chain units carry the esters (multiplicative or ring-assembly names, "
            "P-15.3/P-28) -- not handled here"
        )
    kind, _, body, atoms = best[0]
    return kind, body, atoms


def name_diester_ring_diyl(mol, matches):
    graph = adjacency(mol)
    nitro = nitro_atoms(mol)
    if any((a.GetFormalCharge() != 0 and a.GetIdx() not in nitro) or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")

    selection = select_skeleton(mol, graph, matches)
    if selection is None:
        raise UnsupportedStructure("no ring system or chain carries the ester oxygens")
    kind, body, pool = selection
    options = [(kind, [path], set(path)) for path in body] if kind == "chain" else [(kind, body, pool)]

    best = None
    for option_kind, option_body, option_pool in options:
        matches_on = [m for m in matches if m[3].GetIdx() in option_pool]
        result = _best_for_option(mol, graph, option_kind, option_body, option_pool, matches_on)
        if result is not None and (best is None or result[0] < best[0]):
            best = result
    if best is None:
        raise UnsupportedStructure("no admissible numbering places every attachment and substituent on this skeleton")
    return best[1]


def _best_for_option(mol, graph, kind, body, pool, matches_on):
    alcohol_idxs = [m[3].GetIdx() for m in matches_on]
    ester_oxygens = [m[2].GetIdx() for m in matches_on]
    anions = acid_anions(mol, matches_on)
    found = evaluate_skeleton(
        mol, graph, kind, body, pool, alcohol_idxs, set(ester_oxygens), "yl", anions=anions, matches_on=matches_on
    )
    if found is None:
        return None
    key, group_name, position_of, ring_stereo, side = found
    valence = len(matches_on)
    locants = [position_of[a] for a in alcohol_idxs]
    names = {name for name, _ in anions}
    cite = len(names) > 1 and not (valence == 2 and not ring_stereo and _diyl_is_symmetric(mol, side, ester_oxygens))
    return key, f"{group_name} {cite_anions(anions, locants, cite)}"


def _mixed_valence_text(diyl, free, valence, position_of, orders):
    """'cyclohexan-1-yl-2-ylidene' from the diyl text 'cyclohexane-1,2-diyl': valences cited by increasing bond
    order (P-29.3.2.2)."""
    tail = f"-{_locs(free)}-{_yl(valence)}"
    if not diyl.endswith(tail):
        return None
    by_order = {}
    for atom, order in orders.items():
        by_order.setdefault(int(order), []).append(position_of[atom])
    if any(order not in SUFFIX_OF_ORDER for order in by_order):
        return None
    return attach(diyl[: -len(tail)], {order: sorted(locs) for order, locs in by_order.items()})


CATION_CENTRES = contextvars.ContextVar("cation_centres", default=())


def evaluate_skeleton(
    mol, graph, kind, body, pool, attach, blocked, suffix, anions=None, matches_on=(), orders=None, centers=(),
    key_centers=(), n_names=(),
):
    from ._substituents import BRANCH_STEREO

    if not key_centers and suffix not in ("ide", "uide"):
        key_centers = tuple(c for c in CATION_CENTRES.get() if c[0] in pool)

    token = SUFFIX_ATOMS.set(frozenset(attach))
    stereo_token = BRANCH_STEREO.set({"atoms": {}, "bonds": {}, "used": set()}) if BRANCH_STEREO.get() is None else None
    anion_token = ANION_SUFFIX.set(
        frozenset(attach) | frozenset(a for a, _ in key_centers) if suffix in ("ide", "uide") else frozenset()
    )
    if not centers and not key_centers:
        centers = [c for c in _marked_centers(mol, pool) if c[0] not in attach or suffix in ("yl", "ylidene", "ylidyne")]
    try:
        return _evaluate_skeleton(
            mol, graph, kind, body, pool, attach, blocked, suffix, anions, matches_on, orders, centers, key_centers,
            n_names,
        )
    finally:
        SUFFIX_ATOMS.reset(token)
        if stereo_token is not None:
            BRANCH_STEREO.reset(stereo_token)
        ANION_SUFFIX.reset(anion_token)


def _evaluate_skeleton(
    mol, graph, kind, body, pool, attach, blocked, suffix, anions=None, matches_on=(), orders=None, centers=(),
    key_centers=(), n_names=(),
):
    """Best numbering of a ring system or chain with free valences/suffix at `attach`; returns
    (key, group_name, position_of, ring_stereo, side) or None."""
    valence = len(attach)
    mixed = bool(orders) and len(set(orders.values())) > 1
    halogens = halogen_substituents(mol)
    carbon_graph = carbon_adjacency(mol)

    if kind == "chain":
        numberings = chain_numberings(mol, body, valence)
    elif is_fullerene_cage(mol, pool):
        numberings = cage_numberings(mol, pool, attach)
    elif len(body) == 1:
        numberings = monocycle_numberings(mol, ring_cycle(graph, list(body[0])), set(), valence)
    else:
        numberings = system_numberings(mol, graph, body, pool)

    from ._substituents import BRANCH_STEREO

    stereo_context = BRANCH_STEREO.get()
    cache = {}

    def analyze(skeleton):
        key = frozenset(skeleton)
        if key in cache:
            return cache[key]
        side = _component(graph, next(iter(skeleton)), blocked)
        for atom, code in cip_labels(mol, side):
            if atom not in skeleton:
                stereo_context["atoms"][atom] = code
        if any(a in side for m in matches_on for a in (m[0].GetIdx(), m[1].GetIdx())):
            raise UnsupportedStructure(
                "a lactone or macrocyclic diester is a heterocyclic pseudoketone (P-65.6.3.5), not an ester of a "
                "polyol, and is named by the heterocycle rules"
            )
        seeds = [(a, n) for a in skeleton for n in graph[a] if n not in skeleton and n not in blocked]
        named, shown, covered = functional_names(mol, graph, seeds, set(skeleton) | blocked, halogens)
        unsaturated = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
            for b in mol.GetBonds()
            if b.GetBondTypeAsDouble() != 1.0
            and b.GetBeginAtomIdx() in side
            and b.GetEndAtomIdx() in side
            and not (b.GetBeginAtomIdx() in skeleton and b.GetEndAtomIdx() in skeleton)
            and not b.GetIsAromatic()
            and b.GetBeginAtomIdx() not in covered
            and b.GetEndAtomIdx() not in covered
        ]

        def branch(root, atom):
            if root in named:
                return named[root]
            atoms = _component(graph, root, set(skeleton) | blocked)
            bonds = [b for b in unsaturated if b[0] in atoms and b[1] in atoms]
            if bonds and any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in atoms):
                raise UnsupportedStructure("an unsaturated substituent bearing other groups is not supported yet")
            if bonds:
                return _aromatic._branch_name(
                    graph, carbon_graph, root, atom, set(skeleton) | blocked, shown, unsaturated, mol=mol
                )
            return name_branch(graph, root, atom, shown, mol=mol)

        cache[key] = (side, branch)
        return cache[key]

    stereo_all = None
    candidates = []
    for numbering in numberings:
        position_of = numbering.position_of
        if any(a not in position_of for a in attach):
            continue
        skeleton = pool if kind == "ring" else set(position_of)
        side, branch = analyze(skeleton)
        if stereo_all is None:
            stereo_all = cip_labels(mol, side)
        substituents = {}
        for atom in position_of:
            roots = [n for n in graph[atom] if n not in skeleton and n not in blocked]
            if roots:
                substituents[position_of[atom]] = [branch(root, atom) for root in roots]
        grouped = group_substituents(substituents)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        free = tuple(sorted(position_of[a] for a in attach))
        cite = tuple(position_of[a] for a in sorted(attach, key=lambda a: (orders[a], position_of[a]))) if mixed else ()
        ring_stereo = [(a, c) for a, c in stereo_all if a in skeleton]
        if any(("atom", a) not in stereo_context["used"] for a, _ in stereo_all if a not in skeleton):
            raise UnsupportedStructure("a stereocenter on a substituent is not supported yet")
        acid_key = anion_locant_key(anions, [position_of[a] for a in attach]) if anions else ()
        stereo_key = tuple(
            _DESCRIPTOR_ORDER.get(code, 9) for _, code in sorted((position_of[a], c) for a, c in ring_stereo)
        )
        center_key = tuple(sorted(position_of[a] for a, _ in centers))
        if key_centers:
            center_key = (
                tuple(sorted(position_of[a] for a, _ in key_centers if a in position_of)),
                tuple(sorted(position_of[a] for a, word in key_centers if word == "uide" and a in position_of)),
            )
        key = (
            numbering.pre_key, center_key, free, cite, numbering.unsat_key, locant_set, citation, acid_key, stereo_key,
            _isotope_key(position_of, skeleton),
        )
        candidates.append((key, numbering, grouped, free, ring_stereo, side))
    if not candidates:
        return None
    key, numbering, grouped, free, ring_stereo, side = min(candidates, key=lambda c: c[0])
    position_of = numbering.position_of

    substituted = frozenset(loc for info in grouped.values() for loc in info["locants"])
    parent = numbering.text(free, valence, substituted, suffix)
    if mixed:
        parent = _mixed_valence_text(parent, free, valence, position_of, orders)
        if parent is None:
            return None
    if centers:
        parent = _with_anion_centers(parent, [(position_of[a], word) for a, word in centers])
    if n_names:
        from ._polyfunctional import _with_n_names

        grouped = _with_n_names(grouped, n_names, position_of, len(attach))
    prefixes = format_substituent_prefixes(grouped)
    if prefixes and (parent[0].isdigit() or parent[0] == "Δ"):
        prefixes += "-"
    group_name = prefixes + parent
    PARENT_START.set(len(prefixes) - (1 if prefixes.endswith("-") else 0))
    if ring_stereo:
        labels = ",".join(f"{loc}{code}" for loc, code in sorted((position_of[a], c) for a, c in ring_stereo))
        group_name = f"({labels})-{group_name}"
    return key, group_name, position_of, ring_stereo, side


def _marked_centers(mol, atoms):
    found = []
    for a in atoms:
        atom = mol.GetAtomWithIdx(a)
        if atom.HasProp("_anion_word"):
            found.extend([(a, atom.GetProp("_anion_word"))] * int(atom.GetProp("_anion_charge")))
        elif atom.HasProp("_anion"):
            found.extend([(a, "ide")] * int(atom.GetProp("_anion")))
    return found


def _with_anion_centers(parent, located):
    """Cite 'ide'/'uide' centers (as 'id'/'uid') before the free-valence suffix (P-72.6.3)."""
    words = {word for _, word in located}
    if len(words) != 1:
        raise UnsupportedStructure("ide and uide centers inside one substituent group are not supported yet")
    locants = sorted(loc for loc, _ in located)
    word = words.pop()
    count = len(locants)
    text = f"{','.join(str(x) for x in locants)}-{_COUNT_PREFIX[count]}{word[:-1]}"
    match = re.search(r"(-\d+(?:,\d+)*-)?(?:di|tri)?(?:yl|ylidene|ylidyne)$", parent)
    if match is None:
        raise UnsupportedStructure("this substituent group cannot carry an anionic center yet")
    start = match.start()
    base = parent[:start]
    if base.endswith("e"):
        base = base[:-1]
    return base + "-" + text + match.group(0)


_COUNT_PREFIX = {1: "", 2: "di", 3: "tri"}


def _assembly_group(mol, graph, root, parent, atoms):
    from ._polyfunctional import _arm_atoms, assembly_substituent

    arm = _arm_atoms(graph, root, parent)
    joined = any(
        a not in atoms and mol.GetAtomWithIdx(a).IsInRing() and any(n in atoms for n in graph[a]) for a in arm
    )
    if not joined:
        return None
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    return assembly_substituent(mol, graph, root, parent, halogen_substituents(mol), aromatic)


def ring_substituent_name(mol, graph, root, parent):
    """(name, is_compound) of the ring system entered at `root` from `parent`, as a substituent prefix."""
    from ._glycosyl import glycosyl_branch

    glycosyl = glycosyl_branch(mol, graph, root, parent)
    if glycosyl is not None:
        return glycosyl
    rings, atoms = _system_of(mol, root)
    order = mol.GetBondBetweenAtoms(parent, root).GetBondTypeAsDouble()
    if any(mol.GetAtomWithIdx(a).GetFormalCharge() for a in atoms):
        return _cationic_ring_substituent(mol, graph, root, parent, rings, atoms, order)
    assembly = _assembly_group(mol, graph, root, parent, atoms)
    if assembly is not None:
        return assembly
    suffix = suffix_of(order)
    if suffix is None:
        raise UnsupportedStructure("this ring substituent bond is not supported yet")
    found = evaluate_skeleton(mol, graph, "ring", rings, atoms, [root], {parent}, suffix)
    if found is None:
        raise UnsupportedStructure("this ring substituent has no supported name yet")
    from ._substituents import BRANCH_STEREO

    context = BRANCH_STEREO.get()
    if context:
        for atom, _ in found[3]:
            context["used"].add(("atom", atom))
    name = re.sub(r"(cyclo[a-z]+?)an-1-(yl|ylidene|ylidyne)$", r"\1\2", found[1])
    return name, any(ch.isdigit() or ch in "(-" for ch in name)


def _cationic_ring_substituent(mol, graph, root, parent, rings, atoms, order):
    """A ring-nitrogen cation as a substituent group (P-75.4, P-73.1.1.2): the neutral analogue's group name with the
    'ium' suffix cited before the free-valence suffix, as in 'pyridin-1-ium-4-yl'."""
    from ._polyfunctional import _ring_center_base

    found = _ring_center_base(mol)
    suffix = suffix_of(order)
    if found is None or found[2] != "ium" or suffix is None or found[1] not in atoms or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("a charged ring atom is not supported in this substituent group yet")
    base, center = found[0], found[1]
    named = evaluate_skeleton(base, graph, "ring", rings, atoms, [root], {parent}, suffix)
    if named is None:
        raise UnsupportedStructure("this ring cation has no supported substituent name yet")
    tail = re.search(r"-\d+(?:,\d+)*-(?:di|tri)?yl(?:idene|idyne)?$", named[1])
    if tail is None or center not in named[2]:
        raise UnsupportedStructure("the cationic ring substituent name is not delimited")
    name = f"{named[1][:tail.start()]}-{named[2][center]}-ium{tail.group(0)}"
    return name, True


def ring_carboxylate_name(mol, graph, acyl_idx, ring_atom):
    """Anion name of a ring-attached carboxylate ('benzoate', 'pyridine-3-carboxylate') and the ring atoms."""
    rings, atoms = _system_of(mol, ring_atom)
    found = evaluate_skeleton(mol, graph, "ring", rings, atoms, [ring_atom], {acyl_idx}, "carboxylate")
    if found is None:
        raise UnsupportedStructure("this ring carboxylate has no supported name yet")
    return found[1], atoms
