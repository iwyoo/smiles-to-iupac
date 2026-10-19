"""Mancude rings fused to a natural-product parent: fusion prefixes with primed locants (P-101.5.1)."""

from dataclasses import dataclass, field

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fused_numbering import _locant_key
from ._fullerene_fused import _bridged_alicyclic
from ._fusion_components import identify, skeleton
from ._fusion_name import fused_prefix_data
from ._np_core import loc_key
from ._numerals import multiplying_prefix

_PRIORITY = ["O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B"]
_PRIMES = ["", "′", "″", "‴"]
_DONORS = ("O", "S", "Se", "Te")


@dataclass
class Fused:
    name: str
    hetero_locants: list
    attached: list
    parent_locants: list
    indicated: list
    fusion_h: list = field(default_factory=list)
    hydro: list = field(default_factory=list)
    fusion_hydro: list = field(default_factory=list)
    cite_attached: bool = True
    prime: int = 0
    fusion_locs: tuple = ()

    @property
    def alpha(self):
        return self.name


def ring_cycle(atoms, adj):
    atoms = set(atoms)
    start = sorted(atoms)[0]
    cycle, previous = [start], None
    while len(cycle) < len(atoms):
        step = next((n for n in sorted(adj[cycle[-1]]) if n in atoms and n != previous and n not in cycle), None)
        if step is None:
            return None
        previous = cycle[-1]
        cycle.append(step)
    return cycle if cycle[0] in adj[cycle[-1]] else None


def _mono_numberings(cycle, elem):
    n = len(cycle)
    options = []
    for start in range(n):
        for step in (1, -1):
            order = [cycle[(start + step * k) % n] for k in range(n)]
            numbering = {atom: k + 1 for k, atom in enumerate(order)}
            key = tuple(sorted(numbering[a] for a in cycle if elem[a] != "C"))
            senior = tuple(numbering[a] for el in _PRIORITY for a in cycle if elem[a] == el)
            options.append((key, senior, numbering))
    best = min((o[0], o[1]) for o in options)
    return [o[2] for o in options if (o[0], o[1]) == best]


def _skeleton(atoms, view):
    index = {a: i for i, a in enumerate(atoms)}
    return skeleton(
        [view.elem[a] for a in atoms], [(index[a], index[b]) for a in atoms for b in view.adj[a] if b in index and a < b]
    )


def _ring_atoms(atoms, view):
    """The atoms of `atoms` that lie in a ring; the acyclic ones are substituents of the added component."""
    members = sorted(atoms)
    sub = _skeleton(members, view)
    info = sub.GetRingInfo()
    return {a for i, a in enumerate(members) if info.NumAtomRings(i)}


def _component(sub, view, members):
    """(locants of the heteroatoms, whether the attached-atom locants are cited, [(numbering, prefix, ranking key)]) of the
    ring or ring system `sub` as a fusion prefix (P-25.3.2); a bridged alicyclic system takes its von Baeyer prefix."""
    found = identify(sub)
    if found:
        component = found[0]
        cite = not (component.kind in ("mono", "attached_only") and all(e == "C" for e in component.elements))
        hetero = component.hetero_text.split(",") if component.hetero_text else []
        return hetero, cite, [(n, component.prefix, ()) for n in component.numberings]
    try:
        bracketed, numberings = fused_prefix_data(sub)
    except UnsupportedStructure:
        kekule = Chem.Mol(view.mol)
        Chem.Kekulize(kekule, clearAromaticFlags=True)
        bridged = _bridged_alicyclic(view.mol, kekule, sub, members, {a: i for i, a in enumerate(members)})
        if bridged is None:
            raise UnsupportedStructure("an added ring system that has no fusion prefix") from None
        return [], True, [(n, text, key) for n, key, text in zip(*bridged)]
    return [], True, [(n, bracketed, ()) for n in numberings]


def _locant_text(number, prime):
    return f"{number}{_PRIMES[prime]}"


def name_fused(comp, cand, view, final, parent_centers, hfaces):
    """Describe one ring fused to the skeleton at a bond (or two) of the parent."""
    image = {a: loc for loc, a in cand.mapping.items()}
    attachments = sorted({a for a, _ in comp.links})
    if len(attachments) != 2:
        raise UnsupportedStructure("an added ring attached at other than two skeleton atoms")
    sa, sb = attachments
    skeleton_adj = {
        atom: {n for n in view.adj[atom] if n in image and image[n] in cand.skel.adj[image[atom]]} for atom in image
    }
    path = _shortest(skeleton_adj, sa, sb)
    if path is None or len(path) > 3:
        raise UnsupportedStructure("an added ring that spans more than three skeleton atoms")
    ring = _ring_atoms(set(comp.atoms) | set(path), view)
    members = sorted(ring)
    index = {a: i for i, a in enumerate(members)}
    sub = _skeleton(members, view)
    hetero_text, cite_attached, options = _component(sub, view, members)
    derived = any(key for _, _, key in options)
    candidates = [({a: numbering[index[a]] for a in members}, name, key) for numbering, name, key in options]

    def fusion_key(numbering):
        return tuple(sorted(_locant_key(numbering[a]) for a in path))

    ordered_path = sorted(path, key=lambda a: loc_key(final(image[a])))
    chosen, name, _ = min(
        candidates, key=lambda c: (c[2], fusion_key(c[0]), tuple(_locant_key(c[0][a]) for a in ordered_path))
    )
    indicated, hydro, fusion_roles = ([], [], {}) if derived else _saturated(comp, path, ring, view, chosen, image, final)
    return Fused(
        name=name,
        hetero_locants=hetero_text,
        attached=[chosen[a] for a in ordered_path],
        parent_locants=[final(image[a]) for a in ordered_path],
        indicated=indicated,
        hydro=hydro,
        fusion_hydro=[final(image[a]) for a, role in fusion_roles.items() if role == "hydro"],
        fusion_h=[
            (final(image[a]), a, fusion_roles[a])
            for a in ordered_path
            if a in fusion_roles and image[a] not in parent_centers
        ],
        cite_attached=cite_attached,
        fusion_locs=tuple(image[a] for a in path),
    )


def _shortest(adj, start, goal):
    frontier, parents = [start], {start: None}
    while frontier:
        nxt = []
        for a in frontier:
            for n in adj[a]:
                if n not in parents:
                    parents[n] = a
                    nxt.append(n)
        frontier = nxt
    if goal not in parents:
        return None
    path = [goal]
    while parents[path[-1]] is not None:
        path.append(parents[path[-1]])
    return path[::-1]


def _double_in(atom, ring, view):
    return any(view.order.get(frozenset((atom, n)), 1) == 2 for n in view.adj[atom] if n in ring)


def _saturated(comp, path, ring, view, numbering, image, final):
    """(indicated locants of the component, hydro locants of the component, {fusion atom: role}) of a fused ring.

    Saturated positions of the component, fusion atoms included, that keep a hydrogen atom are cited with indicated
    hydrogen when odd in number and otherwise with 'hydro' prefixes; fusion atoms keep their parent locants
    (P-101.5.1.2, P-101.6.5)."""
    atoms = sorted(ring)
    positions = [
        a for a in atoms
        if view.elem[a] not in _DONORS
        and not _double_in(a, ring, view)
        and view.mol.GetAtomWithIdx(a).GetTotalNumHs() > 0
    ]

    def value(atom):
        text = str(numbering[atom]) if atom in comp.atoms else str(final(image[atom]))
        digits = "".join(c for c in text if c.isdigit())
        return (int(digits or 0), loc_key(text))

    positions.sort(key=value)
    indicated = []
    if len(positions) % 2:
        indicated = [positions.pop(0)]
    roles = {a: "hydro" for a in positions if a in path and a not in comp.atoms}
    roles.update({a: "indicated" for a in indicated if a in path and a not in comp.atoms})
    hydro = [numbering[a] for a in positions if a in comp.atoms]
    return [numbering[a] for a in indicated if a in comp.atoms], hydro, roles


def _fusion_h(atom, view, ring, image, parent_centers):
    if _double_in(atom, ring, view) or image[atom] in parent_centers:
        return False
    return view.mol.GetAtomWithIdx(atom).GetTotalNumHs() > 0


def assign_primes(fused):
    """Order the fused components for citation and give each its prime (P-101.5.1.1)."""
    ordered = []
    groups = {}
    for f in fused:
        groups.setdefault((f.name, tuple(f.hetero_locants)), []).append(f)
    for key in sorted(groups, key=lambda k: k[0]):
        ordered.append(groups[key])
    prime = 0
    for items in ordered:
        for f in items:
            prime += 1
            f.prime = prime
    return ordered


def fusion_text(groups):
    text = ""
    for items in groups:
        brackets = []
        for f in items:
            attached = ",".join(_locant_text(n, f.prime) for n in f.attached) if f.cite_attached else ""
            parent = ",".join(f.parent_locants)
            brackets.append(f"{attached}:{parent}" if attached else parent)
        hetero_text = f"[{','.join(str(h) for h in items[0].hetero_locants)}]" if items[0].hetero_locants else ""
        multiplier = "" if len(items) == 1 else multiplying_prefix(len(items), compound=bool(hetero_text))
        text += f"{multiplier}{hetero_text}{items[0].name}[{';'.join(brackets)}]"
    return text


def indicated_texts(groups):
    found = []
    for items in groups:
        for f in items:
            found += [(f.prime, n) for n in f.indicated]
    return [f"{n}{_PRIMES[prime]}H" for prime, n in sorted(found)]


def hydro_texts(groups):
    """Hydro locants of the fused components as (sort value, text); component locants carry their prime."""
    found = []
    for items in groups:
        for f in items:
            found += [f"{n}{_PRIMES[f.prime]}" for n in f.hydro]
            found += list(f.fusion_hydro)
    return found
