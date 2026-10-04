"""Names of natural products on the semisystematic parents of Table 10.1 (P-101)."""

from collections import Counter
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure
from ._np_build import build
from ._np_core import get_parent, loc_key
from ._np_match import View, embeddings
from ._np_name import _SENIORITY, classify
from ._np_parents import PARENTS
from ._np_rings import components
from ._np_skel import Skel, variants
from ._np_text import stem_info

SYSTEMATIC_PREFERRED = {"tropane", "bornane", "carane", "fenchane", "pinane", "thujane", "p-menthane", "bisabolane"}
_MIN_HEAVY_ATOMS = 9
_MAX_HEAVY_ATOMS = 90
_MAX_RINGS = 12
_EXTRA_AROMATIC_ATOMS = 12
_MIN_RINGS_FOR_OPERATIONS = 3


def _cyclomatic(adj):
    edges = sum(len(n) for n in adj.values()) // 2
    seen, components = set(), 0
    for start in adj:
        if start in seen:
            continue
        components += 1
        stack = [start]
        seen.add(start)
        while stack:
            for n in adj[stack.pop()]:
                if n not in seen:
                    seen.add(n)
                    stack.append(n)
    return edges - len(adj) + components


def _fused_rings(mol):
    """Number of rings in the largest group of rings that share bonds."""
    rings = [set(r) for r in mol.GetRingInfo().BondRings()]
    groups = []
    for ring in rings:
        touching = [g for g in groups if any(ring & other for other in g[0])]
        merged = [ring]
        count = 1
        for g in touching:
            merged += g[0]
            count += g[1]
            groups.remove(g)
        groups.append((merged, count))
    return max((g[1] for g in groups), default=0)


@lru_cache(maxsize=None)
def _parent_facts(name):
    parent = get_parent(name)
    rings = Counter(len(r) for r in parent.mol.GetRingInfo().AtomRings())
    return _cyclomatic(parent.adj), rings, Counter(parent.elem.values()), len(parent.aromatic_atoms), _fused_rings(parent.mol)


@lru_cache(maxsize=None)
def _variants(name, cost, terminal_only=False):
    return tuple(variants(get_parent(name), cost, terminal_only))


def _plausible(name, view, cost):
    """Cheap necessary conditions for a parent (modified by at most `cost` operations) to occur in the molecule."""
    cyc, rings, elements, aromatic, _ = _parent_facts(name)
    if cyc - cost > view.cyclomatic or len(PARENTS) and get_parent(name).adj.__len__() - cost > len(view.adj):
        return False
    if view.aromatic_atoms > aromatic + _EXTRA_AROMATIC_ATOMS:
        return False
    missing = sum(max(0, n - view.rings_by_size.get(size, 0)) for size, n in rings.items())
    if missing > 2 * cost + 2:
        return False
    hetero_missing = sum(max(0, n - view.elements.get(el, 0)) for el, n in elements.items())
    return hetero_missing <= cost + _MAX_MODIFICATIONS


def _rank(cand, view):
    comps = components(view, set(cand.mapping.values()))
    ignored = set().union(*(c.atoms for c in comps)) if comps else set()
    groups = classify(cand, view, ignored)
    principal = next((c for c in _SENIORITY if c in groups.classes), None)
    locs = tuple(sorted(loc_key(m[0]) for m in groups.classes.get(principal, []))) if principal else ()
    everything = tuple(sorted(loc_key(m[0]) for members in groups.classes.values() for m in members))
    branches = tuple(sorted(loc_key(loc) for loc, _ in groups.branches))
    return locs, everything, branches


_MAX_MODIFICATIONS = 2
_MIN_PARENT_ATOMS_FOR_MODIFICATION = 14


def _admissible(cand):
    modifications = len(cand.skel.ops) + len(cand.cyclo) + len(cand.replaced)
    if modifications > _MAX_MODIFICATIONS:
        return False
    if modifications and len(cand.parent.adj) < _MIN_PARENT_ATOMS_FOR_MODIFICATION:
        return False
    facts = _parent_facts(cand.parent.name)
    cyc, fused = facts[0], facts[4]
    if cand.cyclo and fused < _MIN_RINGS_FOR_OPERATIONS:
        return False
    if cand.replaced and not (fused >= 3 or (cyc >= 2 and stem_info(cand.parent.name)[1] == "ane")):
        return False
    carotene = cand.parent.name.endswith("carotene")
    if not carotene and cand.skel.ops and fused < _MIN_RINGS_FOR_OPERATIONS:
        return False
    if any(op[0] == "des" for op in cand.skel.ops) and (len(cand.skel.ops) > 1 or cand.cyclo or cand.replaced):
        return False
    return True


def _terminal_op(parent, op):
    if op[0] == "nor":
        return len(parent.adj[op[1]]) == 1
    return op[1] == "terminal"


def _skeleton_rings(skel):
    mol = Chem.RWMol()
    index = {a: mol.AddAtom(Chem.Atom(6)) for a in skel.adj}
    for a, neighbors in skel.adj.items():
        for b in neighbors:
            if index[a] < index[b]:
                mol.AddBond(index[a], index[b], Chem.BondType.SINGLE)
    return Counter(len(r) for r in Chem.GetSymmSSSR(mol))


def _rings_contained(cand, view):
    """The rings of the modified parent are rings of the molecule too (a bridge or fusion may add more)."""
    extra = _skeleton_rings(cand.skel) - view.rings_by_size
    return not extra


def _best_per_skeleton(cands, view):
    """One candidate per distinct set of mapped atoms: the numbering with the lowest locants for the groups."""
    best = {}
    for cand in cands:
        if not _admissible(cand):
            continue
        if cand.skel.ops and not _rings_contained(cand, view):
            continue
        key = (cand.parent.name, tuple(sorted(map(str, (op[:2] for op in cand.skel.ops)))), frozenset(cand.mapping.values()))
        try:
            rank = _rank(cand, view)
        except UnsupportedStructure:
            continue
        order = sorted(cand.mapping, key=loc_key)
        entry = (rank, tuple(cand.mapping[l] for l in order))
        if key not in best or entry < best[key][0]:
            best[key] = (entry, cand)
    return [c for _, c in best.values()]


def _terminal_only(skel):
    parent = skel.parent
    for op in skel.ops:
        if op[0] == "nor" and len(parent.adj[op[1]]) != 1:
            return False
        if op[0] == "homo" and op[1] != "terminal":
            return False
        if op[0] in ("seco", "cyclo"):
            return False
    return True


def _candidates(skel, view):
    if len(skel.adj) > len(view.adj) or _cyclomatic(skel.adj) > view.cyclomatic:
        return []
    secos = sum(1 for op in skel.ops if op[0] == "seco")
    if secos > 1 or (secos and _cyclomatic(skel.adj) < 2):
        return []
    return embeddings(skel, view, limit=3000)


def _name_once(mol):
    heavy = mol.GetNumAtoms()
    if heavy < _MIN_HEAVY_ATOMS or heavy > _MAX_HEAVY_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    view = View(mol)
    view.cyclomatic = _cyclomatic(view.adj)
    if view.cyclomatic > _MAX_RINGS:
        return None
    view.rings_by_size = Counter(len(r) for r in mol.GetRingInfo().AtomRings())
    view.elements = Counter(view.elem.values())
    view.has_stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )
    view.aromatic_atoms = len({i for bond in view.aromatic_bonds for i in bond})
    best = None
    for cost in (0, 1, 2):
        if best is not None and (best.cost < cost or (best.cost == cost and not best.key[1])):
            break
        found = []
        for name in PARENTS:
            if name in SYSTEMATIC_PREFERRED or not _plausible(name, view, cost):
                continue
            parent = get_parent(name)
            terminal_only = view.cyclomatic < 2
            skels = [Skel.of(parent)] if cost == 0 else _variants(name, cost, terminal_only)
            for skel in skels:
                found += _candidates(skel, view)
        for cand in _best_per_skeleton(found, view):
            if cand.replaced and set(view.elements) == {"C"}:
                continue
            if len(cand.skel.ops) + len(cand.cyclo) + len(cand.replaced) > 1 and not view.has_stereo:
                continue
            try:
                built = build(cand, view)
            except UnsupportedStructure:
                continue
            if best is None or built.key < best.key or (built.key == best.key and built.name < best.name):
                best = built
    return best


def _mirror(mol):
    mirror = Chem.RWMol(mol)
    for atom in mirror.GetAtoms():
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return Chem.MolFromSmiles(Chem.MolToSmiles(mirror))


def _plain(mol):
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))


def name_natural_product(mol):
    """Natural-product name, with 'ent' for a full inversion and 'rac'/'rel' for stereo groups (P-101.8)."""
    groups = mol.GetStereoGroups()
    if groups:
        return _name_with_groups(mol, groups)
    plain = _plain(mol)
    built = _name_once(plain)
    if built is None:
        return None
    if built.implied_total and built.implied_cited == built.implied_total:
        other = _name_once(_mirror(plain))
        if other is not None and other.implied_cited == 0 and other.implied_total == built.implied_total:
            return f"ent-{other.name}"
    return built.name


def _name_with_groups(mol, groups):
    kinds = {g.GetGroupType() for g in groups}
    if len(groups) != 1 or len(kinds) != 1:
        return None
    kind = next(iter(kinds))
    chiral = {a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
    if {a.GetIdx() for a in groups[0].GetAtoms()} != chiral:
        return None
    word = {Chem.StereoGroupType.STEREO_AND: "rac", Chem.StereoGroupType.STEREO_OR: "rel"}.get(kind)
    if word is None:
        return None
    plain = _plain(mol)
    built = _name_once(plain)
    if built is None:
        return None
    if built.first_face == "β":
        other = _name_once(_mirror(plain))
        if other is not None and other.first_face == "α":
            built = other
    return f"{word}-{built.name}"
