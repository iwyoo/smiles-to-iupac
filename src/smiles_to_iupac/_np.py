"""Names of natural products on the semisystematic parents of Table 10.1 (P-101)."""

from collections import Counter
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure
from ._pin import mark
from ._np_build import _homo_position, build
from ._np_core import PARENTS, get_parent, loc_key
from ._np_diff import WORK_LIMIT, lower_bound, read_operations
from ._np_match import View, embeddings, skeleton_has_stereo
from ._np_name import _SENIORITY, classify
from ._np_rings import components
from ._np_skel import Skel, variants
from ._np_text import stem_info

PREFERRED_OPERATIONS = 2
_NO_PIN = "P-101 identifies no preferred IUPAC names for natural-product parent structures"
SYSTEMATIC_PREFERRED = {"tropane", "bornane", "carane", "fenchane", "pinane", "thujane", "p-menthane", "bisabolane"}
_MIN_HEAVY_ATOMS = 9
_MAX_HEAVY_ATOMS = 90
_MAX_RINGS = 12
_EXTRA_AROMATIC_ATOMS = 12
_MIN_RINGS_FOR_OPERATIONS = 3
_MAX_ATOMS_REMOVED_PER_OPERATION = 12


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
def _variants(name, cost, terminal_only=False, special_only=False):
    return tuple(variants(get_parent(name), cost, terminal_only, special_only))


_MAX_SPECIAL_COST = 2


def _may_modify(name):
    parent = get_parent(name)
    if name.endswith(("carotene", "neolignane")):
        return True
    return len(parent.adj) >= _MIN_PARENT_ATOMS_FOR_MODIFICATION and _parent_facts(name)[4] >= _MIN_RINGS_FOR_OPERATIONS


def _skeletons(name, view, cost, terminal_only):
    parent = get_parent(name)
    if cost == 0:
        return [Skel.of(parent)]
    if not _may_modify(name):
        return []
    found = read_operations(parent, view, cost, terminal_only)
    plain = list(_variants(name, cost, terminal_only)) if found is None and cost <= _MAX_SPECIAL_COST else list(found or [])
    special = _variants(name, cost, terminal_only, True) if cost <= _MAX_SPECIAL_COST else ()
    return plain + list(special)


def _plausible(name, view, cost):
    """Cheap necessary conditions for a parent (modified by at most `cost` operations) to occur in the molecule."""
    cyc, rings, elements, aromatic, fused = _parent_facts(name)
    if cost and -(-(fused - 1) // 2) > view.fused:
        return False
    if cyc - cost > view.cyclomatic or len(get_parent(name).adj) - cost * _MAX_ATOMS_REMOVED_PER_OPERATION > len(view.adj):
        return False
    if view.aromatic_atoms > aromatic + _EXTRA_AROMATIC_ATOMS:
        return False
    missing = sum(max(0, n - view.rings_by_size.get(size, 0)) for size, n in rings.items())
    if missing > 2 * cost + 2:
        return False
    return _replacements_needed(name, view) + cost <= _MAX_MODIFICATIONS


def _replacements_needed(name, view):
    """Heteroatoms of the parent that the molecule lacks: each needs a skeletal replacement, which counts as a modification."""
    elements = _parent_facts(name)[2]
    return sum(max(0, n - view.elements.get(el, 0)) for el, n in elements.items() if el != "C")


def _rank(cand, view):
    comps = components(view, set(cand.mapping.values()))
    ignored = set().union(*(c.atoms for c in comps)) if comps else set()
    groups = classify(cand, view, ignored)
    principal = next((c for c in _SENIORITY if c in groups.classes), None)
    locs = tuple(sorted(loc_key(m[0]) for m in groups.classes.get(principal, []))) if principal else ()
    everything = tuple(sorted(loc_key(m[0]) for members in groups.classes.values() for m in members))
    branches = tuple(sorted(loc_key(loc) for loc, _ in groups.branches))
    return locs, everything, branches


_MAX_MODIFICATIONS = 4
_MIN_PARENT_ATOMS_FOR_MODIFICATION = 14


def _admissible(cand, view):
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
    open_chain = cand.parent.name.endswith(("carotene", "neolignane"))
    if cand.parent.name.endswith("neolignane") and (cand.skel.ops or cand.replaced) and view.cyclomatic <= cyc:
        return False
    if not open_chain and cand.skel.ops and fused < _MIN_RINGS_FOR_OPERATIONS:
        return False
    if any(op[0] == "des" for op in cand.skel.ops) and (len(cand.skel.ops) > 1 or cand.cyclo or cand.replaced):
        return False
    return True


def _branch_homo_outside_ring(cand, view):
    return any(
        op[0] == "homo" and op[1] == "terminal" and len(cand.parent.adj[op[2]]) > 1 and cand.mapping[op[3]] not in view.rings
        for op in cand.skel.ops
    )


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
        if not _admissible(cand, view):
            continue
        if cand.skel.ops and not cand.cyclo and not _rings_contained(cand, view):
            continue
        key = (cand.parent.name, tuple(sorted(map(str, (op[:2] for op in cand.skel.ops)))), frozenset(cand.mapping.values()))
        try:
            rank = _rank(cand, view)
        except UnsupportedStructure:
            continue
        order = sorted(cand.mapping, key=loc_key)
        inserted = tuple(sorted((-_homo_position(op) for op in cand.skel.ops if op[0] == "homo")))
        entry = (rank, inserted, tuple(cand.mapping[l] for l in order))
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


def _skeleton_facts(skel):
    if not hasattr(skel, "_facts"):
        skel._facts = (len(skel.adj), _cyclomatic(skel.adj), _skeleton_rings(skel), Counter(skel.elem.values()))
    return skel._facts


_RING_SLACK = 3


def _candidates(skel, view):
    size, cyclomatic, rings, elements = _skeleton_facts(skel)
    if size > len(view.adj) or cyclomatic > view.cyclomatic:
        return []
    secos = sum(1 for op in skel.ops if op[0] == "seco")
    if secos > 1 or (secos and cyclomatic < 2):
        return []
    budget = _MAX_MODIFICATIONS - len(skel.ops)
    missing = sum(max(0, n - view.rings_by_size.get(length, 0)) for length, n in rings.items())
    if missing > 2 * budget + _RING_SLACK:
        return []
    if sum(max(0, n - view.elements.get(e, 0)) for e, n in elements.items() if e != "C") > budget:
        return []
    return embeddings(skel, view, limit=3000)


def _expected_cost(name, view):
    """Operations the parent needs at least: intact-group bound, replacements, and rings the molecule has beyond the parent."""
    extra_rings = max(0, view.cyclomatic - _parent_facts(name)[0])
    return lower_bound(get_parent(name), view, _MAX_MODIFICATIONS) + _replacements_needed(name, view) + extra_rings


def _promising(view, cost, exact):
    """Plausible parents, most promising first: fewest operations to leave a group intact, then closest in size."""
    names = [n for n in PARENTS if n not in SYSTEMATIC_PREFERRED and n not in exact and _plausible(n, view, cost)]
    if cost:
        names = [n for n in names if _may_modify(n)]
        names.sort(key=lambda n: (_expected_cost(n, view), abs(len(get_parent(n).adj) - len(view.adj)), n))
    return names


def _name_once(mol):
    heavy = mol.GetNumAtoms()
    if heavy < _MIN_HEAVY_ATOMS or heavy > _MAX_HEAVY_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    view = View(mol)
    view.cyclomatic = _cyclomatic(view.adj)
    view.fused = _fused_rings(mol)
    if view.cyclomatic > _MAX_RINGS:
        return None
    view.rings_by_size = Counter(len(r) for r in mol.GetRingInfo().AtomRings())
    view.elements = Counter(view.elem.values())
    view.has_stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )
    view.aromatic_atoms = len({i for bond in view.aromatic_bonds for i in bond})
    view.work = WORK_LIMIT
    best = None
    exact = set()
    for cost in range(_MAX_MODIFICATIONS + 1) if view.has_stereo else (0, 1):
        if best is not None and (best.cost < cost or (best.cost == cost and not best.key[1])):
            break
        found = []
        for name in _promising(view, cost, exact):
            terminal_only = view.cyclomatic < 2
            skels = _skeletons(name, view, cost, terminal_only)
            for skel in skels:
                found += _candidates(skel, view)
        embedded = {c.parent.name for c in found} if cost == 0 else set()
        blocked, built_names = set(), set()
        for cand in _best_per_skeleton(found, view):
            configured = skeleton_has_stereo(cand, view)
            if cand.replaced and set(view.elements) == {"C"} and not configured:
                continue
            if len(cand.skel.ops) + len(cand.cyclo) + len(cand.replaced) > 1 and not configured:
                continue
            if not configured and any(op[0] == "des" for op in cand.skel.ops):
                continue
            if _branch_homo_outside_ring(cand, view):
                continue
            try:
                built = build(cand, view)
            except UnsupportedStructure as error:
                if "acyclic part" in str(error):
                    blocked.add(cand.parent.name)
                continue
            built_names.add(cand.parent.name)
            if best is None or built.key < best.key or (built.key == best.key and built.name < best.name):
                best = built
        if cost == 0:
            exact = embedded - (blocked - built_names)
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
    return name_natural_product_ranked(mol)[0]


def name_natural_product_ranked(mol):
    """(name, number of modifying operations incl. replacement): more than two have no preferred semisystematic name (P-101.3.7.2)."""
    name, skeletal = _natural_product_name(mol)
    return (mark(name, _NO_PIN), skeletal) if name else (None, 0)


def _natural_product_name(mol):
    """(name, operations): natural-product name, with 'ent' for a full inversion and 'rac'/'rel' for stereo groups (P-101.8)."""
    groups = mol.GetStereoGroups()
    if groups:
        return _name_with_groups(mol, groups)
    plain = _plain(mol)
    if plain is None:
        return None, 0
    built = _name_once(plain)
    if built is None:
        return None, 0
    if built.implied_total and built.implied_cited == built.implied_total:
        other = _name_once(_mirror(plain))
        if other is not None and other.implied_cited == 0 and other.implied_total == built.implied_total:
            return f"ent-{other.name}", other.skeletal
    return built.name, built.skeletal


def _name_with_groups(mol, groups):
    kinds = {g.GetGroupType() for g in groups}
    if len(groups) != 1 or len(kinds) != 1:
        return None, 0
    kind = next(iter(kinds))
    chiral = {a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
    if {a.GetIdx() for a in groups[0].GetAtoms()} != chiral:
        return None, 0
    word = {Chem.StereoGroupType.STEREO_AND: "rac", Chem.StereoGroupType.STEREO_OR: "rel"}.get(kind)
    if word is None:
        return None, 0
    plain = _plain(mol)
    built = _name_once(plain)
    if built is None:
        return None, 0
    if built.first_face == "β":
        other = _name_once(_mirror(plain))
        if other is not None and other.first_face == "α":
            built = other
    return f"{word}-{built.name}", built.skeletal
