"""Connected regions of a molecule that a modified parent skeleton can occupy (P-101.3.7.1, P-101.3.7.2)."""

from itertools import combinations
from math import comb
from typing import NamedTuple

from rdkit import Chem

_MAX_BRIDGE_SETS = 2000


class Region(NamedTuple):
    size: int
    cyclomatic: int
    stereo: bool


def _stereo_elements(view):
    cached = view.__dict__.get("stereo_elements")
    if cached is None:
        mol = view.mol
        atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE]
        cached = view.stereo_elements = (atoms, bonds)
    return cached


def regions(view, elements, spare):
    """Connected atom sets built from atoms of `elements` and at most `spare` bridging atoms of other elements.

    A skeleton is connected and each of its atoms is either of the parent's elements or a replacement, which counts as
    a modification; an atom with one neighbour cannot be replaced (P-101.3.7.1)."""
    cache = view.__dict__.setdefault("region_cache", {})
    key = (frozenset(elements), spare)
    if key not in cache:
        cache[key] = _regions(view, elements, spare)
    return cache[key]


def _regions(view, elements, spare):
    base = {a for a, e in view.elem.items() if e in elements}
    bridges = [a for a, e in view.elem.items() if e not in elements and len(view.adj[a]) >= 2]
    if sum(comb(len(bridges), k) for k in range(spare + 1)) > _MAX_BRIDGE_SETS:
        base, spare = base | set(bridges), 0
    stereo_atoms, stereo_bonds = _stereo_elements(view)
    seen, found = set(), []
    for k in range(spare + 1):
        for chosen in combinations(bridges, k):
            allowed = base | set(chosen)
            done = set()
            for start in allowed:
                if start in done:
                    continue
                region, stack = {start}, [start]
                while stack:
                    for n in view.adj[stack.pop()]:
                        if n in allowed and n not in region:
                            region.add(n)
                            stack.append(n)
                done |= region
                frozen = frozenset(region)
                if frozen in seen:
                    continue
                seen.add(frozen)
                bonds = sum(1 for a in region for n in view.adj[a] if n in region) // 2
                stereo = bool(region & stereo_atoms) or any(a in region and b in region for a, b in stereo_bonds)
                found.append(Region(len(region), bonds - len(region) + 1, stereo))
    return found
