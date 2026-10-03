"""Stereodescriptors of metal complexes (Red Book IR-9.3, P-93.3): polyhedral
symbol, configuration index (CIP-ranked ligating atoms, priming for chelates)
and C/A or Delta/Lambda chirality, read off idealised coordinates for
T-4, SP-4, TBPY-5 and OC-6, also per metal in polynuclear complexes.
"""

import itertools
import math

from rdkit import Chem

from ._common import UnsupportedStructure

_SP_TRANS = {1: ((0, 2), (1, 3)), 2: ((0, 1), (2, 3)), 3: ((0, 3), (1, 2))}
_TB = {
    1: ((0, 4), 1), 2: ((0, 4), -1), 3: ((0, 3), 1), 4: ((0, 3), -1), 5: ((0, 2), 1), 6: ((0, 2), -1),
    7: ((0, 1), 1), 8: ((0, 1), -1), 9: ((1, 4), 1), 10: ((1, 3), 1), 11: ((1, 4), -1), 12: ((1, 3), -1),
    13: ((1, 2), 1), 14: ((1, 2), -1), 15: ((2, 4), 1), 16: ((2, 3), 1), 17: ((3, 4), 1), 18: ((3, 4), -1),
    19: ((2, 3), -1), 20: ((2, 4), -1),
}
_OH_COORDS = {
    1: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)),
    2: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
    3: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0)),
    4: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, -1.0)),
    5: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0)),
    6: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0)),
    7: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0)),
    8: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)),
    9: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0)),
    10: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0)),
    11: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0)),
    12: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0)),
    13: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0)),
    14: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, -1.0)),
    15: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0)),
    16: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0)),
    17: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0)),
    18: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    19: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0)),
    20: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0)),
    21: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0)),
    22: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0)),
    23: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0)),
    24: ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
    25: ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0)),
    26: ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0)),
    27: ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0)),
    28: ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0)),
    29: ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0)),
    30: ((0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
}


def _ring(radius_z, count, sign):
    return [
        (math.cos(sign * 2 * math.pi * k / count), math.sin(sign * 2 * math.pi * k / count), radius_z)
        for k in range(count)
    ]


def _coordinates(kind, perm, nbrs):
    """{neighbour atom: unit-ish vector}; axis atom 1 sits on +z and the ring
    is anticlockwise (right-handed about +z) for '@'."""
    if kind == "SP":
        coords = {}
        angles = (0, 180, 90, 270)
        k = 0
        for a, b in _SP_TRANS[perm]:
            for atom, ang in ((nbrs[a], angles[k]), (nbrs[b], angles[k + 1])):
                coords[atom] = (math.cos(math.radians(ang)), math.sin(math.radians(ang)), 0.0)
            k += 2
        return coords
    if kind == "TB":
        (x, y), sign = _TB[perm]
        eq = [i for i in range(5) if i not in (x, y)]
        coords = {nbrs[x]: (0.0, 0.0, 1.0), nbrs[y]: (0.0, 0.0, -1.0)}
        for atom, vec in zip((nbrs[i] for i in eq), _ring(0.0, 3, sign)):
            coords[atom] = vec
        return coords
    if kind == "OH":
        return {atom: vec for atom, vec in zip(nbrs, _OH_COORDS[perm])}
    raise UnsupportedStructure("unsupported coordination geometry")


def _tetrahedral_coordinates(nbrs, anticlockwise):
    sign = 1 if anticlockwise else -1
    coords = {nbrs[0]: (0.0, 0.0, 1.0)}
    for atom, vec in zip(nbrs[1:], _ring(-1 / 3, 3, sign)):
        coords[atom] = vec
    return coords


# ---- CIP-style ranking of ligating atoms ------------------------------------


class _Node:
    __slots__ = ("z", "kids")

    def __init__(self, z, kids):
        self.z, self.kids = z, kids

    def key(self):
        return (self.z, tuple(sorted((k.key() for k in self.kids), reverse=True)))


def _build(mol, idx, parent, path, depth):
    atom = mol.GetAtomWithIdx(idx)
    kids = []
    if depth > 0:
        for bond in atom.GetBonds():
            other = bond.GetOtherAtomIdx(idx)
            if other == parent:
                continue
            nz = mol.GetAtomWithIdx(other).GetAtomicNum()
            order = bond.GetBondTypeAsDouble()
            extra = 1 if order == 1.5 else int(order) - 1
            if bond.GetBondType() == Chem.BondType.DATIVE:
                continue
            kids.append(_Node(nz, []) if other in path else _build(mol, other, idx, path | {other}, depth - 1))
            kids.extend(_Node(nz, []) for _ in range(extra))
        kids.extend(_Node(1, []) for _ in range(atom.GetTotalNumHs()))
    return _Node(atom.GetAtomicNum(), kids)


def _signature(root, depth=8):
    levels = [[root]]
    spheres = [(root.z,)]
    for _ in range(depth):
        ordered = [sorted(n.kids, key=lambda c: c.key(), reverse=True) for n in levels[-1]]
        spheres.append(tuple(tuple(sorted((c.z for c in kids), reverse=True)) for kids in ordered))
        nxt = [c for kids in ordered for c in kids]
        if not nxt:
            break
        levels.append(nxt)
    return spheres


def priority_numbers(mol, metal_idx, donors):
    sigs = {d: _signature(_build(mol, d, metal_idx, {metal_idx, d}, 8)) for d in donors}
    uniq: list = []
    for s in sorted(sigs.values(), reverse=True):
        if not uniq or s != uniq[-1]:
            uniq.append(s)
    return {d: uniq.index(sigs[d]) + 1 for d in donors}


# ---- geometry helpers ---------------------------------------------------------


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _anticlockwise_from(coords, viewpoint, atoms):
    """`atoms` ordered anticlockwise as seen by a viewer on `viewpoint`."""
    t = coords[viewpoint]
    norm = math.sqrt(_dot(t, t))
    t = tuple(x / norm for x in t)
    ref = (1.0, 0.0, 0.0) if abs(t[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = _cross(t, _cross(ref, t))
    un = math.sqrt(_dot(u, u))
    u = tuple(x / un for x in u)
    v = _cross(t, u)

    def angle(a):
        c = coords[a]
        return math.atan2(_dot(c, v), _dot(c, u))

    return sorted(atoms, key=angle)


def _chirality_symbol(ring, number, clockwise_symbol, anticlockwise_symbol):
    """(sequence, symbol): the lower of the clockwise / anticlockwise
    priority-number readings, starting from the highest-priority atom;
    symbol '' when both readings are equal."""
    n = len(ring)
    low = min(number[a] for a in ring)
    best = None
    for first in range(n):
        if number[ring[first]] != low:
            continue
        anti = tuple(number[ring[(first + k) % n]] for k in range(n))
        clock = tuple(number[ring[(first - k) % n]] for k in range(n))
        for seq, symbol in ((clock, clockwise_symbol), (anti, anticlockwise_symbol)):
            if best is None or seq < best[0]:
                best = (seq, {symbol})
            elif seq == best[0]:
                best[1].add(symbol)
    return best[0], (next(iter(best[1])) if len(best[1]) == 1 else "")


def _rotations():
    mats = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            m = [[0] * 3 for _ in range(3)]
            for row, col in enumerate(perm):
                m[row][col] = signs[row]
            det = (
                m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
            )
            if det == 1:
                mats.append(m)
    return mats


def _octahedron_is_chiral(coords, number, alternatives):
    """True unless a rotation plus a re-priming of identical ligands maps the
    labelled octahedron onto its mirror image (so priming never creates
    chirality that the molecule does not have)."""

    def snap(vec):
        return tuple(int(round(x)) for x in vec)

    def apply(m, vec):
        return tuple(sum(m[r][c] * vec[c] for c in range(3)) for r in range(3))

    mirror = {(x, -y, z): number[a] for a, (x, y, z) in ((a, snap(coords[a])) for a in coords)}
    rotations = _rotations()
    for other in alternatives:
        positions = {snap(coords[a]): other[a] for a in coords}
        if any({apply(m, p): n for p, n in positions.items()} == mirror for m in rotations):
            return False
    return True


# ---- priming (IR-9.3.5.3) --------------------------------------------------------


def _fmt(key):
    base, primes = key
    return f"{base}{chr(0x2032) * primes}"


def _fragments(mol, metal_idx, donors, metals=()):
    """[(atoms, donors)] for every ligand fragment bound to the metal
    (all `metals` are removed so a bridging ligand stays one fragment)."""
    removed = set(metals) | {metal_idx}
    rw = Chem.RWMol(mol)
    index = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    for i in sorted(removed, reverse=True):
        rw.RemoveAtom(i)
    fragments = []
    for comp in Chem.GetMolFrags(rw.GetMol()):
        atoms = {index[i] for i in comp}
        if atoms & set(donors):
            fragments.append((atoms, [d for d in donors if d in atoms]))
    return fragments


def _fragment_keys(mol, atoms, donors):
    flagged = Chem.RWMol(mol)
    for d in donors:
        flagged.GetAtomWithIdx(d).SetIsotope(900)
    probe = flagged.GetMol()
    smiles = Chem.MolFragmentToSmiles(probe, atomsToUse=sorted(atoms), isomericSmiles=False, canonical=True)
    ranks = Chem.CanonicalRankAtomsInFragment(probe, atomsToUse=sorted(atoms), breakTies=False)
    return smiles, ranks


def _priming_assignments(mol, metal_idx, donors, metals=()):
    """Every way of priming identical polydentate ligands (and the symmetric
    halves of a ligand with three or more donors); the primed number ranks
    just below the unprimed one with the same base number."""
    fragments = [f for f in _fragments(mol, metal_idx, donors, metals) if len(f[1]) > 1]
    classes: dict[str, list] = {}
    sym_classes = []
    for atoms, ds in fragments:
        smiles, ranks = _fragment_keys(mol, atoms, ds)
        classes.setdefault(smiles, []).append((atoms, ds))
        if len(ds) >= 3:
            by_rank: dict = {}
            for d in ds:
                by_rank.setdefault(ranks[d], []).append(d)
            sym_classes.extend(group for group in by_rank.values() if len(group) > 1)
    copy_groups = [group for group in classes.values() if len(group) > 1]
    if copy_groups:
        sym_classes = []
    options = []
    for group in copy_groups:
        options.append([{d: level for level, (_, ds) in zip(perm, group) for d in ds} for perm in itertools.permutations(range(len(group)))])
    for group in sym_classes:
        options.append([{d: level for level, d in zip(perm, group)} for perm in itertools.permutations(range(len(group)))])
    if not options:
        yield {}
        return
    for combo in itertools.product(*options):
        merged: dict = {}
        for part in combo:
            for d, level in part.items():
                merged[d] = merged.get(d, 0) + level
        yield merged


# ---- descriptors ----------------------------------------------------------------


def complex_stereo_prefix(mol, metal_idx):
    """'(OC-6-22)-' style prefix, '' when the metal carries no stereo tag."""
    body = complex_stereo_descriptor(mol, metal_idx)
    return f"({body})-" if body else ""


def complex_stereo_descriptor(mol, metal_idx, metals=()):
    """'OC-6-22' style descriptor of one central atom, '' when unspecified."""
    metal = mol.GetAtomWithIdx(metal_idx)
    tag = metal.GetChiralTag()
    if tag == Chem.ChiralType.CHI_UNSPECIFIED:
        return ""
    info = next((s for s in Chem.FindPotentialStereo(mol) if s.centeredOn == metal_idx and s.specified == Chem.StereoSpecified.Specified), None)
    if info is None:
        return ""
    nbrs = list(info.controllingAtoms)
    if len(nbrs) != metal.GetDegree():
        raise UnsupportedStructure("implicit ligands on a stereo-specified metal are not supported here")
    base = priority_numbers(mol, metal_idx, nbrs)
    tetrahedral = tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    kind = {
        Chem.ChiralType.CHI_SQUAREPLANAR: "SP",
        Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL: "TB",
        Chem.ChiralType.CHI_OCTAHEDRAL: "OH",
    }.get(tag)
    if kind is None and not tetrahedral:
        raise UnsupportedStructure("this coordination geometry is not supported here")
    if tetrahedral:
        coords = _tetrahedral_coordinates(nbrs, info.descriptor == Chem.StereoDescriptor.Tet_CCW)
    else:
        coords = _coordinates(kind, int(metal.GetProp("_chiralPermutation")), nbrs)
    assignments = list(_priming_assignments(mol, metal_idx, nbrs, metals))
    numberings = [{d: (base[d], primes.get(d, 0)) for d in nbrs} for primes in assignments]
    outcomes = []
    for number in numberings:
        if tetrahedral:
            outcomes.append(("T-4", (), (), _tetrahedral(coords, nbrs, number)[1]))
        elif kind == "SP":
            index, sequence, symbol = _square_planar(coords, nbrs, number)
            outcomes.append(("SP-4", index, sequence, symbol))
        elif kind == "TB":
            index, sequence, symbol = _bipyramid(coords, nbrs, number)
            outcomes.append(("TBPY-5", index, sequence, symbol))
        else:
            index, sequence, symbol = _octahedral(coords, nbrs, number, numberings)
            outcomes.append(("OC-6", index, sequence, symbol))
    symbol_name, index_keys, sequence, _ = min(outcomes, key=lambda o: (o[1], o[2]))
    symbols = {o[3] for o in outcomes if o[1] == index_keys and o[2] == sequence}
    if len(symbols) > 1 or "?" in symbols:
        symbol = _skew_line_symbol(mol, metal_idx, nbrs, coords, metals)
    else:
        symbol = symbols.pop()
    index = "".join(_fmt(k) for k in index_keys)
    parts = [symbol_name] + ([index] if index else []) + ([symbol] if symbol else [])
    return "-".join(parts)


def _skew_line_symbol(mol, metal_idx, donors, coords, metals=()):
    """Delta/Lambda from the skew lines of two chelate rings (IR-9.3.4.11-13)."""
    bidentate = [ds for _, ds in _fragments(mol, metal_idx, donors, metals) if len(ds) == 2]
    if len(bidentate) < 2 or len(donors) != 6:
        raise UnsupportedStructure("the chirality of this chelate complex is not determined by the IUPAC priming rules")
    (a1, a2), (b1, b2) = bidentate[0], bidentate[1]
    a = tuple(x - y for x, y in zip(coords[a2], coords[a1]))
    b = tuple(x - y for x, y in zip(coords[b2], coords[b1]))
    normal = _cross(a, b)
    offset = tuple(x - y for x, y in zip(coords[b1], coords[a1]))
    helix = (1 if _dot(normal, offset) > 0 else -1) * (1 if _dot(a, b) > 0 else -1)
    return "\u0394" if helix < 0 else "\u039b"


def _trans_partner(coords):
    partner = {}
    for a, va in coords.items():
        for b, vb in coords.items():
            if a != b and _dot(va, vb) < -0.99:
                partner[a] = b
    return partner


def _square_planar(coords, nbrs, number):
    partner = _trans_partner(coords)
    lowest = min(number.values())
    firsts = [d for d in nbrs if number[d] == lowest]
    return (max(number[partner[d]] for d in firsts),), (), ""


def _bipyramid(coords, nbrs, number):
    partner = _trans_partner(coords)
    axis = [a for a in nbrs if a in partner]
    x, y = sorted(axis, key=lambda a: number[a])[:2]
    index = (number[x], number[y])
    equatorial = [a for a in nbrs if a not in axis]
    sequence, symbol = (), ""
    if number[x] != number[y] and len({number[a] for a in equatorial}) == 3:
        ring = _anticlockwise_from(coords, x, equatorial)
        sequence, symbol = _chirality_symbol(ring, number, "C", "A")
    return index, sequence, symbol


def _octahedral(coords, nbrs, number, alternatives):
    partner = _trans_partner(coords)
    lowest = min(number.values())
    best = None
    for top in [d for d in nbrs if number[d] == lowest]:
        bottom = partner[top]
        plane = [a for a in nbrs if a not in (top, bottom)]
        low = min(number[a] for a in plane)
        second = max(number[partner[a]] for a in plane if number[a] == low)
        key = (number[bottom], second)
        if best is None or key > best[0]:
            best = (key, top, bottom, plane)
    index, top, bottom, plane = best
    sequence, symbol = (), ""
    if _octahedron_is_chiral(coords, number, alternatives):
        ring = _anticlockwise_from(coords, top, plane)
        sequence, symbol = _chirality_symbol(ring, number, "C", "A")
        symbol = symbol or "?"
    return index, sequence, symbol


def _tetrahedral(coords, nbrs, number):
    if len(nbrs) != 4:
        raise UnsupportedStructure("this coordination geometry is not supported here")
    if len(set(number.values())) != 4:
        return (), ""
    ranked = sorted(nbrs, key=lambda a: number[a])
    lowest = ranked[3]
    v1, v2, v3 = (coords[a] for a in ranked[:3])
    normal = tuple(sum(c) for c in zip(_cross(v1, v2), _cross(v2, v3), _cross(v3, v1)))
    return (), "R" if _dot(normal, coords[lowest]) > 0 else "S"
