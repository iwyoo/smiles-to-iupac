"""Stereoparent structures of Table 10.1 (P-101.2.7): numbered graphs, layouts and α/β faces (P-101.2.6)."""

import math
import re
from functools import lru_cache

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor

from ._np_parents import PARENTS

_PRIMES = {"′": 1, "″": 2}
_SUPERSCRIPTS = {"¹": 1, "²": 2, "³": 3}
FACES = {"a": "α", "b": "β", "x": "ξ"}
_STEROIDS = (
    "gonane estrane androstane pregnane cholane cholestane ergostane campestane stigmastane poriferastane gorgostane "
    "cardanolide bufanolide furostan spirostan spirosolane solanidane conanine"
).split()
PLANE_EXCLUDED = {"morphinan": {"15", "16", "17"}}
UNIMPLIED = {name: {"5"} for name in _STEROIDS}
UNIMPLIED.update({"spirostan": {"5", "22", "25"}, "spirosolane": {"5", "22", "25"}, "furostan": {"5", "22"}})


def loc_key(loc):
    primes = 0
    while loc and loc[-1] in _PRIMES:
        primes += _PRIMES[loc[-1]]
        loc = loc[:-1]
    match = re.match(r"(\d+)(.*)", loc)
    base, rest = (int(match.group(1)), match.group(2)) if match else (10**6, loc)
    sup = sum(_SUPERSCRIPTS.get(c, 0) for c in rest)
    letters = "".join(c for c in rest if c not in _SUPERSCRIPTS)
    return (primes, base, sup, letters)


def is_numbered(loc):
    return loc_key(loc)[1] < 10**6


class Parent:
    def __init__(self, name):
        smiles, locants, anchor = PARENTS[name]
        self.name = name
        mol = Chem.MolFromSmiles(smiles)
        for atom, loc in zip(mol.GetAtoms(), locants.split()):
            atom.SetProp("loc", loc)
        self.loc_of = {a.GetIdx(): a.GetProp("loc") for a in mol.GetAtoms()}
        self.aromatic_bonds = {
            frozenset((self.loc_of[b.GetBeginAtomIdx()], self.loc_of[b.GetEndAtomIdx()]))
            for b in mol.GetBonds()
            if b.GetIsAromatic()
        }
        self.aromatic_atoms = {loc for bond in self.aromatic_bonds for loc in bond}
        centers = {e.centeredOn: e.specified for e in Chem.FindPotentialStereo(mol) if e.type == Chem.StereoType.Atom_Tetrahedral}
        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        self.mol = mol
        self.centers = {self.loc_of[i] for i in centers}
        self.implied = {self.loc_of[i] for i, state in centers.items() if state == Chem.StereoSpecified.Specified} - UNIMPLIED.get(name, set())
        self.idx_of = {loc: idx for idx, loc in self.loc_of.items()}
        self.order = sorted(self.idx_of, key=loc_key)
        self.elem = {loc: mol.GetAtomWithIdx(i).GetSymbol() for loc, i in self.idx_of.items()}
        self.adj = {loc: set() for loc in self.idx_of}
        self.bond_order = {}
        for bond in mol.GetBonds():
            a, b = self.loc_of[bond.GetBeginAtomIdx()], self.loc_of[bond.GetEndAtomIdx()]
            self.adj[a].add(b)
            self.adj[b].add(a)
            self.bond_order[frozenset((a, b))] = int(bond.GetBondTypeAsDouble())
        ring_info = mol.GetRingInfo()
        self.ring_atoms = {self.loc_of[i] for ring in ring_info.AtomRings() for i in ring}
        self.ring_bonds = {
            frozenset((self.loc_of[b.GetBeginAtomIdx()], self.loc_of[b.GetEndAtomIdx()]))
            for b in mol.GetBonds()
            if b.IsInRing()
        }
        self.anchor = tuple(anchor.split(":")) if anchor else None
        self.plane = self.ring_atoms - PLANE_EXCLUDED.get(name, set())
        self.plane_bonds = {b for b in self.ring_bonds if b <= self.plane}
        self._layout = None

    def ending_kind(self):
        return "ane" if re.search(r"(ane|an|anine|stane|ostan)$", self.name) else "other"

    def layout(self):
        if self._layout is None:
            self._layout = _compute_layout(self)
        return self._layout


def _segments_cross(p1, p2, p3, p4):
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    return orient(p1, p2, p3) * orient(p1, p2, p4) < -1e-9 and orient(p3, p4, p1) * orient(p3, p4, p2) < -1e-9


def _compute_layout(parent):
    """(coordinates of the ring atoms, orientation sign) or None when no faithful planar drawing exists."""
    ring_atoms = parent.plane
    if not ring_atoms or parent.anchor is None:
        return None
    mol = Chem.RWMol(parent.mol)
    for loc in sorted((l for l in parent.idx_of if l in PLANE_EXCLUDED.get(parent.name, set())), key=lambda l: -parent.idx_of[l]):
        mol.RemoveAtom(parent.idx_of[loc])
    kept = [l for l in sorted(parent.idx_of, key=lambda l: parent.idx_of[l]) if l not in PLANE_EXCLUDED.get(parent.name, set())]
    mol = mol.GetMol()
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xy = {loc: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i, loc in enumerate(kept)}
    ring = sorted(ring_atoms)
    for a, b in ((a, b) for k, a in enumerate(ring) for b in ring[k + 1 :]):
        if math.dist(xy[a], xy[b]) < 0.45:
            return None
    bonds = [tuple(sorted(bd)) for bd in parent.plane_bonds]
    for k, (a, b) in enumerate(bonds):
        for c, d in bonds[k + 1 :]:
            if len({a, b, c, d}) == 4 and _segments_cross(xy[a], xy[b], xy[c], xy[d]):
                return None
    return {loc: xy[loc] for loc in ring_atoms}


_TAG_OF_VOLUME = {}


def _volume_sign(center, neighbors):
    c = np.array(center)
    n = [np.array(p) - c for p in neighbors]
    return 1 if float(np.dot(np.cross(n[1] - n[0], n[2] - n[0]), n[3] - n[0])) > 0 else -1


def _calibrate():
    for smiles in ("N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O", "F[C@](Cl)(Br)I"):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        from rdkit.Chem import AllChem

        AllChem.EmbedMolecule(mol, randomSeed=7)
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            if atom.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW) and atom.GetDegree() == 4:
                pos = [tuple(conf.GetAtomPosition(n.GetIdx())) for n in atom.GetNeighbors()]
                sign = _volume_sign(tuple(conf.GetAtomPosition(atom.GetIdx())), pos)
                _TAG_OF_VOLUME.setdefault(sign, atom.GetChiralTag())
                assert _TAG_OF_VOLUME[sign] == atom.GetChiralTag()


_calibrate()
_CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW


def _tag_sign(tag):
    return 1 if tag == _CW else -1


def center_signature(mol_h, atom_idx, labels):
    """±1: handedness of the atom's neighbours read in the order of their sorted `labels` (a label per neighbour)."""
    atom = mol_h.GetAtomWithIdx(atom_idx)
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    if len(neighbors) != 4 or atom.GetChiralTag() not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
        return None
    sequence = [labels[n] for n in neighbors]
    ranked = sorted(range(4), key=lambda k: sequence[k])
    inversions = sum(1 for i in range(4) for j in range(i + 1, 4) if ranked[i] > ranked[j])
    sign = _tag_sign(atom.GetChiralTag())
    return sign if inversions % 2 == 0 else -sign


def _outward(origin, neighbor_xy):
    vector = np.zeros(2)
    for p in neighbor_xy:
        d = np.array(p) - np.array(origin)
        norm = np.linalg.norm(d)
        if norm > 1e-9:
            vector -= d / norm
    norm = np.linalg.norm(vector)
    return vector / norm if norm > 1e-6 else np.array([1.0, 0.0])


def exo_faces(parent, mol_h, mapping, alias=None):
    """{M atom: 'a'|'b'} for the exocyclic neighbours of every specified stereocentre on a ring atom of `parent`.

    `mapping` maps parent locants to atoms of `mol_h` (hydrogens explicit). The faces follow the parent's drawing;
    `alias(locant, atom)` names the parent atom that a modified neighbour stands for."""
    sign = orientation_sign(parent)
    return None if sign is None else _faces(parent, mol_h, mapping, sign, alias)


def _faces(parent, mol_h, mapping, sign, alias=None):
    layout = parent.layout()
    if layout is None:
        return None
    mapped = {atom: loc for loc, atom in mapping.items()}
    faces = {}
    for loc, atom_idx in mapping.items():
        if loc not in parent.plane:
            continue
        atom = mol_h.GetAtomWithIdx(atom_idx)
        if atom.GetChiralTag() not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW) or atom.GetDegree() != 4:
            continue
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]

        def plane_loc(n):
            other = mapped.get(n)
            if other is not None and frozenset((loc, other)) in parent.plane_bonds:
                return other
            other = alias(loc, n) if alias else None
            if other is not None and frozenset((loc, other)) in parent.plane_bonds:
                return other
            return None

        in_plane = [n for n in neighbors if plane_loc(n) is not None]
        exo = [n for n in neighbors if n not in in_plane]
        if len(in_plane) not in (2, 3) or len(exo) != 4 - len(in_plane):
            continue
        origin = layout[loc]
        plane_xy = {n: layout[plane_loc(n)] for n in in_plane}
        out = _outward(origin, plane_xy.values())
        coords = {}
        for n in in_plane:
            coords[n] = (plane_xy[n][0], plane_xy[n][1], 0.0)
        patterns = [(1,), (-1,)] if len(exo) == 1 else [(1, -1), (-1, 1)]
        target = _tag_sign(atom.GetChiralTag())
        for pattern in patterns:
            for n, z in zip(exo, pattern):
                coords[n] = (origin[0] + 0.5 * out[0], origin[1] + 0.5 * out[1], 0.8 * z)
            positions = [coords[n] for n in neighbors]
            volume = _volume_sign((origin[0], origin[1], 0.0), positions)
            if _tag_sign(_TAG_OF_VOLUME[volume]) == target:
                for n, z in zip(exo, pattern):
                    faces[n] = "b" if z * sign > 0 else "a"
                break
    return faces


@lru_cache(maxsize=None)
def get_parent(name):
    return Parent(name)


@lru_cache(maxsize=None)
def parent_h(name):
    return Chem.AddHs(get_parent(name).mol)


def orientation_sign(parent):
    """+1 when the layout's upper face is the drawing's β face, -1 when mirrored, None without a usable layout."""
    if not hasattr(parent, "_sign"):
        parent._sign = _anchor_sign(parent)
    return parent._sign


def _anchor_sign(parent):
    if parent.layout() is None:
        return None
    locant, face = parent.anchor
    mol_h = parent_h(parent.name)
    mapping = {loc: idx for idx, loc in parent.loc_of.items()}
    faces = _faces(parent, mol_h, mapping, 1)
    atom = mol_h.GetAtomWithIdx(parent.idx_of[locant])
    exo = [n.GetIdx() for n in atom.GetNeighbors() if frozenset((locant, parent.loc_of.get(n.GetIdx(), "?"))) not in parent.plane_bonds]
    if not faces or len(exo) != 1 or exo[0] not in faces:
        return None
    return 1 if faces[exo[0]] == face else -1


def ez_relation(mol, bond, near_a, near_b):
    """'trans' or 'cis' between neighbour `near_a` of one bond atom and `near_b` of the other, None if unspecified."""
    state = bond.GetStereo()
    if state not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOCIS):
        return None
    ends = list(bond.GetStereoAtoms())
    if len(ends) != 2:
        return None
    first, second = (ends[0], ends[1])
    if mol.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), first) is None:
        first, second = second, first
    trans = state in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
    if near_a != first:
        trans = not trans
    if near_b != second:
        trans = not trans
    return "trans" if trans else "cis"
