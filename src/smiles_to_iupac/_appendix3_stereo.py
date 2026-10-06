"""Configuration on an Appendix 3 parent (P-101.2.6, P-101.6.2, P-101.8).

The parent name implies the configuration drawn in Appendix 3, stored as the stereo SMILES of `_stereoparents.STEREOPARENTS`.
A molecule is compared with it centre by centre and only what differs is cited:

- on a planar ring system a ring centre that differs from the parent, or has none in the parent, is cited as α/β
  before the parent name (P-101.2.6.1.1); a substituent on a CH2 of the parent carries α/β before its own locant
  (P-101.7.1). β is the side from which the numbering of the reference ring runs anticlockwise in the drawing;
- a centre that is not a ring atom of a planar system, or whose arrangement is not planar, is cited with CIP R/S
  (P-101.2.6.1.3, P-101.8.4); a parent drawn in perspective (bridged) implies no configuration here, so all its
  centres are cited that way;
- chain centres and double bonds are cited with CIP R/S and E/Z; the chain double bonds of carotenoids and retinoids
  with cis/trans (P-101.6.2, P-101.6.3)."""

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler, rdDepictor
from rdkit.Geometry import Point3D

from ._common import UnsupportedStructure

_TAGS = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
_FACE = {"a": "α", "b": "β"}


@dataclass
class StereoResult:
    front: str = ""
    parent: str = ""
    faces: dict = field(default_factory=dict)
    outside: list = field(default_factory=list)


class ParentStereo:
    """The implied configuration of one parent and the drawing frame that gives α and β a meaning."""

    def __init__(self, skeleton, smiles, ref, anticlockwise, label_fn, sort_key):
        mol = Chem.MolFromSmiles(smiles)
        self.labels = {a.GetIdx(): label_fn(a.GetAtomMapNum()) for a in mol.GetAtoms()}
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        self.mol = mol
        self.molh = Chem.AddHs(mol)
        self.atom_of = {label: idx for idx, label in self.labels.items()}
        self.bonds = {}
        for element in Chem.FindPotentialStereo(mol):
            if element.type == Chem.StereoType.Bond_Double and element.specified == Chem.StereoSpecified.Specified:
                bond = mol.GetBondWithIdx(element.centeredOn)
                self.bonds[frozenset((self.labels[bond.GetBeginAtomIdx()], self.labels[bond.GetEndAtomIdx()]))] = (bond, element)
        self.ring_labels = {skeleton.labels[i] for i in skeleton.ring_atoms}
        self.ref = [label for label in ref.split(",") if label]
        self.fused_labels = _fused_labels(skeleton, self.ref)
        self.planar = bool(self.ref)
        self.layout = _layout(skeleton)
        self.beta_sign = _beta_sign(self.layout, self.ref, anticlockwise) if self.planar else 0
        self.sort_key = sort_key

    def specified(self, label):
        index = self.atom_of.get(label)
        return index is not None and self.molh.GetAtomWithIdx(index).GetChiralTag() in _TAGS


def _fused_labels(skeleton, ref):
    """Labels of the rings joined to the reference ring through shared bonds: a ring that is only spiro-linked does
    not lie in the same plane as far as α and β are concerned."""
    rings = [{skeleton.labels[i] for i in ring} for ring in skeleton.query.GetRingInfo().AtomRings()]
    included = [ring for ring in rings if set(ref) <= ring]
    grew = True
    while grew:
        grew = False
        for ring in rings:
            if ring not in included and any(len(ring & other) >= 2 for other in included):
                included.append(ring)
                grew = True
    return set().union(*included) if included else set()


def _layout(skeleton):
    query = Chem.Mol(skeleton.query)
    rdDepictor.Compute2DCoords(query)
    conf = query.GetConformer()
    return {
        skeleton.labels[a.GetIdx()]: (conf.GetAtomPosition(a.GetIdx()).x, conf.GetAtomPosition(a.GetIdx()).y)
        for a in query.GetAtoms()
    }


def _beta_sign(layout, ref, anticlockwise):
    """+1 when β is the viewer side of the planar layout, else -1."""
    area = 0.0
    for i, label in enumerate(ref):
        x1, y1 = layout[label]
        x2, y2 = layout[ref[(i + 1) % len(ref)]]
        area += x1 * y2 - x2 * y1
    return 1 if (area > 0) == anticlockwise else -1


def _tag_for(molh, center, placement):
    """The chiral tag RDKit assigns at `center` when its neighbours sit at `placement` ({atom: (x, y, z)})."""
    copy = Chem.Mol(molh)
    conf = Chem.Conformer(copy.GetNumAtoms())
    conf.Set3D(True)
    for idx in range(copy.GetNumAtoms()):
        conf.SetAtomPosition(idx, Point3D(100.0 + 10.0 * idx, 0.0, 0.0))
    conf.SetAtomPosition(center, Point3D(0.0, 0.0, 0.0))
    for idx, (x, y, z) in placement.items():
        conf.SetAtomPosition(idx, Point3D(x, y, z))
    copy.RemoveAllConformers()
    copy.AddConformer(conf, assignId=True)
    for atom in copy.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignAtomChiralTagsFromStructure(copy)
    return copy.GetAtomWithIdx(center).GetChiralTag()


def _unit(dx, dy):
    length = (dx * dx + dy * dy) ** 0.5 or 1.0
    return dx / length, dy / length


def faces_at(molh, center, ring_neighbors, vectors, beta_sign):
    """{external atom: 'a'|'b'} at the ring atom `center`; `vectors`: ring neighbour -> layout vector from the centre.
    None when the arrangement is not a plain ring atom with one or two exocyclic neighbours."""
    tag = molh.GetAtomWithIdx(center).GetChiralTag()
    if tag not in _TAGS:
        return None
    externals = [n.GetIdx() for n in molh.GetAtomWithIdx(center).GetNeighbors() if n.GetIdx() not in ring_neighbors]
    units = {n: _unit(*vectors[n]) for n in ring_neighbors}
    placement = {n: (u[0], u[1], 0.0) for n, u in units.items()}
    side = lambda z: "b" if z * beta_sign > 0 else "a"
    if len(ring_neighbors) == 3 and len(externals) == 1:
        for z in (1.0, -1.0):
            trial = {**placement, externals[0]: (0.0, 0.0, z)}
            if _tag_for(molh, center, trial) == tag:
                return {externals[0]: side(z)}
        return None
    if len(ring_neighbors) == 2 and len(externals) == 2:
        first, second = units.values()
        wx, wy = -(first[0] + second[0]), -(first[1] + second[1])
        if abs(wx) + abs(wy) < 1e-6:
            wx, wy = -first[1], first[0]
        wx, wy = _unit(wx, wy)
        for z in (1.0, -1.0):
            trial = {
                **placement,
                externals[0]: (wx * 0.5, wy * 0.5, 0.8 * z),
                externals[1]: (wx * 0.5, wy * 0.5, -0.8 * z),
            }
            if _tag_for(molh, center, trial) == tag:
                return {externals[0]: side(z), externals[1]: side(-z)}
        return None
    return None


def _parity(order):
    return sum(1 for i in range(len(order)) for j in range(i + 1, len(order)) if order[i] > order[j]) % 2


def same_configuration(molh, atom, label_of, parent, label):
    """True when the centre at `atom` has the configuration of the parent's centre `label`, False when inverted, None
    when the two centres cannot be paired."""
    center = molh.GetAtomWithIdx(atom)
    parent_center = parent.molh.GetAtomWithIdx(parent.atom_of[label])
    tag_m, tag_p = center.GetChiralTag(), parent_center.GetChiralTag()
    if tag_m not in _TAGS or tag_p not in _TAGS:
        return None
    parent_neighbors = [n.GetIdx() for n in parent_center.GetNeighbors()]
    parent_external = [i for i in parent_neighbors if i not in parent.labels]
    mol_external = [n.GetIdx() for n in center.GetNeighbors() if n.GetIdx() not in label_of]
    if len(parent_external) != len(mol_external) or len(parent_external) > 1:
        return None
    pairing = []
    for n in center.GetNeighbors():
        index = n.GetIdx()
        target = parent.atom_of.get(label_of[index]) if index in label_of else (parent_external[0] if parent_external else None)
        if target is None:
            return None
        pairing.append(parent_neighbors.index(target))
    return (tag_m == tag_p) == (_parity(pairing) == 0)


def _cip(labelled, kind, index):
    item = labelled.GetAtomWithIdx(index) if kind == "atom" else labelled.GetBondWithIdx(index)
    if not item.HasProp("_CIPCode"):
        raise UnsupportedStructure("a stereo element of the skeleton has no CIP label")
    return item.GetProp("_CIPCode")


def describe(mol, skeleton, mapping, sort_key, parent):
    """The stereodescriptors of `mol` on `skeleton`. Specified elements outside the skeleton are returned in
    `outside` for the caller to place."""
    elements = [e for e in Chem.FindPotentialStereo(mol) if e.specified == Chem.StereoSpecified.Specified]
    if not elements:
        return StereoResult()
    label_of = {atom: label for label, atom in mapping.items()}
    labelled = Chem.Mol(mol)
    if not mol.HasProp("cip_assigned"):
        rdCIPLabeler.AssignCIPLabels(labelled)
    molh = Chem.AddHs(mol)
    result = StereoResult()
    marks, rs, bonds, chain_marks = [], [], [], []
    for element in elements:
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            atom = element.centeredOn
            label = label_of.get(atom)
            if label is None or label.startswith("_"):
                result.outside.append(("atom", atom, _cip(labelled, "atom", atom)))
                continue
            cip = _cip(labelled, "atom", atom)
            _center(molh, label_of, label, atom, cip, parent, result, marks, rs)
        elif element.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(element.centeredOn)
            ends = [label_of.get(bond.GetBeginAtomIdx()), label_of.get(bond.GetEndAtomIdx())]
            if None in ends or any(e.startswith("_") for e in ends):
                result.outside.append(("bond", element.centeredOn, _cip(labelled, "bond", element.centeredOn)))
                continue
            relation = _chain_relation(mol, skeleton, label_of, element, bond)
            if relation is not None:
                chain_marks.append((min(ends, key=sort_key), relation))
                continue
            if _implied_double_bond(mol, label_of, bond, element, parent, ends, sort_key) is True:
                continue
            bonds.append((min(ends, key=sort_key), _cip(labelled, "bond", element.centeredOn)))
        else:
            raise UnsupportedStructure("this kind of stereo element is not supported on a natural-product skeleton")
    if marks:
        marks.sort(key=lambda m: sort_key(m[0]))
        result.parent = ",".join(f"{label}{symbol}" for label, symbol in marks) + "-"
    cites = sorted(rs + bonds, key=lambda item: sort_key(item[0]))
    if cites:
        result.front = "(" + ",".join(f"{label}{code}" for label, code in cites) + ")-"
    if chain_marks:
        result.front = _cis_trans_text(skeleton, chain_marks, sort_key) + result.front
    return result


def _across(mol, bond, element, near, far):
    """'cis' or 'trans' between the neighbour `near` of one end of a stereo double bond and `far` of the other."""
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    controlling = list(element.controllingAtoms)
    valid = mol.GetNumAtoms()
    if controlling[0] >= valid or controlling[2] >= valid:
        return None
    first_end = a if mol.GetBondBetweenAtoms(controlling[0], a) is not None else b
    sides = {first_end: controlling[0:2], (b if first_end == a else a): controlling[2:4]}
    end_of_near = a if near in sides[a] else b if near in sides[b] else None
    if end_of_near is None or far not in sides[b if end_of_near == a else a]:
        return None
    base_near = sides[end_of_near][0]
    base_far = sides[b if end_of_near == a else a][0]
    flips = (near != base_near) + (far != base_far)
    cis = element.descriptor == Chem.StereoDescriptor.Bond_Cis
    return "cis" if cis ^ (flips % 2 == 1) else "trans"


def _chain_relation(mol, skeleton, label_of, element, bond):
    """'cis' or 'trans' of the main-chain neighbours across a polyene double bond of a carotenoid or retinoid
    (P-101.6.2), else None."""
    chain = skeleton.polyene
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    la, lb = label_of.get(a), label_of.get(b)
    if la not in chain or lb not in chain:
        return None
    ia, ib = chain.index(la), chain.index(lb)
    if abs(ia - ib) != 1:
        return None
    low, high = sorted((ia, ib))
    if low == 0 or high + 1 >= len(chain):
        return None
    atom_of = {label: atom for atom, label in label_of.items()}
    return _across(mol, bond, element, atom_of[chain[low - 1]], atom_of[chain[high + 1]])


def _implied_double_bond(mol, label_of, bond, element, parent, ends, sort_key):
    """True when the geometry of a skeleton double bond is the one the parent implies, False when it is the other,
    None when the parent does not specify it."""
    entry = parent.bonds.get(frozenset(ends)) if parent else None
    if entry is None:
        return None
    parent_bond, parent_element = entry
    atom_of = {label: atom for atom, label in label_of.items()}
    picks = []
    for end in ends:
        atom = atom_of[end]
        other = atom_of[ends[1] if end == ends[0] else ends[0]]
        neighbors = sorted(
            (label_of[n.GetIdx()] for n in mol.GetAtomWithIdx(atom).GetNeighbors() if n.GetIdx() != other and n.GetIdx() in label_of),
            key=sort_key,
        )
        if not neighbors:
            return None
        picks.append(neighbors[0])
    mine = _across(mol, bond, element, atom_of[picks[0]], atom_of[picks[1]])
    implied = _across(
        parent.mol, parent_bond, parent_element, parent.atom_of[picks[0]], parent.atom_of[picks[1]]
    )
    return None if mine is None or implied is None else mine == implied


def _cis_trans_text(skeleton, marks, sort_key):
    """'all-trans-' when every chain double bond is specified and trans (P-101.6.3), else the cis ones by locant;
    the others are the usual trans."""
    cis = sorted((label for label, relation in marks if relation == "cis"), key=sort_key)
    if not cis:
        return "all-trans-" if len(marks) == skeleton.polyene_bonds else ""
    prefix = {1: "", 2: "di-", 3: "tri-"}.get(len(cis), "")
    return f"{','.join(cis)}-{prefix}cis-"


def _center(molh, label_of, label, atom, cip, parent, result, marks, rs):
    if parent is None:
        rs.append((label, cip))
        return
    has_center = parent.specified(label)
    verdict = same_configuration(molh, atom, label_of, parent, label) if has_center else None
    if verdict is True:
        return
    if not (parent.planar and label in parent.fused_labels):
        rs.append((label, cip))
        return
    ring_neighbors = [n.GetIdx() for n in molh.GetAtomWithIdx(atom).GetNeighbors() if label_of.get(n.GetIdx()) in parent.ring_labels]
    x0, y0 = parent.layout[label]
    vectors = {n: (parent.layout[label_of[n]][0] - x0, parent.layout[label_of[n]][1] - y0) for n in ring_neighbors}
    found = faces_at(molh, atom, ring_neighbors, vectors, parent.beta_sign)
    if found is None:
        rs.append((label, cip))
        return
    chain = [e for e in found if e in label_of]
    if len(ring_neighbors) == 3 or chain or verdict is False:
        key = chain[0] if chain else next(iter(found))
        marks.append((label, _FACE[found[key]]))
        return
    for external, face in found.items():
        if molh.GetAtomWithIdx(external).GetAtomicNum() != 1:
            result.faces[external] = _FACE[face]


def deviations(mol, mapping, parent):
    """How many specified centres of `mol` have the opposite configuration to the parent's, for choosing between
    parents that differ only in configuration (dammarane and protostane)."""
    if parent is None:
        return 0
    molh = Chem.AddHs(mol)
    label_of = {atom: label for label, atom in mapping.items()}
    count = 0
    for label, atom in mapping.items():
        if parent.specified(label) and molh.GetAtomWithIdx(atom).GetChiralTag() in _TAGS:
            count += same_configuration(molh, atom, label_of, parent, label) is False
    return count
