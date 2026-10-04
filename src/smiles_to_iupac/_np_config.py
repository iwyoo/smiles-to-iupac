"""Configuration of a natural product on its parent: α/β faces, CIP descriptors and E/Z (P-101.2.6, P-101.6.2)."""

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure
from ._np_core import FACES, center_signature, exo_faces, ez_relation, loc_key, parent_h


@dataclass
class Config:
    parent: list = field(default_factory=list)
    side: list = field(default_factory=list)
    faces: dict = field(default_factory=dict)
    exo: dict = field(default_factory=dict)
    hfaces: dict = field(default_factory=dict)
    implied_total: int = 0
    implied_cited: int = 0


def _alias(skel, parent, loc, label):
    """The parent atom a neighbour of `loc` stands for after 'homo' insertions or a 'nor' removal."""
    if label.startswith("h") and label in skel.adj:
        previous, current = loc, label
        while current.startswith("h"):
            step = next((n for n in skel.adj[current] if n != previous), None)
            if step is None:
                return label
            previous, current = current, step
        return current
    if label in parent.adj and label not in parent.adj.get(loc, ()):
        bridge = parent.adj.get(loc, set()) & parent.adj[label]
        removed = [x for x in bridge if x not in skel.adj]
        if removed:
            return removed[0]
    return label


def _natural(parent, loc, mol_h, atom, image, skel=None):
    ph = parent_h(parent.name)
    pidx = parent.idx_of[loc]
    plabels = {n.GetIdx(): parent.loc_of.get(n.GetIdx(), "~H") for n in ph.GetAtomWithIdx(pidx).GetNeighbors()}
    mlabels = {}
    for n in mol_h.GetAtomWithIdx(atom).GetNeighbors():
        label = image.get(n.GetIdx(), "~H")
        if skel is not None and label != "~H":
            label = _alias(skel, parent, loc, label)
        mlabels[n.GetIdx()] = label
    missing = [l for l in plabels.values() if l not in mlabels.values()]
    extra = [l for l in mlabels.values() if l not in plabels.values()]
    if len(missing) == 1 and extra == ["~H"] and skel is not None and missing[0] not in skel.adj:
        mlabels = {n: (missing[0] if l == "~H" else l) for n, l in mlabels.items()}
        plabels = {n: (missing[0] if l == "~H" else l) for n, l in plabels.items()}
        if sorted(plabels.values()) != sorted(mlabels.values()):
            return False
        return center_signature(ph, pidx, {k: v for k, v in plabels.items()}) == center_signature(mol_h, atom, mlabels)
    if sorted(plabels.values()) != sorted(mlabels.values()):
        return False
    return center_signature(ph, pidx, plabels) == center_signature(mol_h, atom, mlabels)


def _cip(mol, atom):
    center = mol.GetAtomWithIdx(atom)
    if not center.HasProp("_CIPCode") or center.GetProp("_CIPCode") not in ("R", "S"):
        raise UnsupportedStructure("a natural-product centre has no CIP R/S label")
    return center.GetProp("_CIPCode")


def configuration(cand, view):
    mol = view.mol
    skel, mapping, parent = cand.skel, cand.mapping, cand.parent
    mapped = set(mapping.values())
    image = {a: loc for loc, a in mapping.items()}
    potential = {e.centeredOn: e for e in Chem.FindPotentialStereo(mol) if e.type == Chem.StereoType.Atom_Tetrahedral}
    specified = {a for a, e in potential.items() if e.specified == Chem.StereoSpecified.Specified}
    config = Config()
    bond_stereo = [
        b for b in mol.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and not b.IsInRing()
        and b.GetBeginAtomIdx() in mapped and b.GetEndAtomIdx() in mapped
        and b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
    ]
    if not (specified & mapped) and not bond_stereo:
        return config
    rdCIPLabeler.AssignCIPLabels(mol)
    for bond in bond_stereo:
        if _implied_bond(cand, view, bond, image):
            continue
        if not bond.HasProp("_CIPCode"):
            raise UnsupportedStructure("a chain double bond has no E/Z label")
        low = min((image[bond.GetBeginAtomIdx()], image[bond.GetEndAtomIdx()]), key=loc_key)
        config.side.append((loc_key(low), f"{low}{bond.GetProp('_CIPCode')}"))
    if not (specified & mapped):
        return config
    mol_h = Chem.AddHs(mol)
    parent_atoms = {loc: a for loc, a in mapping.items() if loc in parent.idx_of}
    def alias(loc, n):
        label = image.get(n)
        if label is None or label in parent.idx_of:
            return None
        other = _alias(skel, parent, loc, label)
        return other if other in parent.idx_of else None

    faces = exo_faces(parent, mol_h, parent_atoms, alias)
    ring = skel.ring_atoms()
    for loc, atom in mapping.items():
        if atom not in potential or loc not in parent.idx_of:
            continue
        in_parent_ring = loc in parent.plane
        exo = [
            n.GetIdx() for n in mol_h.GetAtomWithIdx(atom).GetNeighbors()
            if not (
                (n.GetIdx() in mapped and frozenset((loc, image[n.GetIdx()])) in parent.plane_bonds)
                or (alias(loc, n.GetIdx()) is not None and frozenset((loc, alias(loc, n.GetIdx()))) in parent.plane_bonds)
            )
        ]
        if loc in parent.centers:
            if atom not in specified:
                config.parent.append((loc_key(loc), f"{loc}ξ"))
                continue
            if loc in parent.implied:
                config.implied_total += 1
                if _natural(parent, loc, mol_h, atom, image, skel):
                    continue
                config.implied_cited += 1
            chosen = _exo_choice(exo, mapped, ring, image, mol_h)
            if in_parent_ring and faces is not None and chosen in faces:
                config.parent.append((loc_key(loc), f"{loc}{FACES[faces[chosen]]}"))
            else:
                config.side.append((loc_key(loc), f"{loc}{_cip(mol, atom)}"))
            continue
        if atom not in specified:
            if in_parent_ring:
                for n in exo:
                    if mol_h.GetAtomWithIdx(n).GetAtomicNum() != 1 and n not in mapped:
                        config.faces[n] = "x"
            continue
        if in_parent_ring and faces is not None and any(n in faces for n in exo):
            for n in exo:
                if n in faces and mol_h.GetAtomWithIdx(n).GetAtomicNum() != 1:
                    config.faces[n] = faces[n]
                elif n in faces:
                    config.hfaces[atom] = faces[n]
        else:
            config.side.append((loc_key(loc), f"{loc}{_cip(mol, atom)}"))
    return config


def _exo_choice(exo, mapped, ring, image, mol_h):
    skeleton = [n for n in exo if n in mapped and image[n] not in ring]
    if skeleton:
        return skeleton[0]
    heavy = [n for n in exo if mol_h.GetAtomWithIdx(n).GetAtomicNum() != 1]
    return heavy[0] if heavy else exo[0]


def _implied_bond(cand, view, bond, image):
    """True when the chain double bond has the geometry the parent implies."""
    parent, mapping = cand.parent, cand.mapping
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    la, lb = image[a], image[b]
    if la not in parent.idx_of or lb not in parent.idx_of:
        return False
    pb = parent.mol.GetBondBetweenAtoms(parent.idx_of[la], parent.idx_of[lb])
    if pb is None or pb.GetBondTypeAsDouble() != 2.0 or pb.GetStereo() == Chem.BondStereo.STEREONONE:
        return False
    picks = {}
    for end, other in ((la, lb), (lb, la)):
        near = sorted((n for n in parent.adj[end] if n != other), key=loc_key)
        if not near:
            return False
        picks[end] = near[-1]
    if any(picks[end] not in mapping or picks[end] not in cand.skel.adj.get(end, ()) for end in picks):
        return False
    begin_first = parent.idx_of[la] == pb.GetBeginAtomIdx()
    near_begin, near_end = (picks[la], picks[lb]) if begin_first else (picks[lb], picks[la])
    parent_relation = ez_relation(parent.mol, pb, parent.idx_of[near_begin], parent.idx_of[near_end])
    mol_relation = ez_relation(view.mol, bond, mapping[picks[la]], mapping[picks[lb]])
    return parent_relation is not None and parent_relation == mol_relation
