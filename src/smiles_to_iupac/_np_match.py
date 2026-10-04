"""Embedding of a (modified) stereoparent graph into a molecule."""

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdqueries


class View:
    """Heavy-atom view of a molecule: adjacency, elements, ring atoms and a bond-order-free copy for matching."""

    def __init__(self, mol):
        self.mol = mol
        self.adj = {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors()} for a in mol.GetAtoms()}
        self.elem = {a.GetIdx(): a.GetSymbol() for a in mol.GetAtoms()}
        self.rings = {i for ring in mol.GetRingInfo().AtomRings() for i in ring}
        flat = Chem.RWMol(mol)
        for atom in flat.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        for bond in flat.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
            bond.SetStereo(Chem.BondStereo.STEREONONE)
        self.flat = flat.GetMol()
        kekule = Chem.Mol(mol)
        Chem.Kekulize(kekule, clearAromaticFlags=False)
        self.order = {}
        for bond in kekule.GetBonds():
            self.order[frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))] = int(bond.GetBondTypeAsDouble())
        self.aromatic_bonds = {
            frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds() if b.GetIsAromatic()
        }


@dataclass
class Candidate:
    parent: object
    skel: object
    mapping: dict
    replaced: list = field(default_factory=list)
    cyclo: list = field(default_factory=list)

    @property
    def ops(self):
        return self.skel.ops


_ATOMIC_NUMBER = {"C": 6, "N": 7, "O": 8, "S": 16, "Se": 34, "Te": 52, "P": 15, "Si": 14, "B": 5, "As": 33}
_REPLACEABLE = "#6,#7,#8,#16,#34,#52,#15,#14,#5"


def _query(skel):
    labels = list(skel.adj)
    index = {label: i for i, label in enumerate(labels)}
    ring = skel.ring_atoms()
    query = Chem.RWMol()
    for label in labels:
        number = _ATOMIC_NUMBER.get(skel.elem[label], 6)
        smarts = f"[{_REPLACEABLE};R]" if label in ring else f"[#{number}]"
        query.AddAtom(Chem.AtomFromSmarts(smarts))
    for a in labels:
        for b in skel.adj[a]:
            if index[a] < index[b]:
                query.AddBond(index[a], index[b], Chem.BondType.SINGLE)
    return labels, query.GetMol()


def embeddings(skel, view, limit=4000):
    """Candidates for every embedding of `skel` (induced up to extra bonds, which become 'cyclo' bonds)."""
    if len(skel.adj) > len(view.adj):
        return []
    labels, query = _query(skel)
    ring = skel.ring_atoms()
    found = []
    for match in view.flat.GetSubstructMatches(query, uniquify=False, maxMatches=limit):
        mapping = {label: match[i] for i, label in enumerate(labels)}
        replaced = []
        for label, atom in mapping.items():
            if skel.elem[label] != view.elem[atom]:
                if label in ring and atom in view.rings:
                    replaced.append(label)
                else:
                    break
        else:
            image = {atom: label for label, atom in mapping.items()}
            cyclo = []
            for atom, label in image.items():
                for n in view.adj[atom]:
                    other = image.get(n)
                    if other is not None and n > atom and other not in skel.adj[label]:
                        cyclo.append((label, other))
            found.append(Candidate(skel.parent, skel, mapping, replaced, cyclo))
    return found
