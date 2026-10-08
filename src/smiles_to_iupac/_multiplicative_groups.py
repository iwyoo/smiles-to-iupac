"""Seniority-ranked functional-group vocabulary used to decide whether a
multiplicative name (P-15.3, P-51.3) is the PIN: the principal characteristic
group class (P-41) must sit entirely inside the multiplied parent structures.
Anything outside the vocabulary makes `classify` return None.
"""

from dataclasses import dataclass

from rdkit import Chem

SUFFIX_RANKS = {
    "carboxylic_acid": 10,
    "sulfonic_acid": 11,
    "sulfinic_acid": 12,
    "ester": 21,
    "sulfonate_ester": 21,
    "acyl_halide": 22,
    "sulfonyl_halide": 22,
    "amide": 23,
    "sulfonamide": 23,
    "nitrile": 25,
    "aldehyde": 26,
    "ketone": 27,
    "alcohol": 28,
    "thiol": 28,
    "hydroperoxide": 29,
    "amine": 30,
    "imine": 32,
}

JUNIOR_RANK = 80

_ALLOWED_ELEMENTS = {1, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 31, 32, 33, 34, 35, 49, 50, 51, 52, 53, 81, 82, 83}

_NOT_CARBONYL = "[#6;!$([#6]=[O,S,N])]"
_CHALCOGEN2 = "[OX2,SX2,SeX2,TeX2;!R]"

_CLASS_PATTERNS = [
    ("carboxylic_acid", "[CX3](=O)[OX2H1]", 0),
    ("sulfonic_acid", "[SX4](=O)(=O)[OX2H1]", 0),
    ("sulfinic_acid", "[SX3](=O)[OX2H1]", 0),
    ("ester", "[CX3](=O)[OX2;!R][#6]", 0),
    ("sulfonate_ester", "[SX4](=O)(=O)[OX2;!R][#6]", 0),
    ("acyl_halide", "[CX3](=O)[F,Cl,Br,I]", 0),
    ("sulfonyl_halide", "[SX4](=O)(=O)[F,Cl,Br,I]", 0),
    ("amide", "[CX3](=O)[NX3;!R]", 0),
    ("sulfonamide", "[SX4](=O)(=O)[NX3;!R]", 0),
    ("nitrile", "[CX2]#[NX1]", 0),
    ("aldehyde", "[CX3H1](=O)[#6]", 0),
    ("ketone", "[#6][CX3](=O)[#6]", 1),
    ("alcohol", f"[OX2H1;!R]{_NOT_CARBONYL}", 0),
    ("thiol", f"[SX2H1;!R]{_NOT_CARBONYL}", 0),
    ("hydroperoxide", "[OX2H1][OX2;!R][#6]", 0),
    ("amine", f"[NX3;!R;!$(N~[!#6;!#1]);!$(N[CX3]=[O,S,N]);!$(N-[#6]#*)]{_NOT_CARBONYL}", 0),
    ("imine", "[CX3;!R;!$(C-[O,N,S,Se,Te])]=[NX2;!R;!$(N-[!#6;!#8])]", 0),
]

_JUNIOR_PATTERNS = [
    f"{_NOT_CARBONYL}{_CHALCOGEN2}{_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[OX2;!R][OX2;!R]{_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[SX2,SeX2,TeX2;!R][SX2,SeX2,TeX2;!R]{_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[SX3;!R](=O){_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[SX4;!R](=O)(=O){_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[PX4,AsX4,SbX4;!R](=O)[OX2H1,#6]",
    "[#6][N+](=O)[O-]",
    "[#6][NX2]=O",
    "[#6][NX2;!R]=[NX2;!R][#6]",
    "[#6;!R]=[NX2;!R][NX3H2;!R]",
]

_CLASS_QUERIES = [(name, Chem.MolFromSmarts(smarts), anchor) for name, smarts, anchor in _CLASS_PATTERNS]
_JUNIOR_QUERIES = [Chem.MolFromSmarts(smarts) for smarts in _JUNIOR_PATTERNS]
_THIOKETONE_QUERY = Chem.MolFromSmarts("[#6][CX3](=[SX1,SeX1,TeX1])[#6]")


@dataclass(frozen=True)
class Group:
    name: str
    rank: int
    anchor: int
    atoms: frozenset
    ring_atom: object


# P-44.1.2.2 / P-41 class order of ring heteroatoms other than N: O S Se Te P As Sb Bi Si Ge Sn Pb B Al Ga In Tl
_NON_NITROGEN_RANK = {z: i for i, z in enumerate((8, 16, 34, 52, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81))}


def ring_seniority_key(mol, ring_atoms):
    """P-44.2.1 (a)-(g) then P-44.4.1.1 sort key of a ring system: a smaller
    key is the senior parent structure."""
    atoms = set(ring_atoms)
    heteroatoms = [mol.GetAtomWithIdx(a).GetAtomicNum() for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    has_nitrogen = 7 in heteroatoms
    order = _NON_NITROGEN_RANK
    best_non_nitrogen = min((order.get(z, len(order)) for z in heteroatoms if z != 7), default=len(order)) if heteroatoms and not has_nitrogen else 0
    rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if set(ring) <= atoms)
    kekulized = Chem.Mol(mol)
    Chem.Kekulize(kekulized, clearAromaticFlags=True)
    multiple = sum(
        1
        for b in kekulized.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms and b.GetBondTypeAsDouble() > 1
    )
    earlier = tuple(-sum(1 for z in heteroatoms if z == e) for e in (*_NON_NITROGEN_RANK, 7))
    return (
        0 if heteroatoms else 1,
        0 if has_nitrogen else 1,
        best_non_nitrogen,
        -rings,
        -len(atoms),
        -len(heteroatoms),
        earlier,
        -multiple,
    )


def _single_bonded_linker(atom):
    return atom.GetDegree() >= 2 and all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds())


def _hydroxy_on_nitrogen(atom):
    """The oxygen of an N-OH group, cited as 'hydroxy' on a linking nitrogen (hydroxyazanediyl, P-68.3.1.1.1.5)."""
    return (
        atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and atom.GetNeighbors()[0].GetAtomicNum() == 7
    )


def classify(mol):
    """Return the list of `Group`s found in `mol`, or None when `mol` holds
    an atom, charge, or group this vocabulary can't place (so no claim about
    its principal class can be made)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ELEMENTS or atom.GetIsotope() or atom.GetNumRadicalElectrons():
            return None
    covered = set()
    groups = []
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    for name, query, anchor in _CLASS_QUERIES:
        for match in mol.GetSubstructMatches(query):
            if any(a in ring_atoms for i, a in enumerate(match) if i != anchor and mol.GetAtomWithIdx(a).GetAtomicNum() != 6):
                return None
            anchor_atom = match[anchor]
            ring_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(anchor_atom).GetNeighbors() if n.GetIdx() in ring_atoms]
            if anchor_atom in ring_atoms and name == "ketone":
                ring_atom = anchor_atom
            elif len(ring_neighbors) == 1:
                ring_atom = ring_neighbors[0]
            else:
                ring_atom = None
            group_atoms = frozenset(match)
            if name == "ketone":
                group_atoms = frozenset({match[1], match[2]})
            groups.append(Group(name, SUFFIX_RANKS[name], anchor_atom, group_atoms, ring_atom))
            covered.update(group_atoms)
    junior = _JUNIOR_QUERIES + [_THIOKETONE_QUERY] if any(g.rank < SUFFIX_RANKS["ketone"] for g in groups) else _JUNIOR_QUERIES
    for query in junior:
        for match in mol.GetSubstructMatches(query):
            covered.update(match)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        z = atom.GetAtomicNum()
        if idx in covered or atom.GetIsAromatic():
            continue
        if (
            z in (7, 8, 16, 34, 52)
            and idx not in ring_atoms
            and not _single_bonded_linker(atom)
            and not _hydroxy_on_nitrogen(atom)
        ):
            return None
        if z == 6 and any(
            b.GetBondTypeAsDouble() >= 2 and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34, 52) for b in atom.GetBonds()
        ):
            return None
        if atom.GetFormalCharge() and idx not in covered:
            return None
    return groups
