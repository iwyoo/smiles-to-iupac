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
    "sulfonamide": 24.6,
    "chalcogen_amide": 24,
    "amidine": 24.5,
    "nitrile": 25,
    "aldehyde": 26,
    "ketone": 27,
    "thione": 27.1,
    "selone": 27.2,
    "tellone": 27.3,
    "alcohol": 28,
    "thiol": 28,
    "hydroperoxide": 29,
    "amine": 30,
    "imine": 32,
}

JUNIOR_RANK = 80

KETONE_LIKE = {"ketone", "thione", "selone", "tellone"}

_ALLOWED_ELEMENTS = {1, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 31, 32, 33, 34, 35, 49, 50, 51, 52, 53, 81, 82, 83}

_NOT_CARBONYL = "[#6;!$([#6]=[O,S,N])]"
_CHALCOGEN2 = "[OX2,SX2,SeX2,TeX2;!R]"

_CLASS_PATTERNS = [
    ("carboxylic_acid", "[CX3](=O)[OX2H1]", 0),
    ("sulfonic_acid", "[SX4](=O)(=O)[OX2H1]", 0),
    ("sulfinic_acid", "[SX3](=O)[OX2H1]", 0),
    ("ester", "[CX3](=O)[OX2;!R][#6]", 0),
    ("ester", "[CX3](=[O,S,Se,Te])[SX2,SeX2,TeX2;!R][#6]", 0),
    ("ester", "[CX3](=[S,Se,Te])[OX2;!R][#6]", 0),
    ("sulfonate_ester", "[SX4](=O)(=O)[OX2;!R][#6]", 0),
    ("acyl_halide", "[CX3](=O)[F,Cl,Br,I;X1]", 0),
    ("sulfonyl_halide", "[SX4](=O)(=O)[F,Cl,Br,I;X1]", 0),
    ("amide", "[CX3](=O)[NX3;!R]", 0),
    ("chalcogen_amide", "[CX3](=[S,Se,Te])[NX3;!R]", 0),
    ("amidine", "[CX3;!$(C(=N)(-N)-[!#6])](=[NX2H1;!R])[NX3H2;!R]", 0),
    ("sulfonamide", "[SX4](=O)(=O)[NX3;!R]", 0),
    ("nitrile", "[CX2]#[NX1]", 0),
    ("aldehyde", "[CX3H1](=O)[#6]", 0),
    ("ketone", "[#6][CX3](=O)[#6]", 1),
    ("ketone", "[#6X3;R;$([#6]@[#7,#8,#16,#34,#52])]=O", 0),
    ("ketone", "[CX3;!R;$(C-[#7;R])]=O", 0),
    ("thione", "[#6][CX3](=[SX1])[#6]", 1),
    ("thione", "[#6X3;R;$([#6]@[#7,#8,#16,#34,#52])]=[SX1]", 0),
    ("selone", "[#6][CX3](=[SeX1])[#6]", 1),
    ("selone", "[#6X3;R;$([#6]@[#7,#8,#16,#34,#52])]=[SeX1]", 0),
    ("tellone", "[#6][CX3](=[TeX1])[#6]", 1),
    ("tellone", "[#6X3;R;$([#6]@[#7,#8,#16,#34,#52])]=[TeX1]", 0),
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
    f"{_NOT_CARBONYL}[SX3,SeX3,TeX3;!R](=O){_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[SX4,SeX4,TeX4;!R](=O)(=O){_NOT_CARBONYL}",
    f"{_NOT_CARBONYL}[PX4,AsX4,SbX4;!R](=[O,S,Se,Te])[OX2,NX3H2,F,Cl,Br,I,#6]",
    "[#6,#7][N+](=O)[O-]",
    "[#6,#7][NX2]=O",
    "[#6][NX2;!R]=[NX2;!R][#6]",
    "[#6;!R]=[NX2;!R][NX3H2;!R]",
]

_CLASS_QUERIES = [(name, Chem.MolFromSmarts(smarts), anchor) for name, smarts, anchor in _CLASS_PATTERNS]
_AMMONIUM_QUERY = ("amine", Chem.MolFromSmarts(f"[NX4+;!R;!$(N~[!#6;!#1])]{_NOT_CARBONYL}"), 0)
_JUNIOR_QUERIES = [Chem.MolFromSmarts(smarts) for smarts in _JUNIOR_PATTERNS]
_THIOKETONE_QUERY = Chem.MolFromSmarts("[#6][CX3](=[SX1,SeX1,TeX1])[#6]")
_HIDDEN_THIOAMIDE_QUERY = Chem.MolFromSmarts("[CX3;!R;$(C-[#7;R])]=[SX1,SeX1,TeX1]")


@dataclass(frozen=True)
class Group:
    name: str
    rank: int
    anchor: int
    atoms: frozenset
    ring_atom: object


def _single_bonded_linker(atom):
    if atom.GetAtomicNum() in (16, 34, 52) and atom.GetFormalCharge() == 0:
        # sulfinyl and sulfonyl linkers and their Se and Te analogues: two single bonds to the units and terminal =O atoms (P-63.6)
        oxo = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8]
        singles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 1.0]
        if oxo and len(singles) == 2 and len(oxo) + len(singles) == atom.GetDegree():
            return all(b.GetOtherAtom(atom).GetDegree() == 1 for b in oxo)
    return atom.GetDegree() >= 2 and all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds())


def _hydrazinediylidene_nitrogen(atom):
    """A nitrogen of an =N-N= linker whose both nitrogens are doubly bonded to carbon (hydrazinediylidene)."""
    if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge() or atom.GetTotalNumHs():
        return False
    doubles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 2.0]
    singles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 1.0]
    if len(doubles) != 1 or len(singles) != 1 or doubles[0].GetOtherAtom(atom).GetAtomicNum() != 6:
        return False
    partner = singles[0].GetOtherAtom(atom)
    return (
        partner.GetAtomicNum() == 7
        and partner.GetDegree() == 2
        and not partner.GetTotalNumHs()
        and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(partner).GetAtomicNum() == 6 for b in partner.GetBonds())
    )


_MONONUCLEAR_CENTRES = {5, 7, 13, 14, 15, 31, 32, 33, 49, 50, 51, 81, 82, 83}


def _hydroxy_on_nitrogen(atom):
    """The oxygen of a hydroxy group on a linking mononuclear hydride atom, cited as 'hydroxy' (hydroxyazanediyl,
    hydroxystibanediyl; P-68.3.1.1.1.5, P-29.3.1)."""
    return (
        atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and atom.GetNeighbors()[0].GetAtomicNum() in _MONONUCLEAR_CENTRES
    )


def terminal_ligand(atom, centre):
    """('O' | 'S' | 'Se' | 'Te' | 'NH' | 'NNH2', atoms) when `atom` is a terminal =X position of a sulfonyl-type centre
    (P-65.3.2.3 infixes), else None."""
    bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), centre.GetIdx())
    if bond.GetBondTypeAsDouble() != 2.0 or atom.GetFormalCharge():
        return None
    z = atom.GetAtomicNum()
    if atom.GetDegree() == 1:
        symbol = {8: "O", 16: "S", 34: "Se", 52: "Te"}.get(z)
        if symbol:
            return symbol, (atom.GetIdx(),)
        if z == 7 and atom.GetTotalNumHs() == 1:
            return "NH", (atom.GetIdx(),)
    if z == 7 and atom.GetDegree() == 2 and not atom.GetTotalNumHs():
        partner = next(n for n in atom.GetNeighbors() if n.GetIdx() != centre.GetIdx())
        bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), partner.GetIdx())
        if partner.GetAtomicNum() == 6 and bond.GetBondTypeAsDouble() == 1.0:
            return "NH", (atom.GetIdx(),)
        if partner.GetAtomicNum() == 7 and partner.GetDegree() == 1 and partner.GetTotalNumHs() == 2 and bond.GetBondTypeAsDouble() == 1.0:
            return "NNH2", (atom.GetIdx(), partner.GetIdx())
    return None


def chalcogen_linker_ligands(atom):
    """The =X ligands of a sulfonyl- or sulfinyl-type centre (S, Se, Te) with two single bonds and one or two terminal
    =X positions, as [(symbol, atoms)], else None (P-65.3.2.3, P-67.1.4.4.1)."""
    if atom.GetAtomicNum() not in (16, 34, 52) or atom.GetFormalCharge() or atom.IsInRing() or atom.GetIsAromatic():
        return None
    singles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 1.0]
    doubles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 2.0]
    if len(singles) != 2 or len(doubles) not in (1, 2) or len(singles) + len(doubles) != atom.GetDegree():
        return None
    ligands = [terminal_ligand(b.GetOtherAtom(atom), atom) for b in doubles]
    return None if None in ligands else ligands


def classify(mol):
    """Return the list of `Group`s found in `mol`, or None when `mol` holds
    an atom, charge, or group this vocabulary can't place (so no claim about
    its principal class can be made)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ELEMENTS or atom.GetIsotope() or atom.GetNumRadicalElectrons():
            return None
    from ._amine import CATIONIC_AMINES

    covered = set()
    groups = []
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    for name, query, anchor in _CLASS_QUERIES + ([_AMMONIUM_QUERY] if CATIONIC_AMINES.get() else []):
        for match in mol.GetSubstructMatches(query):
            if any(a in ring_atoms for i, a in enumerate(match) if i != anchor and mol.GetAtomWithIdx(a).GetAtomicNum() != 6):
                return None
            anchor_atom = match[anchor]
            ring_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(anchor_atom).GetNeighbors() if n.GetIdx() in ring_atoms]
            if anchor_atom in ring_atoms and name in KETONE_LIKE:
                ring_atom = anchor_atom
            elif len(ring_neighbors) == 1:
                ring_atom = ring_neighbors[0]
            else:
                ring_atom = None
            group_atoms = frozenset(match)
            if name in KETONE_LIKE and len(match) == 4:
                group_atoms = frozenset({match[1], match[2]})
            groups.append(Group(name, SUFFIX_RANKS[name], anchor_atom, group_atoms, ring_atom))
            covered.update(group_atoms)
    junior = _JUNIOR_QUERIES + [_THIOKETONE_QUERY] if any(g.rank < SUFFIX_RANKS["ketone"] for g in groups) else _JUNIOR_QUERIES
    if any(g.rank <= SUFFIX_RANKS["tellone"] for g in groups):
        junior = junior + [_HIDDEN_THIOAMIDE_QUERY]
    for query in junior:
        for match in mol.GetSubstructMatches(query):
            covered.update(match)
    from ._formazan import FORMAZAN_SKELETON

    for chain in mol.GetSubstructMatches(FORMAZAN_SKELETON):
        covered.update(chain)
    for atom in mol.GetAtoms():
        ligands = chalcogen_linker_ligands(atom)
        if ligands is not None:
            covered.add(atom.GetIdx())
            covered.update(a for _, atoms in ligands for a in atoms)
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
            and not _hydrazinediylidene_nitrogen(atom)
        ):
            return None
        if z == 6 and any(
            b.GetBondTypeAsDouble() >= 2
            and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34, 52)
            and not _hydrazinediylidene_nitrogen(b.GetOtherAtom(atom))
            for b in atom.GetBonds()
        ):
            return None
        if atom.GetFormalCharge() and idx not in covered:
            return None
    return groups
