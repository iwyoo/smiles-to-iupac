"""Parent hydrides whose bonding number differs from the standard one (P-21.1.2.1) and homogeneous chalcogen chains
(P-21.2.2, P-21.2.4): the symbol λn is affixed to the hydride name ('λ5-arsane', '2λ4-trisulfane'); a mononuclear center
may carry single-bonded carbon substituents ('dimethyl-λ4-sulfane')."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency
from ._hetero_prefixes import AMINE_PATTERN
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, name_branch

_HYDRIDES = {
    5: ("borane", 3), 13: ("alumane", 3), 31: ("gallane", 3), 49: ("indigane", 3), 81: ("thallane", 3),
    14: ("silane", 4), 32: ("germane", 4), 50: ("stannane", 4), 82: ("plumbane", 4),
    7: ("azane", 3), 15: ("phosphane", 3), 33: ("arsane", 3), 51: ("stibane", 3), 83: ("bismuthane", 3),
    8: ("oxidane", 2), 16: ("sulfane", 2), 34: ("selane", 2), 52: ("tellane", 2),
    9: ("fluorane", 1), 17: ("chlorane", 1), 35: ("bromane", 1), 53: ("iodane", 1),
}
_CHALCOGENS = (8, 16, 34, 52)
_HYDROXY_ON_CARBON = Chem.MolFromSmarts("[OX2H1][#6]")
_SENIOR_TO_ALCOHOL = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=O)[OX2H1]", "[#16,#34,#52;X3,X4](=O)[OX2H1]", "[CX3](=O)[OX2][#6]", "[CX3](=[O,S,Se,Te])[NX3]",
        "[CX2]#[NX1]", "[CX3H1](=O)", "[#6][CX3](=O)[#6]", "[CX3](=O)[F,Cl,Br,I]",
    )
]
_LOWER_PRINCIPAL = [Chem.MolFromSmarts("[SX2H1][#6]"), AMINE_PATTERN]


def _plain(atom):
    return not atom.GetFormalCharge() and not atom.GetIsotope()


def _center(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    bridges = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8 and a.GetDegree() == 2 and _substituent_atom(a)}
    centers = [
        a for a in mol.GetAtoms()
        if a.GetAtomicNum() in _HYDRIDES and not _terminal_standard(a) and a.GetIdx() not in bridges and not _oxo_acyl(a, bridges)
    ]
    if len(centers) != 1 or not _plain(centers[0]) or centers[0].IsInRing() or centers[0].GetIsAromatic():
        return None
    center = centers[0]
    excess = center.GetTotalValence() - _HYDRIDES[center.GetAtomicNum()][1]
    if excess == 0 or excess % 2 or (excess < 0 and center.GetDegree()):
        return None
    if any(b.GetBondTypeAsDouble() != 1.0 or not _substituent_atom(b.GetOtherAtom(center)) for b in center.GetBonds()):
        return None
    if _outranked(mol, center):
        return None
    return center


def _oxo_acyl(atom, bridges):
    """The sulfur of a sulfonyl-type acyl group joined to a centre through an oxygen of `bridges`."""
    return any(n.GetIdx() in bridges for n in atom.GetNeighbors()) and any(
        b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds()
    )


def _terminal_standard(atom):
    return (
        atom.GetDegree() == 1
        and atom.GetAtomicNum() in (8, 9, 17, 35, 53)
        and atom.GetTotalValence() == _HYDRIDES[atom.GetAtomicNum()][1]
    )


def _outranked(mol, center):
    """A characteristic group on a carbon substituent that is senior to the centre's own hydroxy groups (or to nothing,
    when it has none) makes the carbon skeleton the parent and the centre a prefix (P-41, P-44.1.1)."""
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_ALCOHOL):
        return True
    hydroxy_on_centre = sum(
        n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1 and n.GetDegree() == 1 for n in center.GetNeighbors()
    )
    hydroxy_on_carbon = len(mol.GetSubstructMatches(_HYDROXY_ON_CARBON))
    if hydroxy_on_carbon > hydroxy_on_centre:
        return True
    return not hydroxy_on_centre and any(mol.HasSubstructMatch(query) for query in _LOWER_PRINCIPAL)


def _substituent_atom(atom):
    """A carbon, halogen, hydroxy oxygen or alkoxy/acyloxy oxygen bonded to a hydride centre of nonstandard bonding
    number: iodine is not a pseudoester element, so an acyl-O-iodane is no ester (P-65.6.3.1.2)."""
    z = atom.GetAtomicNum()
    if z == 6:
        return True
    if atom.GetFormalCharge():
        return False
    if z == 8 and atom.GetDegree() == 2 and not atom.GetTotalNumHs() and not atom.IsInRing():
        return any(n.GetAtomicNum() in (6, 16) for n in atom.GetNeighbors())
    return atom.GetDegree() == 1 and (z in (9, 17, 35, 53) or (z == 8 and atom.GetTotalNumHs() == 1))


def _chalcogen_chain(mol):
    atoms = list(mol.GetAtoms())
    if len(atoms) < 2 or len({a.GetAtomicNum() for a in atoms}) != 1 or atoms[0].GetAtomicNum() not in _CHALCOGENS:
        return None
    if mol.GetRingInfo().NumRings() or any(not _plain(a) or a.GetDegree() > 2 for a in atoms):
        return None
    if len(Chem.GetMolFrags(mol)) != 1 or any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()):
        return None
    if any(a.GetTotalValence() % 2 or a.GetNumRadicalElectrons() for a in atoms):
        return None
    return atoms


def has_nonstandard_hydride_shape(mol) -> bool:
    if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        return False
    if _center(mol) is not None:
        return True
    return _chalcogen_chain(mol) is not None


def _chain_name(mol, atoms):
    graph = adjacency(mol)
    start = next(a.GetIdx() for a in atoms if a.GetDegree() <= 1)
    order = [start]
    while len(order) < len(atoms):
        order.append(next(n for n in graph[order[-1]] if n not in order))
    stem = _HYDRIDES[atoms[0].GetAtomicNum()][0]
    best = None
    for sequence in (order, order[::-1]):
        locants = tuple(i + 1 for i, a in enumerate(sequence) if mol.GetAtomWithIdx(a).GetTotalValence() != 2)
        if best is None or locants < best[0]:
            best = (locants, sequence)
    labels = ",".join(f"{i}λ{mol.GetAtomWithIdx(best[1][i - 1]).GetTotalValence()}" for i in best[0])
    return f"{labels + '-' if labels else ''}{multiplying_prefix(len(atoms))}{stem}"


def name_nonstandard_hydride(mol) -> str:
    center = _center(mol)
    if center is None:
        chain = _chalcogen_chain(mol)
        if chain is None:
            raise UnsupportedStructure("not a nonstandard-bonding-number hydride")
        return _chain_name(mol, chain)
    if any(a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified structures are not supported yet")
    graph = adjacency(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    hydroxyls = [n for n in graph[center.GetIdx()] if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1]
    halogens = {a.GetIdx(): HALOGEN_PREFIXES[a.GetAtomicNum()] for a in mol.GetAtoms() if a.GetAtomicNum() in HALOGEN_PREFIXES}
    entries = [
        name_branch(graph, n, center.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
        for n in graph[center.GetIdx()]
        if n not in hydroxyls
    ]
    stem, _ = _HYDRIDES[center.GetAtomicNum()]
    prefixes = f"{format_mononuclear_prefixes(entries)}-" if entries else ""
    suffix = ""
    if hydroxyls:
        suffix = {1: "ol", 2: "diol", 3: "triol"}[len(hydroxyls)]
        stem = stem[:-1] if len(hydroxyls) == 1 else stem
    return f"{prefixes}λ{center.GetTotalValence()}-{stem}{suffix}"
