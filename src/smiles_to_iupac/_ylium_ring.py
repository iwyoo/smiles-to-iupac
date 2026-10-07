"""Cationic rings whose centre is a heteroatom with one more skeletal bond than in the neutral heterocycle (P-73.3): the
retained names of Table 7.5 (pyrylium, thiopyrylium, xanthylium ...), the lambda names built on benzopyran and on the
furan family ('1λ4-benzopyran-1-ylium', '3H-1λ4-thiophen-1-ylium') and 'quinolizin-5-ylium'. Each skeleton is a
template whose atoms carry the locants of the neutral parent hydride; the embeddings of the template in the molecule are
the candidate numberings, and the lowest locants go to the prefixes.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    substituent_locant_set_and_citation,
)
from ._substituents import format_substituent_prefixes, name_branch

_INFIX = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}
_PYRAN = {8: "pyran", 16: "thiopyran", 34: "selenopyran", 52: "telluropyran"}
_FIVE = {8: "furan", 16: "thiophen", 34: "selenophen", 52: "tellurophen"}
_XANTHENE = ["10", "4a", "4", "3", "2", "1", "9a", "9", "8a", "8", "7", "6", "5", "10a"]
_BENZO = [1, 2, 3, 4, "4a", 5, 6, 7, 8, "8a"]


def _skeletons(z):
    e = f"[#{z}+]"
    infix = _INFIX.get(z, "")
    skeletons = []
    if z in _INFIX:
        skeletons += [
            (f"{e}1ccccc1", [1, 2, 3, 4, 5, 6], 0, f"{infix}pyrylium"),
            (f"{e}1c2ccccc2cc2ccccc12", _XANTHENE, 0, f"{infix}xanthylium"),
            (f"{e}1cccc2ccccc12", _BENZO, 0, f"1λ4-benzo{_PYRAN[z]}-1-ylium"),
            (f"c1{e}ccc2ccccc12", _BENZO, 1, f"2λ4-benzo{_PYRAN[z]}-2-ylium"),
            (f"{e}1-[CX4]-[#6]=[#6]-[#6]=1", [1, 2, 3, 4, 5], 0, f"2H-1λ4-{_FIVE[z]}-1-ylium"),
            (f"{e}1=[#6]-[CX4]-[#6]=[#6]-1", [1, 2, 3, 4, 5], 0, f"3H-1λ4-{_FIVE[z]}-1-ylium"),
        ]
    if z == 7:
        skeletons.append(("c1ccc[#7+]2ccccc12", [1, 2, 3, 4, 5, 6, 7, 8, 9, "9a"], 4, "5λ5-quinolizin-5-ylium"))
    return skeletons


def _sigma(atom):
    return atom.GetDegree() + atom.GetTotalNumHs()


def has_ylium_ring_shape(mol) -> bool:
    centres = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(centres) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return False
    atom = centres[0]
    return (
        atom.GetFormalCharge() == 1
        and atom.IsInRing()
        and not atom.GetIsotope()
        and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        and ((atom.GetAtomicNum() in _INFIX and _sigma(atom) == 2) or (atom.GetAtomicNum() == 7 and atom.GetIsAromatic() and _sigma(atom) == 3 and atom.GetDegree() == 3))
    )


def _prefixes(mol, graph, ring_atoms, locant_of):
    from ._functional_prefixes import functional_names

    halogens = halogen_substituents(mol)
    seeds = [(a, n) for a in ring_atoms for n in graph[a] if n not in ring_atoms]
    named, shown, _ = functional_names(mol, graph, seeds, set(ring_atoms), halogens)
    entries = {}
    for atom, root in seeds:
        entry = named[root] if root in named else name_branch(graph, root, atom, shown, mol=mol)
        entries.setdefault(locant_of[atom], []).append(entry)
    return group_substituents(entries)


def name_ylium_ring(mol) -> str:
    from ._multiplicative import _ring_systems

    centre = next(a for a in mol.GetAtoms() if a.GetFormalCharge())
    graph = adjacency(mol)
    system = next(s for s in _ring_systems(mol) if centre.GetIdx() in s)
    best = None
    for smarts, locants, centre_index, core in _skeletons(centre.GetAtomicNum()):
        query = Chem.MolFromSmarts(smarts)
        if query is None or query.GetNumAtoms() != len(system):
            continue
        for match in mol.GetSubstructMatches(query, uniquify=False, useChirality=False):
            if set(match) != set(system) or match[centre_index] != centre.GetIdx():
                continue
            locant_of = {atom: locants[i] for i, atom in enumerate(match)}
            grouped = _prefixes(mol, graph, system, locant_of)
            locant_set, _, citation = substituent_locant_set_and_citation(grouped)
            key = (locant_set, citation)
            if best is None or key < best[0]:
                best = (key, grouped, core)
    if best is None:
        raise UnsupportedStructure("this cationic ring is not one of the ylium rings with a retained or lambda name")
    _, grouped, core = best
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return f"{prefix}-{core}" if prefix and core[0].isdigit() else prefix + core
