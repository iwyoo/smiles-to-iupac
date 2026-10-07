"""Ring heteroatom cations named by 'ylium' on the lambda hydride (P-73.3): the retained names of Table 7.5 and the
templates for benzopyran, furan and quinolizine come first; any other ortho- or peri-fused ring system is named from
its neutral lambda hydride (cationic centres regain a hydride) with indicated hydrogen, lambda and 'ylium' locants,
and the lowest locants go to the centres, then to the prefixes.
"""

from rdkit import Chem

from ._fused_numbering import _locant_key
from ._fusion_name import fused_parent_data, marked_name

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


_STANDARD_BONDING = {7: 3, 15: 3, 33: 3, 51: 3, 83: 3, 8: 2, 16: 2, 34: 2, 52: 2, 9: 1, 17: 1, 35: 1, 53: 1}
_MULTIPLIER = {1: "", 2: "di", 3: "tri", 4: "tetra"}


def _is_ylium_centre(atom) -> bool:
    standard = _STANDARD_BONDING.get(atom.GetAtomicNum())
    if standard is None or atom.GetFormalCharge() != 1 or not atom.IsInRing() or atom.GetIsotope():
        return False
    if atom.GetNumRadicalElectrons() or atom.GetTotalNumHs() or any(not b.IsInRing() for b in atom.GetBonds()):
        return False
    valence = atom.GetTotalValence()
    hydride = valence + 1 - standard
    quaternary = valence == atom.GetDegree() and atom.GetDegree() >= 4
    return hydride > 0 and hydride % 2 == 0 and not quaternary


def has_ylium_ring_shape(mol) -> bool:
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    return (
        bool(charged)
        and len(Chem.GetMolFrags(mol)) == 1
        and all(_is_ylium_centre(a) for a in charged)
        and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
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


def _hydride_parent(mol, system, centres):
    """The ring system alone as the neutral lambda hydride: each cationic centre regains its hydride and every
    substituent position is shown as hydrogen (P-73.3.1)."""
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    order = sorted(system)
    index = {a: i for i, a in enumerate(order)}
    parent = Chem.RWMol()
    for a in order:
        source = kekule.GetAtomWithIdx(a)
        atom = Chem.Atom(source.GetAtomicNum())
        outside = sum(1 for n in source.GetNeighbors() if n.GetIdx() not in system)
        atom.SetNumExplicitHs(source.GetTotalNumHs() + outside + (a in centres))
        atom.SetNoImplicit(True)
        parent.AddAtom(atom)
    for bond in kekule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in index and b in index and bond.IsInRing():
            parent.AddBond(index[a], index[b], bond.GetBondType())
    result = parent.GetMol()
    result.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(result)
    return result, order


def _name_fused_ylium(mol, graph, system, centres):
    parent, order = _hydride_parent(mol, system, centres)
    name, options, indicated, lam, delta = fused_parent_data(parent)
    index = {a: i for i, a in enumerate(order)}
    best = None
    for numbering in options:
        locant_of = {a: numbering[index[a]] for a in order}
        grouped = _prefixes(mol, graph, system, locant_of)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (
            sorted(_locant_key(numbering[i]) for i in indicated),
            sorted((-lam[i], _locant_key(numbering[i])) for i in lam),
            locant_set,
            citation,
        )
        if best is None or key < best[0]:
            best = (key, numbering, grouped)
    _, numbering, grouped = best
    base = marked_name(name, numbering, indicated, lam, delta)
    locants = sorted((numbering[index[a]] for a in centres), key=_locant_key)
    multiplier = _MULTIPLIER[len(locants)]
    stem = base[:-1] if base.endswith("e") and not multiplier else base
    return grouped, f"{stem}-{','.join(locants)}-{multiplier}ylium"


def name_ylium_ring(mol) -> str:
    from ._multiplicative import _ring_systems

    centres = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge()]
    graph = adjacency(mol)
    system = next(s for s in _ring_systems(mol) if centres[0] in s)
    if not all(c in system for c in centres):
        raise UnsupportedStructure("the cationic centres lie in different ring systems")
    skeletons = _skeletons(mol.GetAtomWithIdx(centres[0]).GetAtomicNum()) if len(centres) == 1 else []
    best = None
    for smarts, locants, centre_index, core in skeletons:
        query = Chem.MolFromSmarts(smarts)
        if query is None or query.GetNumAtoms() != len(system):
            continue
        for match in mol.GetSubstructMatches(query, uniquify=False, useChirality=False):
            if set(match) != set(system) or match[centre_index] != centres[0]:
                continue
            locant_of = {atom: locants[i] for i, atom in enumerate(match)}
            grouped = _prefixes(mol, graph, system, locant_of)
            locant_set, _, citation = substituent_locant_set_and_citation(grouped)
            key = (locant_set, citation)
            if best is None or key < best[0]:
                best = (key, grouped, core)
    if best is not None:
        _, grouped, core = best
    else:
        try:
            grouped, core = _name_fused_ylium(mol, graph, system, set(centres))
        except UnsupportedStructure as error:
            raise UnsupportedStructure("this cationic ring is not one of the ylium rings with a retained or lambda name") from error
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return f"{prefix}-{core}" if prefix and core[0].isdigit() else prefix + core
