"""Monospiro ring systems of two ring components, at least one of them polycyclic (P-24.3, P-24.5, P-24.8): the
component names are cited in alphanumerical order inside 'spiro[...]' (or once inside 'spirobi[...]'), the second
component takes primed locants, the spiro atom gets the lowest locants, and indicated hydrogen and 'hydro' prefixes are
those of the complete skeleton made mancude after the union. A spiro atom with a nonstandard bonding number carries its
lambda number, and a cationic one ends in 'ylium'.
"""

import re
from itertools import combinations, product

import networkx as nx
from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, multiplied_word, ring_cycle
from ._fused_numbering import HETERO_RANK
from ._fusion_name import Context, fusion_name, system_numbering_options
from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes

_PRIME = "′"
_STANDARD = {5: 3, 7: 3, 8: 2, 15: 3, 16: 2, 33: 3, 34: 2, 51: 3, 52: 2, 83: 3}
_GROUP_14 = {6, 14, 32, 50, 82}
_HALOGENS = {9, 17, 35, 53}


def _lk(text):
    match = re.match(r"(\d+)([a-z]*)(′*)", text)
    return int(match.group(1)), len(match.group(3)), match.group(2)


def _components(mol):
    info = mol.GetRingInfo()
    rings = [set(r) for r in info.AtomRings()]
    bond_rings = [set(r) for r in info.BondRings()]
    parent = list(range(len(rings)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in combinations(range(len(rings)), 2):
        if bond_rings[i] & bond_rings[j]:
            parent[find(i)] = find(j)
    groups = {}
    for i in range(len(rings)):
        groups.setdefault(find(i), []).append(i)
    return [
        {"atoms": set().union(*(rings[i] for i in g)), "bonds": set().union(*(bond_rings[i] for i in g)), "rings": len(g)}
        for g in groups.values()
    ]


def _spiro_split(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    comps = _components(mol)
    if len(comps) != 2 or max(c["rings"] for c in comps) < 2:
        return None
    shared = comps[0]["atoms"] & comps[1]["atoms"]
    if len(shared) != 1 or any(a.IsInRing() and a.GetIdx() not in comps[0]["atoms"] | comps[1]["atoms"] for a in mol.GetAtoms()):
        return None
    return comps, next(iter(shared))


def has_spiro_union_shape(mol) -> bool:
    return _spiro_split(mol) is not None


def _skeleton(mol, comp):
    atoms = sorted(comp["atoms"])
    index = {a: i for i, a in enumerate(atoms)}
    sub = Chem.RWMol()
    for a in atoms:
        sub.AddAtom(Chem.Atom(mol.GetAtomWithIdx(a).GetAtomicNum()))
    for b in comp["bonds"]:
        bond = mol.GetBondWithIdx(b)
        sub.AddBond(index[bond.GetBeginAtomIdx()], index[bond.GetEndAtomIdx()], Chem.BondType.SINGLE)
    result = sub.GetMol()
    Chem.SanitizeMol(result)
    Chem.GetSymmSSSR(result)
    return result, atoms


def _bracket_locants(name):
    match = re.match(r"(\d+(?:,\d+)*)-(.+)", name)
    return f"[{match.group(1)}]{match.group(2)}" if match else name


def _monocycle(sub, atoms, mol):
    double = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in sub.GetBonds()
        if mol.GetBondBetweenAtoms(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()]).GetBondTypeAsDouble() == 2.0
    }
    if any(mol.GetBondBetweenAtoms(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()]).GetBondTypeAsDouble() not in (1.0, 2.0) for b in sub.GetBonds()):
        raise UnsupportedStructure("this bond order is not supported in a spiro component")
    hetero = any(a.GetAtomicNum() != 6 for a in sub.GetAtoms())
    if hetero:
        if double or not has_hetero_monocyclic_name(sub):
            raise UnsupportedStructure("this monocyclic spiro component has no supported name")
        name = name_hetero_monocyclic(sub)
    else:
        name = "cyclo" + alkane_name(sub.GetNumAtoms())
    cycle = ring_cycle(adjacency(sub), list(range(sub.GetNumAtoms())))
    size = len(cycle)
    options = []
    for start in range(size):
        for step in (1, -1):
            order = [cycle[(start + step * k) % size] for k in range(size)]
            locants = {a: k + 1 for k, a in enumerate(order)}
            symbols = [(HETERO_RANK.get(sub.GetAtomWithIdx(a).GetSymbol(), 99), locants[a]) for a in order if sub.GetAtomWithIdx(a).GetAtomicNum() != 6]
            ene = tuple(sorted(min(locants[a] for a in pair) for pair in double))
            key = (sorted(l for _, l in symbols), [l for _, l in sorted(symbols)])
            options.append((key, {a: str(l) for a, l in locants.items()}, ene))
    best = min(k for k, _, _ in options)
    return name, [(n, e) for k, n, e in options if k == best]


def _ene_name(name, ene):
    if not ene:
        return name
    stem = name[:-3]
    if len(ene) == 1:
        return f"{stem}-{ene[0]}-ene"
    return f"{stem}a-{','.join(map(str, ene))}-{multiplied_word(len(ene), 'ene')}"


def _component(mol, comp):
    sub, atoms = _skeleton(mol, comp)
    if comp["rings"] == 1:
        name, numberings = _monocycle(sub, atoms, mol)
    else:
        name, root = fusion_name(sub)
        numberings = [(n, ()) for n in system_numbering_options(Context(sub), name, root)]
        name = _bracket_locants(name)
    return name, [({atoms[i]: loc for i, loc in n.items()}, ene) for n, ene in numberings]


def _capable(mol, comp, spiro):
    result = set()
    for a in comp["atoms"]:
        if a == spiro:
            continue
        atom = mol.GetAtomWithIdx(a)
        ring_degree = sum(1 for n in atom.GetNeighbors() if mol.GetBondBetweenAtoms(a, n.GetIdx()).IsInRing())
        z = atom.GetAtomicNum()
        if z in _GROUP_14:
            capacity = 1 if ring_degree <= 3 else 0
        elif z in _STANDARD:
            if atom.GetTotalValence() - atom.GetFormalCharge() > _STANDARD[z]:
                raise UnsupportedStructure("a nonstandard bonding number away from the spiro atom is not supported yet")
            capacity = 1 if _STANDARD[z] - ring_degree == 1 else 0
        else:
            raise UnsupportedStructure("this ring atom is not supported in a spiro component yet")
        if capacity:
            result.add(a)
    return result


def _hydrogen_choices(mol, capable, polycyclic_bonds, monocyclic_bonds):
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    doubled = set()
    for bond in kekule.GetBonds():
        if bond.GetBondTypeAsDouble() == 2.0:
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bond.GetIdx() in monocyclic_bonds:
                continue
            if a not in capable or b not in capable or bond.GetIdx() not in polycyclic_bonds:
                raise UnsupportedStructure("a double bond outside the polycyclic spiro components is not supported yet")
            doubled |= {a, b}
        elif bond.GetBondTypeAsDouble() not in (1.0,):
            raise UnsupportedStructure("this bond order is not supported in a spiro system")
    saturated = capable - doubled
    graph = nx.Graph()
    graph.add_nodes_from(capable)
    for bond_idx in polycyclic_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in capable and b in capable:
            graph.add_edge(a, b)

    def unmatched(removed):
        sub = graph.subgraph(capable - set(removed))
        return len(sub) - 2 * len(nx.max_weight_matching(sub, maxcardinality=True))

    required = unmatched(())
    choices = [set(i) for i in combinations(sorted(saturated), required) if unmatched(i) == 0]
    if not choices:
        raise UnsupportedStructure("the indicated hydrogen of this spiro system cannot be placed")
    return saturated, choices


def _name_key(name):
    stripped = re.sub(r"-[\d,]+-", "", re.sub(r"^\[?\d+(?:,\d+)*\]?-?", "", name))
    lead = re.match(r"\[?(\d+(?:,\d+)*)\]?", name)
    locants = tuple(sorted(int(x) for x in lead.group(1).split(","))) if lead and re.match(r"\[?\d", name) else ()
    return alpha_sort_key(stripped), locants, name


def _spiro_atom(mol, spiro):
    atom = mol.GetAtomWithIdx(spiro)
    z = atom.GetAtomicNum()
    charge = atom.GetFormalCharge()
    if z == 6:
        if charge or atom.GetTotalNumHs() or atom.GetDegree() != 4:
            raise UnsupportedStructure("this spiro carbon is not supported")
        return None, False
    if z not in _STANDARD or charge not in (0, 1):
        raise UnsupportedStructure("this spiro atom is not supported")
    bonding = atom.GetTotalValence() + charge
    if bonding <= _STANDARD[z] or (bonding - _STANDARD[z]) % 2:
        raise UnsupportedStructure("this spiro heteroatom has no lambda bonding number")
    return bonding, bool(charge)


def name_spiro_union(mol) -> str:
    from ._ylium_ring import _prefixes

    split = _spiro_split(mol)
    comps, spiro = split
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopes, radicals and stereodescriptors of a spiro union are not supported yet")
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        raise UnsupportedStructure("double-bond stereo of a spiro union is not supported yet")
    lam, cationic = _spiro_atom(mol, spiro)
    if any(a.GetFormalCharge() and a.GetIdx() != spiro for a in mol.GetAtoms()):
        raise UnsupportedStructure("a charge away from the spiro atom is not supported yet")
    ring_atoms = comps[0]["atoms"] | comps[1]["atoms"]
    graph = adjacency(mol)
    if not cationic:
        for a in mol.GetAtoms():
            if a.GetIdx() not in ring_atoms and a.GetAtomicNum() != 6 and a.GetAtomicNum() not in _HALOGENS:
                raise UnsupportedStructure("a spiro union with a principal characteristic group is not supported yet")

    named = [_component(mol, c) for c in comps]
    capable = set()
    polycyclic_bonds = set()
    monocyclic_bonds = set()
    for comp in comps:
        if comp["rings"] > 1:
            capable |= _capable(mol, comp, spiro)
            polycyclic_bonds |= comp["bonds"]
        else:
            monocyclic_bonds |= comp["bonds"]
    if lam and not cationic:
        if lam - 4 > 1:
            raise UnsupportedStructure("a spiro atom with several multiple bonds is not supported yet")
        if lam - 4 == 1:
            capable.add(spiro)
    saturated, choices = _hydrogen_choices(mol, capable, polycyclic_bonds, monocyclic_bonds)

    order = sorted(range(2), key=lambda i: _name_key(named[i][0]))
    identical = named[0][0] == named[1][0]
    orientations = [order] if not identical else [[0, 1], [1, 0]]
    best = None
    for first, second in orientations:
        if named[second][1][0][1]:
            raise UnsupportedStructure("the locants of an unsaturated second component are not supported yet")
        for (n1, ene1), (n2, _) in product(named[first][1], named[second][1]):
            locant_of = {a: l for a, l in n1.items()}
            for a, l in n2.items():
                if a != spiro:
                    locant_of[a] = l + _PRIME
            spiro_locants = (n1[spiro], n2[spiro] + _PRIME)
            grouped = _prefixes(mol, graph, ring_atoms, locant_of)
            locant_set = sorted((_lk(l) for info in grouped.values() for l in info["locants"]))
            citation = tuple(
                _lk(l) for name in sorted(grouped, key=alpha_sort_key) for l in sorted(grouped[name]["locants"], key=_lk)
            )
            choice = min(
                (
                    (sorted(_lk(locant_of[a]) for a in picked), sorted(_lk(locant_of[a]) for a in saturated - picked), picked)
                    for picked in choices
                ),
                key=lambda c: (c[0], c[1]),
            )
            unsaturation = sorted(choice[1] + [_lk(str(x)) for x in ene1])
            key = (sorted(_lk(l) for l in spiro_locants), choice[0], unsaturation, locant_set, citation)
            if best is None or key < best[0]:
                best = (key, (first, second), spiro_locants, locant_of, choice[2], grouped, ene1)
    _, (first, second), spiro_locants, locant_of, picked, grouped, ene = best
    first_name = _ene_name(named[first][0], ene)

    indicated = ",".join(f"{locant_of[a]}H" for a in sorted(picked, key=lambda a: _lk(locant_of[a])))
    hydro_atoms = sorted(saturated - picked, key=lambda a: _lk(locant_of[a]))
    hydro = f"{','.join(locant_of[a] for a in hydro_atoms)}-{multiplied_word(len(hydro_atoms), 'hydro')}" if hydro_atoms else ""
    prefixes = format_substituent_prefixes(grouped) if grouped else ""

    lam_mark = f"λ{lam}" if lam else ""
    if identical:
        front = f"{spiro_locants[0]}{lam_mark},{spiro_locants[1]}-spirobi[{first_name}]"
    else:
        front = f"{spiro_locants[0]}{lam_mark}-" if lam else ""
        front += f"spiro[{first_name}-{spiro_locants[0]},{spiro_locants[1]}-{named[second][0]}]"
    if cationic:
        last = named[second][0]
        if last.endswith("e"):
            front = front[: -len(last) - 1] + last[:-1] + "]"
        front += f"-{spiro_locants[0]}-ylium"
    core = front

    out = prefixes
    for token in (hydro, indicated, core):
        if not token:
            continue
        joined = out and (token[0].isdigit() or token[0] in "[(" or out.endswith("H"))
        out += ("-" if joined else "") + token
    return out
