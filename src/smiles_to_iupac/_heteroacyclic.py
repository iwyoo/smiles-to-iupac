"""Skeletal-replacement ('a') parents for acyclic chains with principal groups or
substituents (P-15.4.3, P-51.4.1): four or more heterounits (O, S, Se, Te, NH)
in an unbranched chain that ends in carbon, none of them part of the principal
group, give '3,6,9,12-tetraoxatetradecanedioic acid'.
"""

import itertools

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import (
    adjacency,
    group_substituents,
    halogen_substituents,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    substituent_locant_set_and_citation,
)
from ._hetero_prefixes import EXTENDED_PREFIXES, is_functional_carbon
from ._acid_groups import acid_group_at
from ._common import UnsupportedStructure, nonstandard_bonding
from ._polyfunctional import (
    CHAIN_SUFFIX_WORD,
    _SENIORITY,
    _TERMINAL,
    _group_of,
    _is_ester_like,
    _paths,
    _ring_occurrences,
)
from ._numerals import multiplying_prefix
from ._substituents import BRANCH_STEREO, format_substituent_prefixes, name_branch

_A_WORD = {
    8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma",
    14: "sila", 32: "germa", 50: "stanna", 82: "plumba", 5: "bora", 13: "alumina", 31: "galla", 49: "inda", 81: "thalla",
}
_A_ORDER = [8, 16, 34, 52, 7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81]
_STANDARD_VALENCE = {7: 3, 15: 3, 33: 3, 51: 3, 83: 3, 14: 4, 32: 4, 50: 4, 82: 4, 5: 3, 13: 3, 31: 3, 49: 3, 81: 3}
_CHAIN_ENDS = {6, *(z for z in _STANDARD_VALENCE if z != 7)}
_MINIMUM_UNITS = 4
_ORDER = [*_SENIORITY[: _SENIORITY.index("amide")], "halide", *_SENIORITY[_SENIORITY.index("amide") :]]
_TERMINAL_CLASSES = _TERMINAL | {"halide"}
_HALIDE_WORD = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_SUFFIX_WORD = {
    **{cls: word for cls, word in CHAIN_SUFFIX_WORD.items() if cls != "ide"},
    "amidine": "imidamide",
    "acid": "oic acid",
    "amide": "amide",
    "nitrile": "nitrile",
    "aldehyde": "al",
    "ketone": "one",
    "alcohol": "ol",
    "thiol": "thiol",
    "amine": "amine",
    "aminium": "aminium",
}


def _phosphorus_unit(mol, atom):
    """A neutral P(=O) with four neighbours: a chain unit 'λ5-phospha' whose other two bonds are substituents."""
    return (
        atom.GetAtomicNum() == 15
        and not atom.IsInRing()
        and not atom.GetFormalCharge()
        and not atom.GetIsotope()
        and atom.GetDegree() == 4
        and any(
            n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        )
    )


def _chalcogen_unit_valence(mol, atom):
    """The bonding number 4 or 6 of an S, Se or Te chain unit that carries one or two terminal doubly bonded chalcogen
    atoms and two chain neighbours, else None (P-14.1.3, P-21.2.3)."""
    if atom.GetAtomicNum() not in (16, 34, 52) or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope():
        return None
    terminal = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() in (8, 16, 34, 52)
        and n.GetDegree() == 1
        and not n.GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if not terminal or len(terminal) != atom.GetDegree() - 2 or len(terminal) > 2 or atom.GetTotalNumHs():
        return None
    return 2 + 2 * len(terminal)


def _chain_acyl_carbon(mol, atom):
    """A carbon with one terminal doubly bonded O, S, Se, Te or NH and two chain neighbours, at least one of them a
    heteroatom: an 'oxo', 'sulfanylidene' or 'imino' substituent of the chain (P-51.4.1.1)."""
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
        return False
    terminal = [
        n
        for n in atom.GetNeighbors()
        if n.GetDegree() == 1
        and not n.GetFormalCharge()
        and (n.GetAtomicNum() in (8, 16, 34, 52) or (n.GetAtomicNum() == 7 and n.GetTotalNumHs() == 1))
        and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in atom.GetNeighbors() if len(terminal) == 1 and n.GetIdx() != terminal[0].GetIdx()]
    return (
        len(rest) == 2
        and all(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0 for n in rest)
        and all(n.GetAtomicNum() in (6, 7, 8, 16) and not n.IsInRing() for n in rest)
        and any(n.GetAtomicNum() in (7, 8, 16) and n.GetDegree() >= 2 for n in rest)
    )


def _aminium_groups(mol):
    """{carbon: {nitrogen}} for an acyclic quaternary N+ bonded to four carbons, any one of which may be the chain."""
    found = {}
    for atom in mol.GetAtoms():
        if (
            atom.GetAtomicNum() == 7
            and atom.GetFormalCharge() == 1
            and not atom.IsInRing()
            and atom.GetDegree() == 4
            and atom.GetTotalNumHs() == 0
            and all(n.GetAtomicNum() == 6 and not n.IsInRing() for n in atom.GetNeighbors())
        ):
            for n in atom.GetNeighbors():
                found.setdefault(n.GetIdx(), set()).add(atom.GetIdx())
    return found


def _carbonless_groups(mol):
    """{carbon: (class, owned atoms)} for an acid, amide, amidine or acid halide group whose carbon has no carbon neighbour:
    the carbon is a chain member bonded to a heteroatom (P-51.4.1.3)."""
    found = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetFormalCharge():
            continue
        neighbors = atom.GetNeighbors()
        if len(neighbors) == 2 and any(_chain_heteroatom(n) for n in neighbors):
            chain = [n for n in neighbors if _chain_heteroatom(n)]
            other = [n for n in neighbors if n not in chain]
            if len(chain) == 1 and len(other) == 1 and other[0].GetDegree() == 1:
                order = mol.GetBondBetweenAtoms(atom.GetIdx(), other[0].GetIdx()).GetBondTypeAsDouble()
                if order == 3.0 and other[0].GetAtomicNum() == 7 and not atom.GetTotalNumHs():
                    found[atom.GetIdx()] = ("nitrile", {other[0].GetIdx()})
                elif order == 2.0 and other[0].GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
                    found[atom.GetIdx()] = ("aldehyde", {other[0].GetIdx()})
            continue
        if atom.GetTotalNumHs() or len(neighbors) != 3 or any(n.GetAtomicNum() == 6 for n in neighbors):
            continue
        chain = [n for n in neighbors if _chain_heteroatom(n) and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(chain) == 2 and sum(n.GetAtomicNum() == 7 for n in chain) == 1 and all(
            n.GetAtomicNum() in _STANDARD_VALENCE for n in chain
        ):
            chain = [n for n in chain if n.GetAtomicNum() != 7]
        if len(chain) != 1:
            continue
        variant = acid_group_at(mol, atom.GetIdx(), hetero_attach=True)
        if variant is not None:
            found[atom.GetIdx()] = ("acid" if variant.spec.plain else variant.spec.key, variant.owned)
            continue
        others = [n for n in neighbors if n.GetIdx() != chain[0].GetIdx()]
        double = [n for n in others if mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0 and n.GetDegree() == 1]
        single = [n for n in others if n not in double and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double) != 1 or len(single) != 1:
            continue
        top, other = double[0], single[0]
        if top.GetAtomicNum() == 8 and other.GetAtomicNum() in _HALIDE_WORD and other.GetDegree() == 1:
            found[atom.GetIdx()] = ("halide", {top.GetIdx(), other.GetIdx()})
        elif top.GetAtomicNum() == 8 and other.GetAtomicNum() == 7 and not other.IsInRing() and (
            chain[0].GetAtomicNum() != 7 or not _chain_heteroatom(other)
        ):
            found[atom.GetIdx()] = ("amide", {top.GetIdx(), other.GetIdx()})
        elif top.GetAtomicNum() == 7 and top.GetTotalNumHs() == 1 and other.GetAtomicNum() == 7 and not other.IsInRing() and not _chain_heteroatom(other):
            found[atom.GetIdx()] = ("amidine", {top.GetIdx(), other.GetIdx()})
    return found


def _carbamoyl_carbons(mol):
    """{carbon: (oxygen, [nitrogens])} for every C(=O)(N)(N) whose two nitrogens may each serve as the chain atom, the other
    being the amide nitrogen: an amide expressed as the principal group is senior to the urea (P-51.4.1.3)."""
    found = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
            continue
        oxygen = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
        nitrogens = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 7 and not n.IsInRing()]
        if len(oxygen) == 1 and len(nitrogens) == 2:
            found[atom.GetIdx()] = (oxygen[0].GetIdx(), nitrogens)
    return found


def has_carbonless_acyl(mol):
    """A carbonyl-type carbon bonded to nitrogen and to no carbon (a carbamoyl or amidine centre): such a carbon only
    belongs to a skeletal replacement chain, so that route must run before the acid-derivative namers."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
            continue
        if any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors()) or not any(n.GetAtomicNum() == 7 for n in atom.GetNeighbors()):
            continue
        if any(
            n.GetDegree() == 1
            and n.GetAtomicNum() in (7, 8, 16)
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        ):
            return True
    return False


def _imino_chain_nitrogen(atom):
    """The =N- of a chain C=N-X: two heavy neighbours, a double bond to carbon and no hydrogen."""
    if atom.GetAtomicNum() != 7 or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetTotalNumHs():
        return False
    mol = atom.GetOwningMol()
    orders = sorted(
        (mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble(), n.GetAtomicNum()) for n in atom.GetNeighbors()
    )
    return len(orders) == 2 and orders[1] == (2.0, 6) and orders[0][0] == 1.0


def _chain_heteroatom(atom, phosphorus=False):
    z = atom.GetAtomicNum()
    if phosphorus and _phosphorus_unit(atom.GetOwningMol(), atom):
        return True
    if _imino_chain_nitrogen(atom):
        return True
    if z in _STANDARD_VALENCE:
        return (
            not atom.IsInRing()
            and not atom.GetFormalCharge()
            and not atom.GetIsotope()
            and not atom.GetIsAromatic()
            and atom.GetDegree() <= _STANDARD_VALENCE[z]
            and (z != 7 or atom.GetDegree() >= 2)
            and atom.GetTotalNumHs() == _STANDARD_VALENCE[z] - atom.GetDegree()
        )
    if _chalcogen_unit_valence(atom.GetOwningMol(), atom):
        return True
    if atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetDegree() != 2:
        return False
    if z in (16, 34, 52) and nonstandard_bonding(atom) and all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds()):
        return True
    if z in (8, 16, 34, 52):
        return atom.GetTotalNumHs() == 0
    return z == 7 and atom.GetTotalNumHs() == 1 and not atom.GetIsAromatic()


def _heterounit_count(mol, path, hetero_positions):
    """The heterounits of a chain (P-51.4.1.1): two identical adjacent heteroatoms (disulfanediyl, disilane-1,2-diyl) or
    two atoms of a senior element around one of a junior element (disiloxane-1,3-diyl, but not oxysilanediyloxy) are one
    unit; every other heteroatom is a unit of its own."""
    elements = {i: mol.GetAtomWithIdx(path[i]).GetAtomicNum() for i in hetero_positions}
    taken = set()
    units = 0
    for i in hetero_positions:
        if i in taken:
            continue
        taken.add(i)
        units += 1
        if (
            elements.get(i + 1) == elements[i]
            and i + 1 not in taken
            and not _chalcogen_unit_valence(mol, mol.GetAtomWithIdx(path[i]))
            and not _chalcogen_unit_valence(mol, mol.GetAtomWithIdx(path[i + 1]))
        ):
            taken.add(i + 1)
        elif (
            i + 2 in elements
            and elements[i + 2] == elements[i]
            and i + 1 in elements
            and elements[i + 1] != elements[i]
            and _A_ORDER.index(elements[i]) > _A_ORDER.index(elements[i + 1])
            and i + 1 not in taken
            and i + 2 not in taken
        ):
            taken.update((i + 1, i + 2))
    return units


def _is_parent_hydride_run(elements):
    """A run of heteroatoms that is a parent hydride by itself: one element (trisulfane) or two alternating (P-21.2.2,
    P-21.2.3.1); any other run needs skeletal replacement in a carbon chain."""
    return len(set(elements)) == 1 or (
        len(set(elements)) == 2 and all(a != b for a, b in zip(elements, elements[1:]))
    )


def _acceptable_path(mol, path):
    """Whether the path is a chain a skeletal replacement name may use: it ends in carbon or a metalloid, carries four or
    more heterounits, and is not one block of heteroatoms that is a parent hydride by itself."""
    if mol.GetAtomWithIdx(path[0]).GetAtomicNum() not in _CHAIN_ENDS or mol.GetAtomWithIdx(path[-1]).GetAtomicNum() not in _CHAIN_ENDS:
        return False
    hetero_positions = [i for i, a in enumerate(path) if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if _heterounit_count(mol, path, hetero_positions) < _MINIMUM_UNITS:
        return False
    if len(hetero_positions) == len(path) or (
        hetero_positions == list(range(hetero_positions[0], hetero_positions[-1] + 1))
        and len(hetero_positions) >= 3
        and _is_parent_hydride_run([mol.GetAtomWithIdx(path[i]).GetAtomicNum() for i in hetero_positions])
        and not any(_chalcogen_unit_valence(mol, mol.GetAtomWithIdx(path[i])) for i in hetero_positions)
    ):
        return False
    return not _long_identical_run(mol, path, hetero_positions)


def _long_identical_run(mol, path, hetero_positions):
    """Three or more identical adjacent heteroatoms form a parent hydride of their own (trisilane, trisulfane)."""
    return any(
        b - a == 1
        and c - b == 1
        and mol.GetAtomWithIdx(path[a]).GetAtomicNum() == mol.GetAtomWithIdx(path[b]).GetAtomicNum() == mol.GetAtomWithIdx(path[c]).GetAtomicNum()
        and mol.GetAtomWithIdx(path[a]).GetAtomicNum() not in (15, 16, 34, 52)
        for a, b, c in zip(hetero_positions, hetero_positions[1:], hetero_positions[2:])
    )


def name_heteroacyclic(mol):
    """The skeletal-replacement name, or None when `mol` does not qualify."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    aminium = {} if any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()) else _aminium_groups(mol)
    hetero_atoms = [a for a in mol.GetAtoms() if _chain_heteroatom(a, bool(aminium))]
    if len(hetero_atoms) < _MINIMUM_UNITS:
        return None
    from ._polyfunctional import boranyl_amine_parent_applies

    if boranyl_amine_parent_applies(mol):
        return None
    graph = adjacency(mol)
    if any(a.GetNumRadicalElectrons() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    groups = {}
    carbonless = _carbonless_groups(mol)
    for atom in mol.GetAtoms():
        found = None if atom.GetIdx() in carbonless else _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
    for carbon, (cls, atoms) in carbonless.items():
        groups.setdefault(cls, {})[carbon] = atoms
    carbamoyl = _carbamoyl_carbons(mol)
    if aminium:
        groups["aminium"] = aminium
    ring_groups = _ring_occurrences(mol)
    classes = set(groups) | {c for c, _, _ in ring_groups} | ({"amide"} if carbamoyl else set())
    principal = "aminium" if aminium else next((c for c in _ORDER if c in classes), None)
    if principal is not None and any(c == principal for c, _, _ in ring_groups):
        return None
    if not aminium and any(
        _is_ester_like(mol, a.GetIdx()) and not _chain_acyl_carbon(mol, a) and a.GetIdx() not in carbonless
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
    ):
        return None

    senior_to_amine = principal is not None and principal not in ("amine", "imine")
    principal_atoms = groups.get(principal, {}) if principal else {}
    owned = set().union(*principal_atoms.values()) if principal_atoms else set()
    eligible = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if (
            a.GetAtomicNum() == 6
            and not a.IsInRing()
            and (
                a.GetIdx() in principal_atoms
                or not is_functional_carbon(mol, a.GetIdx())
                or _chain_acyl_carbon(mol, a)
                or (aminium and _is_ester_like(mol, a.GetIdx()))
            )
        )
        or (
            _chain_heteroatom(a, bool(aminium))
            and a.GetIdx() not in owned
            and (senior_to_amine or not (a.GetAtomicNum() == 7 and a.GetDegree() == 3))
        )
    }
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    atom_codes = {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.HasProp("_CIPCode")}
    bond_codes = {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode") for b in probe.GetBonds() if b.HasProp("_CIPCode")
    }
    best = None
    for path in _chain_paths(graph, eligible, carbamoyl if principal == "amide" else {}):
        if not _acceptable_path(mol, path):
            continue
        for chain in (path, path[::-1]):
            candidate = _evaluate(mol, graph, chain, principal, principal_atoms, owned, atom_codes, bond_codes, carbamoyl)
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    if best is None:
        return None
    return best[1]


def _chain_paths(graph, eligible, carbamoyl):
    """The maximal chains of `eligible`, plus those that end at a carbamoyl carbon with one of its two nitrogens cut off
    as the amide nitrogen."""
    seen = set()
    for choice in itertools.product(*([None, *v[1]] for v in carbamoyl.values())):
        adj = {a: list(ns) for a, ns in graph.items()}
        for carbon, cut in zip(carbamoyl, choice):
            if cut is not None:
                adj[carbon].remove(cut)
                adj[cut].remove(carbon)
        for path in _paths(adj, eligible):
            if tuple(path) not in seen:
                seen.add(tuple(path))
                yield path


def _chain_stereo(chain, position_of, atom_codes, bond_codes):
    entries = [(position_of[a], atom_codes[a]) for a in chain if a in atom_codes]
    for (a, b), code in bond_codes.items():
        if a in position_of and b in position_of and abs(position_of[a] - position_of[b]) == 1:
            entries.append((min(position_of[a], position_of[b]), code))
    entries.sort()
    return f"({','.join(f'{p}{c}' for p, c in entries)})-" if entries else ""


def _evaluate(
    mol, graph, chain, principal, principal_atoms, owned, atom_codes, bond_codes, carbamoyl=None, attach=None, blocked=frozenset(),
    attaches=(),
):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    if attach is not None and attach not in chain_set:
        return None
    if attaches and (len(attaches) != 2 or any(a not in chain_set for a in attaches)):
        return None
    if attaches:
        attach = attaches[0]
    if principal == "amide" and carbamoyl:
        principal_atoms = dict(principal_atoms)
        for end, neighbor in ((chain[0], chain[1]), (chain[-1], chain[-2])):
            if end in carbamoyl and end not in principal_atoms:
                amide_nitrogen = [n for n in carbamoyl[end][1] if n != neighbor]
                if len(amide_nitrogen) == 1 and amide_nitrogen[0] not in chain_set:
                    principal_atoms[end] = {carbamoyl[end][0], amide_nitrogen[0]}
        owned = owned | set().union(*principal_atoms.values()) if principal_atoms else owned
    on_chain = [a for a in principal_atoms if a in chain_set]
    if principal is not None and not on_chain:
        return None
    owned = set().union(*(principal_atoms[a] for a in on_chain)) if on_chain else set()
    if principal in _TERMINAL_CLASSES and any(position_of[a] not in (1, len(chain)) for a in on_chain):
        return None
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    ene, yne = [], []
    for a, b in zip(chain, chain[1:]):
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order == 2.0:
            ene.append(position_of[a])
        elif order == 3.0:
            yne.append(position_of[a])
    context = {"atoms": atom_codes, "bonds": bond_codes, "used": set()}
    token = BRANCH_STEREO.set(context)
    extended = EXTENDED_PREFIXES.set(True)
    entries = {}
    nitrogen_entries = {}
    try:
        for atom in chain:
            if (
                mol.GetAtomWithIdx(atom).GetAtomicNum() != 6
                and mol.GetAtomWithIdx(atom).GetAtomicNum() not in _STANDARD_VALENCE
                and not _chalcogen_unit_valence(mol, mol.GetAtomWithIdx(atom))
            ):
                continue
            for neighbor in graph[atom]:
                if neighbor in chain_set or neighbor in owned or neighbor in blocked:
                    continue
                name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                entries.setdefault(position_of[atom], []).append((name, compound))
        if principal == "aminium":
            for carbon in on_chain:
                for nitrogen in principal_atoms[carbon]:
                    for neighbor in graph[nitrogen]:
                        if neighbor != carbon:
                            name, compound = name_branch(graph, neighbor, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                            nitrogen_entries.setdefault("N", []).append((name, compound))
        elif principal in ("amide", "amidine"):
            primes = ("N", "N'", "N''", "N'''")
            for rank, carbon in enumerate(sorted(on_chain, key=position_of.get)):
                nitrogens = sorted(
                    (n for n in principal_atoms[carbon] if mol.GetAtomWithIdx(n).GetAtomicNum() == 7),
                    key=lambda n: mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble(),
                )
                for offset, nitrogen in enumerate(nitrogens):
                    label = primes[rank + offset] if principal == "amidine" and len(on_chain) == 1 else primes[rank]
                    if principal == "amidine" and len(on_chain) > 1:
                        label = primes[2 * rank + offset]
                    for neighbor in graph[nitrogen]:
                        if neighbor != carbon:
                            name, compound = name_branch(graph, neighbor, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                            nitrogen_entries.setdefault(label, []).append((name, compound))
    except UnsupportedStructure:
        return None
    finally:
        EXTENDED_PREFIXES.reset(extended)
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in atom_codes if a not in chain_set) or any(
        ("bond", b) not in context["used"] for b in bond_codes if not (b[0] in chain_set and b[1] in chain_set)
    ):
        return None
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    grouped_all = group_substituents({**entries, **nitrogen_entries})

    by_element = {}
    for atom in chain:
        z = mol.GetAtomWithIdx(atom).GetAtomicNum()
        if z != 6:
            by_element.setdefault(z, []).append(position_of[atom])
    hetero_set = lowest_locant_set(p for ps in by_element.values() for p in ps)
    hetero_order = tuple(tuple(by_element.get(z, ())) for z in _A_ORDER)
    suffix_locants = sorted(position_of[a] for a in on_chain)
    count = len(on_chain)
    length = len(chain)

    lambda_of = {position_of[a]: 5 for a in chain if _phosphorus_unit(mol, mol.GetAtomWithIdx(a))}
    lambda_of.update(
        {position_of[a]: v for a in chain if (v := _chalcogen_unit_valence(mol, mol.GetAtomWithIdx(a)))}
    )
    lambda_of.update(
        {position_of[a]: v for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and (v := nonstandard_bonding(mol.GetAtomWithIdx(a)))}
    )

    def a_unit(z):
        locants = sorted(by_element[z])
        multiplier = multiplying_prefix(len(locants)) if len(locants) > 1 else ""
        cited = ",".join(f"{p}λ{lambda_of[p]}" if p in lambda_of else str(p) for p in locants)
        return f"{cited}-{multiplier}{_A_WORD[z]}"

    a_text = "-".join(a_unit(z) for z in _A_ORDER if z in by_element)
    prefix = format_substituent_prefixes(grouped_all)
    halide_word = ""
    if attaches:
        body = name_from_substituents(
            length, ene, yne, "diyl", sorted(position_of[a] for a in attaches), force_own_locant=True
        )
    elif attach is not None:
        body = name_from_substituents(length, ene, yne, "yl", [position_of[attach]], force_own_locant=True)
    elif principal is None:
        body = name_from_substituents(length, ene, yne, "e")
    else:
        if principal == "halide":
            halogens_cited = {mol.GetAtomWithIdx(i).GetAtomicNum() for a in on_chain for i in principal_atoms[a] if mol.GetAtomWithIdx(i).GetAtomicNum() in _HALIDE_WORD}
            if len(halogens_cited) != 1:
                return None
            halide_word = " " + ("di" if count == 2 else "") + _HALIDE_WORD[halogens_cited.pop()]
            suffix_word = multiplied_word(count, "oyl")
        else:
            if principal not in _SUFFIX_WORD:
                return None
            suffix_word = multiplied_word(count, _SUFFIX_WORD[principal])
        body = name_from_substituents(length, ene, yne, suffix_word, suffix_locants) + halide_word
    name = prefix + ("-" if prefix and a_text and not prefix.endswith("-") else "") + a_text + body
    stereo = _chain_stereo(chain, position_of, atom_codes, bond_codes)
    if attach is not None:
        key = (
            -sum(len(ps) for ps in by_element.values()),
            -length,
            tuple(-len(by_element.get(z, ())) for z in _A_ORDER),
            -(len(ene) + len(yne)),
            -len(ene),
            -len(lambda_of),
            tuple(-v for v in sorted(lambda_of.values(), reverse=True)),
            hetero_set,
            hetero_order,
            position_of[attach],
            lowest_locant_set(ene + yne),
            lowest_locant_set(ene),
            tuple(sorted(lambda_of)),
            tuple(-lambda_of[p] for p in sorted(lambda_of)),
            -total_count,
            locant_set,
            citation,
            name,
        )
        return key, stereo + name
    key = (
        -count,
        -sum(len(ps) for ps in by_element.values()),
        -length,
        hetero_set,
        hetero_order,
        tuple(sorted(lambda_of)),
        tuple(-lambda_of[p] for p in sorted(lambda_of)),
        tuple(suffix_locants),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        locant_set,
        citation,
        name,
    )
    return key, stereo + name


def skeletal_linker(mol, graph, span, attaches, blocked):
    """The 'a' name of the divalent chain `span` joining two multiplied units ('3,6,9,12-tetraoxatetradecane-1,14-diyl',
    P-51.3.1), or None when no chain with four heterounits runs between the two attachment atoms."""
    atoms = [mol.GetAtomWithIdx(a) for a in span]
    if len(attaches) != 2 or any(a.IsInRing() or a.GetFormalCharge() or a.GetIsotope() for a in atoms):
        return None
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in atoms):
        return None
    if any(
        mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() != 1.0 for a in attaches for n in graph[a] if n in blocked
    ):
        return None
    eligible = {
        a.GetIdx()
        for a in atoms
        if (a.GetAtomicNum() == 6 and (not is_functional_carbon(mol, a.GetIdx()) or _chain_acyl_carbon(mol, a)))
        or _chain_heteroatom(a)
    }
    adj = {a: [n for n in graph[a] if n in span] for a in span}
    best = None
    for path in _paths(adj, eligible):
        if set(path) != set(span) or {path[0], path[-1]} != set(attaches) or not _acceptable_path(mol, path):
            continue
        for chain in (path, path[::-1]):
            candidate = _evaluate(mol, graph, chain, None, {}, set(), {}, {}, attaches=tuple(attaches), blocked=blocked)
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    return best[1] if best else None


def skeletal_substituent(mol, graph, root, coming_from):
    """The skeletal replacement group entered at carbon `root` when its principal substituent chain carries four or
    more heterounits (P-46.0, P-46.1): the chain with the most heteroatoms, then the longest, then the first of the
    criteria (a)-(m) that decides; None when no such chain exists."""
    if mol.GetAtomWithIdx(root).GetAtomicNum() != 6 or mol.GetAtomWithIdx(root).IsInRing():
        return None
    arm = {root}
    stack = [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != coming_from and n not in arm:
                arm.add(n)
                stack.append(n)
    atoms = [mol.GetAtomWithIdx(a) for a in arm]
    if (
        sum(1 for a in atoms if _chain_heteroatom(a)) < _MINIMUM_UNITS
        or any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in atoms)
        or any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in atoms)
        or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0
    ):
        return None
    eligible = {
        a.GetIdx()
        for a in atoms
        if (a.GetAtomicNum() == 6 and not a.IsInRing() and (not is_functional_carbon(mol, a.GetIdx()) or _chain_acyl_carbon(mol, a)))
        or _chain_heteroatom(a)
    }
    if root not in eligible:
        return None
    adj = {a: [n for n in graph[a] if n in arm] for a in arm}
    best = None
    for path in _paths(adj, eligible):
        if root not in path or not _acceptable_path(mol, path):
            continue
        for chain in (path, path[::-1]):
            candidate = _evaluate(mol, graph, chain, None, {}, set(), {}, {}, attach=root, blocked=frozenset({coming_from}))
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    return (best[1], True) if best else None


def _terminal_ester_groups(mol):
    """[(acyl carbon, ester oxygen, organyl root)] for every -C(=O)-O-R end group whose organyl R carries only carbon
    and halogen atoms and whose carbon continues the chain through exactly one carbon neighbour."""
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
            continue
        oxo = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        oxygens = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8 and n.GetDegree() == 2 and not n.IsInRing() and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(oxo) != 1 or len(oxygens) != 1 or len(carbons) != 1:
            continue
        roots = [n for n in oxygens[0].GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        if len(roots) != 1 or roots[0].GetAtomicNum() != 6:
            continue
        seen, stack = {roots[0].GetIdx()}, [roots[0].GetIdx()]
        while stack:
            for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                if n.GetIdx() != oxygens[0].GetIdx() and n.GetIdx() not in seen:
                    seen.add(n.GetIdx())
                    stack.append(n.GetIdx())
        if atom.GetIdx() in seen or any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (6, 9, 17, 35, 53) for i in seen):
            continue
        found.append((atom.GetIdx(), oxygens[0].GetIdx(), roots[0].GetIdx(), seen))
    return found


def name_heteroacyclic_ester(mol):
    """'dimethyl 3,8,10,15-tetraoxo-4,7,11,14-tetraoxaheptadecane-1,17-dioate' (P-51.4.1.1, P-65.6.3.2): the free acid
    of the chain named by skeletal replacement, with the organyl groups of its terminal esters cited first; None when
    the molecule is not of that shape."""
    from ._acid_derivatives import _multiplied, anion_name

    terminal = _terminal_ester_groups(mol)
    if not terminal or len(Chem.GetMolFrags(mol)) != 1 or Chem.FindMolChiralCenters(mol, includeUnassigned=False):
        return None
    editable = Chem.RWMol(mol)
    for _, oxygen, _, far in terminal:
        atom = editable.GetAtomWithIdx(oxygen)
        atom.SetNumExplicitHs(1)
        atom.SetNoImplicit(True)
    removed = set().union(*(far for _, _, _, far in terminal))
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    name = name_heteroacyclic(acid)
    if name is None or not name.endswith("oic acid"):
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    counts = {}
    for _, oxygen, root, _ in terminal:
        organyl, compound = name_branch(graph, root, oxygen, halogens, mol=mol)
        entry = counts.setdefault(organyl, [0, compound])
        entry[0] += 1
    pendants = _multiplied([(n, compound, k, []) for n, (k, compound) in counts.items()])
    return f"{pendants} {anion_name(name)}"
