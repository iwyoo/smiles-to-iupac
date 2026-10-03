"""Skeletal-replacement ('a') parents for acyclic chains with principal groups or
substituents (P-15.4.3, P-51.4.1): four or more heterounits (O, S, Se, Te, NH)
in an unbranched chain that ends in carbon, none of them part of the principal
group, give '3,6,9,12-tetraoxatetradecanedioic acid'.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    substituent_locant_set_and_citation,
)
from ._hetero_prefixes import is_functional_carbon
from ._polyfunctional import (
    _SENIORITY,
    _TERMINAL,
    _group_of,
    _is_ester_like,
    _paths,
    _ring_occurrences,
    _stereo_free,
)
from ._numerals import multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch

_A_WORD = {8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza"}
_A_ORDER = [8, 16, 34, 52, 7]
_MINIMUM_UNITS = 4
_SUFFIX_WORD = {
    "acid": "oic acid",
    "amide": "amide",
    "nitrile": "nitrile",
    "aldehyde": "al",
    "ketone": "one",
    "alcohol": "ol",
    "thiol": "thiol",
    "amine": "amine",
}


def _chain_heteroatom(atom):
    z = atom.GetAtomicNum()
    if atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetDegree() != 2:
        return False
    if z in (8, 16, 34, 52):
        return atom.GetTotalNumHs() == 0
    return z == 7 and atom.GetTotalNumHs() == 1 and not atom.GetIsAromatic()


def name_heteroacyclic(mol):
    """The skeletal-replacement name, or None when `mol` does not qualify."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    hetero_atoms = [a for a in mol.GetAtoms() if _chain_heteroatom(a)]
    if len(hetero_atoms) < _MINIMUM_UNITS:
        return None
    graph = adjacency(mol)
    if any(a.GetNumRadicalElectrons() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    groups = {}
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
    ring_groups = _ring_occurrences(mol)
    classes = set(groups) | {c for c, _, _ in ring_groups}
    principal = next((c for c in _SENIORITY if c in classes), None)
    if principal is not None and any(c == principal for c, _, _ in ring_groups):
        return None
    if any(_is_ester_like(mol, a.GetIdx()) for a in mol.GetAtoms() if a.GetAtomicNum() == 6):
        return None

    principal_atoms = groups.get(principal, {}) if principal else {}
    owned = set().union(*principal_atoms.values()) if principal_atoms else set()
    eligible = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if (
            a.GetAtomicNum() == 6
            and not a.IsInRing()
            and (a.GetIdx() in principal_atoms or not is_functional_carbon(mol, a.GetIdx()))
        )
        or (_chain_heteroatom(a) and a.GetIdx() not in owned)
    }
    best = None
    for path in _paths(graph, eligible):
        if mol.GetAtomWithIdx(path[0]).GetAtomicNum() != 6 or mol.GetAtomWithIdx(path[-1]).GetAtomicNum() != 6:
            continue
        hetero_positions = [i for i, a in enumerate(path) if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
        if len(hetero_positions) < _MINIMUM_UNITS:
            continue
        if any(b - a == 1 for a, b in zip(hetero_positions, hetero_positions[1:])):
            continue
        for chain in (path, path[::-1]):
            candidate = _evaluate(mol, graph, chain, principal, principal_atoms, owned)
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    if best is None:
        return None
    if not _stereo_free(mol):
        raise UnsupportedStructure("stereodescriptors with a skeletal-replacement parent are not supported yet")
    return best[1]


def _evaluate(mol, graph, chain, principal, principal_atoms, owned):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    on_chain = [a for a in principal_atoms if a in chain_set]
    if principal is not None and not on_chain:
        return None
    if principal in _TERMINAL and any(position_of[a] not in (1, len(chain)) for a in on_chain):
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
    entries = {}
    for atom in chain:
        if mol.GetAtomWithIdx(atom).GetAtomicNum() != 6:
            continue
        for neighbor in graph[atom]:
            if neighbor in chain_set or neighbor in owned:
                continue
            name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            entries.setdefault(position_of[atom], []).append((name, compound))
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)

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

    a_text = "-".join(
        f"{','.join(str(p) for p in sorted(by_element[z]))}-{(multiplying_prefix(len(by_element[z])) if len(by_element[z]) > 1 else '')}{_A_WORD[z]}"
        for z in _A_ORDER
        if z in by_element
    )
    prefix = format_substituent_prefixes(grouped)
    if principal is None:
        body = name_from_substituents(length, ene, yne, "e")
    elif principal in _TERMINAL:
        body = name_from_substituents(length, ene, yne, multiplied_word(count, _SUFFIX_WORD[principal]))
    else:
        body = name_from_substituents(length, ene, yne, multiplied_word(count, _SUFFIX_WORD[principal]), suffix_locants)
    name = prefix + ("-" if prefix and a_text and not prefix.endswith("-") else "") + a_text + body
    key = (
        -count,
        hetero_set,
        hetero_order,
        tuple(suffix_locants),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        locant_set,
        citation,
        name,
    )
    return key, name
