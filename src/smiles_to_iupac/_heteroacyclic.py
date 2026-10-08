"""Skeletal-replacement ('a') parents for acyclic chains with principal groups or
substituents (P-15.4.3, P-51.4.1): four or more heterounits (O, S, Se, Te, NH)
in an unbranched chain that ends in carbon, none of them part of the principal
group, give '3,6,9,12-tetraoxatetradecanedioic acid'.
"""

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
from ._hetero_prefixes import is_functional_carbon
from ._polyfunctional import (
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
    14: "sila", 32: "germa", 50: "stanna", 82: "plumba", 5: "bora", 13: "aluma", 31: "galla", 49: "indiga", 81: "thalla",
}
_A_ORDER = [8, 16, 34, 52, 7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81]
_STANDARD_VALENCE = {15: 3, 33: 3, 51: 3, 83: 3, 14: 4, 32: 4, 50: 4, 82: 4, 5: 3, 13: 3, 31: 3, 49: 3, 81: 3}
_CHAIN_ENDS = {6, *_STANDARD_VALENCE}
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


def _chain_acyl_carbon(mol, atom):
    """A carbonyl carbon bonded to two chain atoms (C and N or O): an 'oxo' substituent of the chain (P-51.4.1.1)."""
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
        return False
    oxo = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in atom.GetNeighbors() if n.GetIdx() != oxo[0].GetIdx()] if len(oxo) == 1 else []
    return (
        len(rest) == 2
        and sorted(n.GetAtomicNum() for n in rest) in ([6, 7], [6, 8])
        and all(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0 for n in rest)
        and any(n.GetAtomicNum() in (7, 8) and n.GetDegree() == 2 and not n.IsInRing() for n in rest)
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


def _chain_heteroatom(atom, phosphorus=False):
    z = atom.GetAtomicNum()
    if phosphorus and _phosphorus_unit(atom.GetOwningMol(), atom):
        return True
    if z in _STANDARD_VALENCE:
        return (
            not atom.IsInRing()
            and not atom.GetFormalCharge()
            and not atom.GetIsotope()
            and not atom.GetIsAromatic()
            and atom.GetDegree() <= _STANDARD_VALENCE[z]
            and atom.GetTotalNumHs() == _STANDARD_VALENCE[z] - atom.GetDegree()
        )
    if atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetDegree() != 2:
        return False
    if z in (8, 16, 34, 52):
        return atom.GetTotalNumHs() == 0
    return z == 7 and atom.GetTotalNumHs() == 1 and not atom.GetIsAromatic()


def _heterounit_count(mol, path, hetero_positions):
    """The heteroatom units of a chain: adjacent atoms of one sulfur, selenium or tellurium kind (a disulfide) count
    once, so that 1-(methyldiselanyl)-2-(methyldisulfanyl)ethane keeps its substitutive name (P-63.3.1)."""
    units = len(hetero_positions)
    for a, b in zip(hetero_positions, hetero_positions[1:]):
        z = mol.GetAtomWithIdx(path[a]).GetAtomicNum()
        if b - a == 1 and z == mol.GetAtomWithIdx(path[b]).GetAtomicNum() and z in (16, 34, 52):
            units -= 1
    return units


def _is_parent_hydride_run(elements):
    """A run of heteroatoms that is a parent hydride by itself: one element (trisulfane) or two alternating (P-21.2.2,
    P-21.2.3.1); any other run needs skeletal replacement in a carbon chain."""
    return len(set(elements)) == 1 or (
        len(set(elements)) == 2 and all(a != b for a, b in zip(elements, elements[1:]))
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
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
    if aminium:
        groups["aminium"] = aminium
    ring_groups = _ring_occurrences(mol)
    classes = set(groups) | {c for c, _, _ in ring_groups}
    principal = "aminium" if aminium else next((c for c in _SENIORITY if c in classes), None)
    if principal is not None and any(c == principal for c, _, _ in ring_groups):
        return None
    if not aminium and any(
        _is_ester_like(mol, a.GetIdx()) and not _chain_acyl_carbon(mol, a) for a in mol.GetAtoms() if a.GetAtomicNum() == 6
    ):
        return None

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
        or (_chain_heteroatom(a, bool(aminium)) and a.GetIdx() not in owned)
    }
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    atom_codes = {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.HasProp("_CIPCode")}
    bond_codes = {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode") for b in probe.GetBonds() if b.HasProp("_CIPCode")
    }
    best = None
    for path in _paths(graph, eligible):
        if mol.GetAtomWithIdx(path[0]).GetAtomicNum() not in _CHAIN_ENDS or mol.GetAtomWithIdx(path[-1]).GetAtomicNum() not in _CHAIN_ENDS:
            continue
        hetero_positions = [i for i, a in enumerate(path) if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
        if _heterounit_count(mol, path, hetero_positions) < _MINIMUM_UNITS:
            continue
        if len(hetero_positions) == len(path) or (
            hetero_positions == list(range(hetero_positions[0], hetero_positions[-1] + 1))
            and len(hetero_positions) >= 3
            and _is_parent_hydride_run([mol.GetAtomWithIdx(path[i]).GetAtomicNum() for i in hetero_positions])
        ):
            # one block of heteroatoms is an alternating or homogeneous parent hydride (P-21.2.2, P-21.2.3.1)
            continue
        if any(
            b - a == 1
            and (
                mol.GetAtomWithIdx(path[a]).GetAtomicNum() == mol.GetAtomWithIdx(path[b]).GetAtomicNum()
                and mol.GetAtomWithIdx(path[a]).GetAtomicNum() not in (15, 16, 34, 52)
            )
            for a, b in zip(hetero_positions, hetero_positions[1:])
        ):
            continue
        for chain in (path, path[::-1]):
            candidate = _evaluate(mol, graph, chain, principal, principal_atoms, owned, atom_codes, bond_codes)
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    if best is None:
        return None
    return best[1]


def _chain_stereo(chain, position_of, atom_codes, bond_codes):
    entries = [(position_of[a], atom_codes[a]) for a in chain if a in atom_codes]
    for (a, b), code in bond_codes.items():
        if a in position_of and b in position_of and abs(position_of[a] - position_of[b]) == 1:
            entries.append((min(position_of[a], position_of[b]), code))
    entries.sort()
    return f"({','.join(f'{p}{c}' for p, c in entries)})-" if entries else ""


def _evaluate(mol, graph, chain, principal, principal_atoms, owned, atom_codes, bond_codes):
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
    context = {"atoms": atom_codes, "bonds": bond_codes, "used": set()}
    token = BRANCH_STEREO.set(context)
    entries = {}
    nitrogen_entries = {}
    try:
        for atom in chain:
            if mol.GetAtomWithIdx(atom).GetAtomicNum() != 6 and mol.GetAtomWithIdx(atom).GetAtomicNum() not in _STANDARD_VALENCE:
                continue
            for neighbor in graph[atom]:
                if neighbor in chain_set or neighbor in owned:
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
    finally:
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

    lambda5_positions = {position_of[a] for a in chain if _phosphorus_unit(mol, mol.GetAtomWithIdx(a))}

    def a_unit(z):
        locants = sorted(by_element[z])
        multiplier = multiplying_prefix(len(locants)) if len(locants) > 1 else ""
        cited = ",".join(f"{p}λ5" if p in lambda5_positions else str(p) for p in locants)
        return f"{cited}-{multiplier}{_A_WORD[z]}"

    a_text = "-".join(a_unit(z) for z in _A_ORDER if z in by_element)
    prefix = format_substituent_prefixes(grouped_all)
    if principal is None:
        body = name_from_substituents(length, ene, yne, "e")
    else:
        body = name_from_substituents(length, ene, yne, multiplied_word(count, _SUFFIX_WORD[principal]), suffix_locants)
    name = prefix + ("-" if prefix and a_text and not prefix.endswith("-") else "") + a_text + body
    stereo = _chain_stereo(chain, position_of, atom_codes, bond_codes)
    key = (
        -count,
        -sum(len(ps) for ps in by_element.values()),
        -length,
        hetero_set,
        hetero_order,
        tuple(suffix_locants),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        locant_set,
        citation,
        name,
    )
    return key, stereo + name


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
