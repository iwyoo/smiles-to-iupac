"""Naming of saturated or mancude monocycles whose ring chalcogen atoms (S/Se/Te) carry
doubly bonded chalcogen atoms (ring sulfoxides, sulfones, sultines, sultones),
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-64.4.2: a ring -SO-/-SO2- group is expressed substitutively by the
  suffix 'one' on the heterocycle whose chalcogen atoms are designated
  λ4/λ6 (P-14.1.3), e.g. 'thiolane 1,1-dioxide' -> '1λ6-thiolane-1,1-dione'.
  P-65.6.3.5.2 gives the same construction for sultones and sultines
  ('1,2λ6-oxathiolane-2,2-dione', '1,2λ4-oxathiolane-2-thione').
- P-22.2.2.1: the ring is a Hantzsch-Widman name built from the prefixes
  'oxa'/'thia'/'selena'/'tellura'/'aza' (decreasing seniority, multiplied,
  final 'a' elided before a vowel) and the ring-size endings; the λ number is
  written directly after the locant of its atom.
- P-22.2.7.1, P-74.2.2.1.8: in a mancude ring the λ atom is a saturated atom,
  so its locant is also the indicated hydrogen ('1H-1λ4-thiophen-1-one',
  '1H-1λ6-thiophene-1,1-dione'). Only rings whose every other carbon/nitrogen
  atom takes part in a ring double bond are handled; thiophene, selenophene
  and tellurophene are the retained names (P-22.2.1).
- Numbering: lowest locants to all heteroatoms, then to O before S before
  Se before Te (P-31.1.4.2.4), then to the atoms with nonstandard bonding
  numbers, the higher number first on a tie (P-44.4.1.3.2), then to the
  principal characteristic groups ('one'/'thione'/...), then to the
  detachable prefixes (P-31.1.4.3.4).
- A doubly bonded chalcogen of a senior kind is the suffix; the others are
  prefixes ('oxo', 'sulfanylidene', ...).
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    multiplied_word,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import numerical_term
from ._substituents import format_substituent_prefixes, name_branch

_RING_HETERO = {8: "O", 16: "S", 34: "Se", 52: "Te", 15: "P"}
_LAMBDA_HETERO = {16, 34, 52, 15}
_CHALCOGENS = {8, 16, 34, 52}
_SENIORITY = ("O", "S", "Se", "Te", "N", "P")
_RING_PREFIX = {"O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha"}
_STEM_ENDING = {3: "irane", 4: "etane", 5: "olane", 6: "ane", 7: "epane", 8: "ocane", 9: "onane", 10: "ecane"}
_MANCUDE_ENDING = {3: "irene", 4: "ete", 5: "ole", 6: "ine", 7: "epine", 8: "ocine", 9: "onine", 10: "ecine"}
_RETAINED = {"S": "thiophene", "Se": "selenophene", "Te": "tellurophene"}
_SUFFIX = {"O": "one", "S": "thione", "Se": "selone", "Te": "tellone"}
_PREFIX = {"O": "oxo", "S": "sulfanylidene", "Se": "selanylidene", "Te": "tellanylidene"}


def _exocyclic_chalcogens(mol, atom_idx, ring_set):
    """The doubly bonded, neutral, terminal chalcogen neighbors of a ring atom."""
    found = []
    for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
        other = bond.GetOtherAtom(mol.GetAtomWithIdx(atom_idx))
        if other.GetIdx() in ring_set or bond.GetBondTypeAsDouble() != 2.0:
            continue
        if other.GetAtomicNum() in _CHALCOGENS and other.GetDegree() == 1 and other.GetFormalCharge() == 0:
            found.append(other.GetIdx())
    return found


def _match(mol):
    """(ring_order, ring_hetero, exo_chalcogen_by_ring_atom, ring_substituent_roots)
    for a saturated monocycle of carbon and O/S/Se/Te atoms where at least one
    S/Se/Te carries doubly bonded chalcogens; None for any other shape."""
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1 or not 3 <= len(ring_info.AtomRings()[0]) <= 10:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_set = set(ring_atoms)
    ring_double = {
        i: sum(
            1
            for b in mol.GetAtomWithIdx(i).GetBonds()
            if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtomIdx(i) in ring_set
        )
        for i in ring_atoms
    }
    mancude = any(ring_double.values())
    ring_hetero = {}
    exo_chalcogens = {}
    roots = {}
    graph = adjacency(mol)
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetIsAromatic() or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        exo = [n for n in graph[idx] if n not in ring_set]
        if atom.GetAtomicNum() == 6:
            terminal = _exocyclic_chalcogens(mol, idx, ring_set)
            if mancude and (terminal or ring_double[idx] != 1):
                return None
            if terminal and len(exo) != 1:
                return None
            if terminal:
                exo_chalcogens[idx] = terminal
            else:
                roots[idx] = exo
        elif atom.GetAtomicNum() == 7 and mancude:
            if ring_double[idx] != 1 or exo or atom.GetTotalNumHs() != 0:
                return None
            ring_hetero[idx] = "N"
        elif atom.GetAtomicNum() in _RING_HETERO:
            if ring_double[idx]:
                return None
            ring_hetero[idx] = _RING_HETERO[atom.GetAtomicNum()]
            terminal = _exocyclic_chalcogens(mol, idx, ring_set)
            if atom.GetAtomicNum() == 15:
                if mancude or len(terminal) != 1:
                    return None
                exo_chalcogens[idx] = terminal
                roots[idx] = [n for n in exo if n != terminal[0]]
                continue
            if len(terminal) != len(exo) or (terminal and atom.GetAtomicNum() not in _LAMBDA_HETERO):
                return None
            if terminal:
                exo_chalcogens[idx] = terminal
        else:
            return None
    if not any(i in exo_chalcogens for i in ring_hetero):
        return None
    for i in range(len(ring_atoms)):
        a, b = ring_atoms[i], ring_atoms[(i + 1) % len(ring_atoms)]
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond is not None and not mancude and bond.GetBondTypeAsDouble() != 1.0:
            return None
    exo_atoms = {x for xs in exo_chalcogens.values() for x in xs}
    phosphorus_branch = set()
    stack = [r for i in ring_hetero if ring_hetero[i] == "P" for r in roots[i]]
    while stack:
        i = stack.pop()
        if i not in phosphorus_branch:
            phosphorus_branch.add(i)
            stack.extend(n for n in graph[i] if n not in ring_set)
    for atom in mol.GetAtoms():
        if atom.GetIdx() in ring_set or atom.GetIdx() in exo_atoms:
            continue
        if atom.GetAtomicNum() == 8 and atom.GetIdx() in phosphorus_branch:
            continue
        if atom.GetAtomicNum() not in (6, 9, 17, 35, 53):
            return None
    if any(
        b.GetBondTypeAsDouble() != 1.0
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() not in exo_atoms
        and b.GetEndAtomIdx() not in exo_atoms
        and not (mancude and b.IsInRing())
    ):
        return None
    return ring_cycle(graph, ring_atoms), ring_hetero, exo_chalcogens, roots


def has_ring_lambda_heterone_shape(mol) -> bool:
    return _match(mol) is not None


def _stem(elements, size, mancude):
    """Hantzsch-Widman stem (no locants) for the ring heteroatom `elements` on
    a ring of `size` atoms."""
    if mancude and len(elements) == 1 and size == 5 and elements[0] in _RETAINED:
        return _RETAINED[elements[0]]
    counts = {e: elements.count(e) for e in _SENIORITY if e in elements}
    prefix = ""
    for e, n in counts.items():
        term = (numerical_term(n) if n > 1 else "") + _RING_PREFIX[e]
        prefix = (prefix[:-1] if prefix and term[0] in "aeiou" else prefix) + term
    ending = (_MANCUDE_ENDING if mancude else _STEM_ENDING)[size]
    if not mancude and size == 6 and "P" in elements:
        ending = "inane"
    return (prefix[:-1] if ending[0] in "aeiou" else prefix) + ending


def _lambda(element, exo):
    return (3 if element == "P" else 2) + 2 * len(exo)


def _numbering_key(position_of, ring_hetero, exo_chalcogens, suffix_kind):
    hetero_locants = sorted(position_of[a] for a in ring_hetero)
    seniority = tuple(
        position_of[a]
        for element in _SENIORITY
        for a in sorted((x for x, e in ring_hetero.items() if e == element), key=position_of.get)
    )
    lam = {position_of[a]: _lambda(ring_hetero[a], exo_chalcogens[a]) for a in ring_hetero if a in exo_chalcogens}
    lam_locants = sorted(lam)
    lam_numbers = tuple(-lam[loc] for loc in lam_locants)
    suffix_locants = sorted(position_of[a] for a, xs in exo_chalcogens.items() for x in xs if x[1] == suffix_kind)
    return hetero_locants, seniority, lam_locants, lam_numbers, suffix_locants


def name_ring_lambda_heterone(mol) -> str:
    ring_order, ring_hetero, exo_chalcogens_raw, roots = _match(mol)
    if specified_stereocenters(mol) is not None or any(
        atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()
    ):
        raise UnsupportedStructure("stereodescriptors on a ring λ-heterone are not supported yet")

    exo_chalcogens = {
        atom: [(x, _RING_HETERO[mol.GetAtomWithIdx(x).GetAtomicNum()]) for x in xs]
        for atom, xs in exo_chalcogens_raw.items()
    }
    kinds = {kind for xs in exo_chalcogens.values() for _, kind in xs}
    suffix_kind = next(e for e in _SENIORITY if e in kinds)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    branch_names = {
        (atom, root): name_branch(graph, root, atom, halogens, mol=mol)
        for atom, rts in roots.items()
        for root in rts
    }

    size = len(ring_order)
    elements_in_order = [ring_hetero[a] for a in ring_order if a in ring_hetero]
    mancude = any(b.GetBondTypeAsDouble() == 2.0 and b.IsInRing() for b in mol.GetBonds())
    stem = _stem(elements_in_order, size, mancude)

    best = None
    for base in (ring_order, list(reversed(ring_order))):
        for start in range(size):
            candidate = base[start:] + base[:start]
            position_of = {a: i + 1 for i, a in enumerate(candidate)}
            hetero_key = _numbering_key(position_of, ring_hetero, exo_chalcogens, suffix_kind)
            prefixes = {}
            for atom, xs in exo_chalcogens.items():
                for _, kind in xs:
                    if kind != suffix_kind:
                        info = prefixes.setdefault(_PREFIX[kind], {"locants": [], "compound": False})
                        info["locants"].append(position_of[atom])
            for (atom, _root), (name, compound) in branch_names.items():
                info = prefixes.setdefault(name, {"locants": [], "compound": compound})
                info["locants"].append(position_of[atom])
            locant_set, _, citation = substituent_locant_set_and_citation(prefixes)
            key = (hetero_key, locant_set, citation)
            if best is None or key < best[0]:
                best = (key, position_of, prefixes)

    _, position_of, prefixes = best
    lam_by_position = {
        position_of[a]: _lambda(ring_hetero[a], exo_chalcogens[a]) for a in ring_hetero if a in exo_chalcogens
    }
    cited = []
    for element in _SENIORITY:
        for a in sorted((x for x, e in ring_hetero.items() if e == element), key=position_of.get):
            loc = position_of[a]
            cited.append(f"{loc}λ{lam_by_position[loc]}" if loc in lam_by_position else str(loc))
    suffix_locants = sorted(
        position_of[a] for a, xs in exo_chalcogens.items() for _, kind in xs if kind == suffix_kind
    )

    word = multiplied_word(len(suffix_locants), _SUFFIX[suffix_kind])
    base = stem[:-1] if word[0] in "aeiouy" else stem
    indicated = f"{','.join(f'{loc}H' for loc in sorted(lam_by_position))}-" if mancude else ""
    name = f"{indicated}{','.join(cited)}-{base}-{','.join(map(str, sorted(suffix_locants)))}-{word}"
    prefix = format_substituent_prefixes(prefixes)
    return f"{prefix}-{name}" if prefix else name
