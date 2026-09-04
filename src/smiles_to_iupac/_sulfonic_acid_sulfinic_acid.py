"""Naming of a molecule combining exactly one sulfonic acid (-SO3H) with
one or more sulfinic acid (-S(=O)OH) groups on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41/P-43 (`_seniority.py`): sulfonic acids (Table 4.4 class 4) outrank
  sulfinic acids (class 9), so a coexisting sulfinic acid is demoted to
  the 'sulfino' substituent prefix (P-65.3.1) instead of its own
  '-sulfinic acid' suffix -- e.g. 'OS(=O)CCS(=O)(=O)O' ->
  '2-sulfinoethane-1-sulfonic acid'. This mirrors
  `_sulfonic_acid_thiol.py`, the first module built on `_seniority.py`,
  with the demoted group swapped from thiol/'sulfanyl' to sulfinic
  acid/'sulfino'.
- Otherwise mirrors `_sulfonic_acid.py`'s acyclic-chain path exactly: the
  P-14.3.4.2(a)/(b) locant-omission rules, and the -SO3H locant minimized
  before substituent-prefix locants (P-45.2).

Scope, deliberately narrow (mirrors `_sulfonic_acid_thiol.py`): a single
sulfonic acid plus one or more sulfinic acids, all on one acyclic
saturated chain, with halogen substituents allowed. Explicitly out of
scope (raise `UnsupportedStructure`): any chain unsaturation (ene/yne),
any ring, more than one sulfonic acid, a coexisting hydroxyl/ether/other
heteroatom, a specified stereocenter, and any sulfonic acid/sulfinic acid
not captured by a single longest chain.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    path_between,
)
from ._numerals import alkane_name
from ._seniority import senior_class
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}

assert senior_class("sulfonic_acid", "sulfinic_acid") == "sulfonic_acid"


def _sulfonic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonic acid group: bonded to exactly
    one carbon, two double-bonded (terminal) oxygens, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_sulfonic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 2 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def _sulfinic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfinic acid group: bonded to exactly
    one carbon, one double-bonded (terminal) oxygen, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_sulfinic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_sulfonic_acid_sulfinic_acid_shape(mol) -> bool:
    return bool(_sulfonic_sulfur_atoms(mol)) and bool(_sulfinic_sulfur_atoms(mol))


def _validate_and_collect(mol):
    sulfonic_atoms = _sulfonic_sulfur_atoms(mol)
    if len(sulfonic_atoms) != 1:
        raise UnsupportedStructure(
            "exactly one sulfonic acid is required; this module only "
            "handles a single sulfonic acid combined with one or more "
            "sulfinic acids"
        )
    (sulfonic_sulfur,) = sulfonic_atoms
    (so3h_carbon,) = (n for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    sulfonic_oxygens = {n.GetIdx() for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 8}

    sulfinics = _sulfinic_sulfur_atoms(mol)
    if not sulfinics:
        raise UnsupportedStructure(
            "no sulfinic acid found; this module only handles a sulfonic "
            "acid combined with at least one sulfinic acid (see "
            "_sulfonic_acid.py for a plain sulfonic acid)"
        )
    sulfinic_idxs = {s.GetIdx() for s in sulfinics}
    sulfinic_oxygens = {
        o.GetIdx() for s in sulfinics for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }

    accounted_sulfur_idxs = {sulfonic_sulfur.GetIdx()} | sulfinic_idxs
    accounted_oxygen_idxs = sulfonic_oxygens | sulfinic_oxygens

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the sulfonic/sulfinic acid groups' "
                "own oxygens (P-65.3.1) and halogen substituents "
                "(P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in accounted_oxygen_idxs:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the sulfonic/sulfinic "
                    "acid groups is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == 16:
            if atom.GetIdx() not in accounted_sulfur_idxs:
                raise UnsupportedStructure(
                    "a sulfur atom not shaped like the sulfonic acid group "
                    "or a plain sulfinic acid is out of scope for this module"
                )
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    return sulfonic_sulfur.GetIdx(), so3h_carbon.GetIdx(), sulfinic_idxs


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, so3h_locant, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "sulfonic acid"
    if chain_length == 2 and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -SO3H) in total omits the locant.
        return alkane_name(2) + "sulfonic acid"

    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{alkane_name(chain_length)}-{so3h_locant}-sulfonic acid"


def _candidate_key(chain_length, so3h_locant, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, so3h_locant, grouped)
    return (so3h_locant, locant_set, citation_locants, name), name


def _longest_chains(graph):
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(path_between(parents[u], u, v))
    return chains


def _substituents_for_chain(graph, chain, names, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, names) for root in branch_roots]
    return substituents


def name_sulfonic_acid_sulfinic_acid(mol) -> str:
    sulfonic_sulfur_idx, so3h_carbon, sulfinic_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a sulfonic acid/sulfinic acid combination on a ring is out "
            "of scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    acid_sulfur_idxs = {sulfonic_sulfur_idx} | sulfinic_idxs
    if any(a not in acid_sulfur_idxs and b not in acid_sulfur_idxs for a, b, _ in all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a sulfonic acid/"
            "sulfinic acid combination is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{s: "sulfino" for s in sulfinic_idxs}}
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    sulfinic_carbons = {
        n.GetIdx()
        for s in sulfinic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    eligible = [
        chain
        for chain in chains
        if so3h_carbon in chain and sulfinic_carbons.issubset(set(chain))
    ]
    if not eligible:
        raise UnsupportedStructure(
            "not every sulfonic acid/sulfinic acid-bearing carbon lies on "
            "a single longest carbon chain; a shorter principal chain is "
            "not supported yet"
        )

    excluded = {sulfonic_sulfur_idx}
    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so3h_locant = position_of[so3h_carbon]
            substituents = _substituents_for_chain(graph, candidate, names, excluded)
            key, name = _candidate_key(chain_length, so3h_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
