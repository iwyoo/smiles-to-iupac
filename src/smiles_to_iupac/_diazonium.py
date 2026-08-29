"""Naming of diazonium cations (R-N#N+), per the IUPAC 2013
Recommendations ("the Blue Book"):

- Chapter P-7 (https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf, P-73): a
  diazonium cation is named by appending the suffix 'diazonium' directly
  to the parent hydride name, with no elision of the parent's terminal
  'e' (unlike the vowel-initial '-ol'/'-one' suffixes) -- 'methane' +
  'diazonium' -> 'methanediazonium'. This is structurally the same
  chain-search-and-locant machinery `_sulfonic_acid.py` already uses for
  its own '-SO3H' suffix, just with the diazonium nitrogen attached
  directly to the chain carbon instead of through an intervening sulfur.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_sulfonic_acid.py` apply (mononuclear parent, or a homogeneous
  two-carbon chain with exactly one substituent in total), e.g.
  'ethanediazonium'.
- Confirmed via PubChem structure match: `C[N+]#N` -> "methanediazonium",
  `CC[N+]#N` -> "ethanediazonium", `CCC[N+]#N` -> "propane-1-diazonium",
  `CC(C)[N+]#N` -> "propane-2-diazonium" (a branched chain, supported the
  same way `_sulfonic_acid.py` supports one), `ClCC[N+]#N` ->
  "2-chloroethanediazonium" (halogen coexistence).

Scope, deliberately narrow, mirroring `_sulfonic_acid.py`'s own chain
scope (acyclic only for this first pass -- no monocyclic diazonium shape
has been verified here): a single -N#N+ on an acyclic chain (branched,
unbranched, or unsaturated), with no other heteroatom anywhere in the
molecule except the diazonium group's own two nitrogens. Explicitly out
of scope (raise `UnsupportedStructure`): any ring, two or more diazonium
groups, and a diazonium group on a carbon that is also part of a C=C/C#C
bond.
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
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _diazonium_nitrogens(mol):
    """Terminal, +1-charged nitrogens shaped like a diazonium group: triple
    bonded to a second, neutral, degree-1 nitrogen, and singly bonded to
    exactly one carbon."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1 or atom.GetDegree() != 2:
            continue
        if atom.GetIsotope() != 0:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(nitrogens) != 1:
            continue
        (carbon,) = carbons
        (terminal_n,) = nitrogens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), terminal_n.GetIdx()).GetBondTypeAsDouble() != _YNE_ORDER:
            continue
        if terminal_n.GetDegree() != 1 or terminal_n.GetFormalCharge() != 0 or terminal_n.GetIsotope() != 0:
            continue
        matches.append(atom)
    return matches


def has_diazonium_shape(mol) -> bool:
    return bool(_diazonium_nitrogens(mol))


def _validate_and_collect_diazonium(mol):
    nitrogens = _diazonium_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure(
            "no diazonium (-N#N+) group found; this module only handles "
            "diazonium cations"
        )
    if len(nitrogens) > 1:
        raise UnsupportedStructure(
            "more than one diazonium group is out of scope for this module"
        )
    (nitrogen,) = nitrogens
    diazonium_atom_idxs = {nitrogen.GetIdx()}
    (terminal_n,) = (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 7)
    diazonium_atom_idxs.add(terminal_n.GetIdx())

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        elif atom.GetIdx() not in diazonium_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a diazonium group (P-73) and "
                "halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetIdx() not in diazonium_atom_idxs and (atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0):
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported yet (this module only handles "
            "acyclic chains)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (carbon,) = (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 6)
    return carbon.GetIdx(), diazonium_atom_idxs


def _reject_enediazonium_carbon(graph, diazonium_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if diazonium_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a diazonium group on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _multiplied_word(count, base):
    """P-16.3.3: a multiplying prefix's terminal 'a' is elided before a
    suffix beginning with 'a' or 'o' (see `_common.py`'s `multiplied_word`
    docstring for the confirmed examples this mirrors)."""
    if count == 0:
        return ""
    if count == 1:
        return base
    prefix = numerical_term(count)
    if prefix.endswith("a") and base[:1] in "ao":
        prefix = prefix[:-1]
    return prefix + base


def _suffix_body(ene_locants, yne_locants, diazonium_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append(([diazonium_locant], "diazonium"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, diazonium_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "diazonium"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole diazonium group) in total omits the
        # locant, e.g. 'ethanediazonium'.
        return alkane_name(2) + "diazonium"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, diazonium_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, diazonium_locant, ene_locants, yne_locants, substituents):
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, diazonium_locant, ene_locants, yne_locants, grouped)
    return (
        (
            diazonium_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


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


def _bond_locant(chain, bond_atoms):
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def _bond_locants(chain, bonds):
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == _ENE_ORDER else yne).append(locant)
    return ene, yne


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def name_diazonium(mol) -> str:
    diazonium_carbon, excluded = _validate_and_collect_diazonium(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    # The diazonium N#N triple bond is the one non-single bond with both
    # atoms in `excluded`; anything else non-single must be a chain
    # ene/yne bond.
    bonds = [
        b for b in all_non_single
        if b[2] in (_ENE_ORDER, _YNE_ORDER) and b[0] not in excluded and b[1] not in excluded
    ]
    diazonium_bonds = [b for b in all_non_single if b[0] in excluded and b[1] in excluded]
    if len(bonds) + len(diazonium_bonds) != len(all_non_single) or len(diazonium_bonds) != 1:
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enediazonium_carbon(graph, diazonium_carbon, bonds)

    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if diazonium_carbon not in chain:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the diazonium-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            diazonium_locant = position_of[diazonium_carbon]
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, diazonium_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
