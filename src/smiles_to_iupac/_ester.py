"""Naming of esters (the '-oate' suffix, R-COO-R') on acyclic saturated or
unsaturated acyl chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-65.6.3, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-65.6.3; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): an ester R-CO-O-R' is named
  as two words, "R'yl R-oate" — the alcohol part (R') cited first as a plain
  substituent-group name, then the acyl part (R) named the same way a
  carboxylic acid's R-COOH would be but with 'oate' instead of 'oic acid'.
  'oate' ranks junior to 'oic acid' and senior to 'amide' in Table 3.3.
- This module's scope (see also `_carboxylic_acid.py`/`_aldehyde.py`, the
  most similar existing modules): the acyl part (R) may be any acyclic
  saturated or unsaturated chain, with halogen substituents, the same as
  `_carboxylic_acid.py`'s R. The alcohol part (R') is restricted to a plain,
  unbranched, unsubstituted, saturated alkyl group attached at its own chain
  terminus (e.g. 'methyl', 'ethyl', 'propyl') — a branched, substituted,
  unsaturated, or ring-bearing R' is deferred (see "Explicitly out of
  scope" below).
- The acyl carbon is always a chain terminus (its remaining two bonds, after
  the carbonyl and ester oxygens, allow at most one more substituent, which
  must be another chain carbon, or nothing for a formate ester), so it is
  never a genuine locant choice: it is always C1 of the acyl chain, the same
  way `_carboxylic_acid.py`'s -COOH carbon is, and that locant is never
  cited (P-14.3.3).
- Since this module only ever handles a single, isolated ester group (see
  scope below), the 'oate' word is never multiplied ('dioate' etc. is out
  of scope), unlike `_carboxylic_acid.py`'s 'dioic acid'.
- P-31.0 / P-31.1.1.1-.2: construction of the 'ene'/'yne' portion of the
  acyl part's name reuses the same mechanics as `_carboxylic_acid.py`.
- P-35.2.1: halogen substituents on the acyl chain are prefix-only and
  coexist freely with the 'oate' suffix, reusing `halogen_substituents`/
  `format_substituent_prefixes` unchanged.
- P-29.3.2.1: the alcohol part's name is a plain alkyl substituent-group
  name (P-13.2.1's "R'yl" role, not a locanted prefix), built with
  `alkyl_name` directly since it's restricted to an unbranched chain here.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule (a ring-attached ester or lactone uses a
  different naming construction, P-65.6.3.2/P-65.6.3.3; acyclic only, per
  this module's scope).
- More than one ester group (a diester), or any oxygen that isn't part of
  the single ester's carbonyl/ester-oxygen pair (an ether, alcohol, or
  second carbonyl elsewhere).
- A branched, substituted, unsaturated, or cyclic alcohol part (R') — only
  a plain unbranched saturated alkyl R' is supported in this first pass.
- Any other heteroatom (N, S, ...).
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    bfs,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    linear_branch,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
)
from ._numerals import alkane_name, alkyl_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def has_ester_shape(mol) -> bool:
    """True if some carbon carries both a doubly-bonded, monovalent oxygen
    and a singly-bonded oxygen that is itself bonded to a second carbon (a
    -C(=O)-O-C- pattern), regardless of whether the rest of the molecule is
    in scope. Used by `core.py` to route ahead of the carboxylic-acid/
    aldehyde/ketone dispatch, since an ester carbon would otherwise look
    carboxylic-acid- or carbonyl-shaped to those modules."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        has_carbonyl = any(
            o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            for o in oxygens
        )
        has_ester_oxygen = any(
            o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
            for o in oxygens
        )
        if has_carbonyl and has_ester_oxygen:
            return True
    return False


def _find_ester_group(mol):
    """Locate the molecule's single ester group and return
    (acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon) atoms, after
    checking the molecule has exactly one such group and no other oxygens
    (see module docstring)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_oxygens = [
            o
            for o in oxygens
            if o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        ]
        if carbonyls and ester_oxygens:
            matches.append((atom, carbonyls, ester_oxygens))

    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one ester group is required; zero or multiple ester "
            "groups (e.g. a diester) are not supported yet (P-65.6.3)"
        )
    acyl_carbon, carbonyls, ester_oxygens = matches[0]
    if len(carbonyls) != 1 or len(ester_oxygens) != 1:
        raise UnsupportedStructure(
            "an ester carbon with more than one carbonyl or ester oxygen "
            "does not match a simple ester group"
        )
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if total_oxygens != 2:
        raise UnsupportedStructure(
            "an oxygen outside the single ester group's carbonyl/ester pair "
            "(e.g. an ether or a hydroxyl) is out of scope for this module"
        )

    acyl_carbon_neighbors = [n for n in acyl_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(acyl_carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "an ester carbon with more than one carbon neighbor besides its "
            "two ester oxygens is not a valid ester group"
        )

    carbonyl_oxygen = carbonyls[0]
    ester_oxygen = ester_oxygens[0]
    alcohol_carbon = next(n for n in ester_oxygen.GetNeighbors() if n.GetIdx() != acyl_carbon.GetIdx())
    return acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon


def _alcohol_component(full_graph, alcohol_carbon_idx, ester_oxygen_idx):
    component = set()
    stack = [alcohol_carbon_idx]
    while stack:
        node = stack.pop()
        if node in component:
            continue
        component.add(node)
        for neighbor in full_graph[node]:
            if neighbor != ester_oxygen_idx and neighbor not in component:
                stack.append(neighbor)
    return component


def _name_alcohol_part(mol, alcohol_carbon, ester_oxygen_idx):
    full_graph = adjacency(mol)
    component = _alcohol_component(full_graph, alcohol_carbon.GetIdx(), ester_oxygen_idx)
    for idx in component:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "the alcohol part (R') must be a plain alkyl group; "
                "heteroatoms/halogens there are not supported yet (P-65.6.3)"
            )
    for a, b, _ in non_single_bonds(mol):
        if a in component and b in component:
            raise UnsupportedStructure(
                "unsaturation in the alcohol part (R') is not supported yet "
                "(P-65.6.3)"
            )

    carbon_graph = carbon_adjacency(mol)
    length = linear_branch(carbon_graph, alcohol_carbon.GetIdx(), None)
    if length is None:
        raise UnsupportedStructure(
            "a branched alcohol part (R') is not supported yet (P-65.6.3)"
        )
    return alkyl_name(length)


def _suffix_body(ene_locants, yne_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'oate' ending
    (e.g. '2-enoate'). The ester group's own count is always 1 (see module
    docstring), so 'oate' is never multiplied, unlike
    `_carboxylic_acid.py`'s 'oic'/'dioic'."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["oate"]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    if segments:
        locant_parts = [
            f"{','.join(str(loc) for loc in locants)}-{word}"
            for (locants, _), word in zip(segments, words[:-1])
        ]
        body = "-".join(locant_parts) + words[-1]
    else:
        body = words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if (ene_locants or yne_locants) else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
    return (
        (
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _component_subgraph(graph, start):
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def _substituents_for_chain(graph, chain, halogens, excluded_oxygens):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_oxygens]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyl_part(mol, acyl_carbon, carbonyl_oxygen_idx, ester_oxygen_idx):
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    halogens = halogen_substituents(mol)
    acyl_carbon_idx = acyl_carbon.GetIdx()
    excluded_oxygens = {carbonyl_oxygen_idx, ester_oxygen_idx}

    acyl_graph = _component_subgraph(carbon_graph, acyl_carbon_idx)
    chains = longest_chains(acyl_graph)
    chain_length = len(chains[0])

    all_non_single = [
        b for b in non_single_bonds(mol) if b[0] not in excluded_oxygens and b[1] not in excluded_oxygens
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    eligible = [chain for chain in chains if not bonds or bond_locants(chain, bonds) is not None]
    if not eligible:
        raise UnsupportedStructure(
            "the acyl chain's unsaturation does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acyl_carbon_idx:
                # The ester carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(full_graph, candidate, halogens, excluded_oxygens)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    if best_name is None:
        raise UnsupportedStructure(
            "the ester's acyl carbon does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )
    return best_name


def name_ester(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-attached ester or lactone uses a different naming "
            "construction (P-65.6.3.2/P-65.6.3.3), out of scope for this "
            "acyclic-only module"
        )
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the ester's own oxygens (P-65.6.3) "
                "and halogen substituents (P-35.2.1) are not supported yet"
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
        elif atomic_num != 8 and atom.GetDegree() != 1:
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

    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = _name_acyl_part(mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx())
    return f"{alcohol_name} {acyl_name}"
