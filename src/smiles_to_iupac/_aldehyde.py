"""Naming of aldehydes (the '-al' suffix, terminal -CHO) on acyclic saturated
or unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-33.3, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'al' is the preselected suffix for -CHO, cited in a combined chain/suffix
  name the same way 'ol'/'amine'/'one' are for -OH/-NH2/=O (`_alcohol.py`/
  `_amine.py`/`_ketone.py`), e.g. 'propanal', 'pentanedial'. Table 3.3 ranks
  'al' senior to 'one'/'ol'/'amine'; a carbonyl carbon shaped like a ketone
  (two carbon neighbors) or a carboxylic acid/ester/amide (a second oxygen on
  the same carbon) is rejected outright rather than silently named as if it
  were an aldehyde, and this module never attempts suffix-vs-suffix
  seniority competition against a coexisting senior group either (carboxylic
  acid/ester/amide/nitrile) since only C, halogen, and aldehyde-shaped
  oxygen atoms are accepted at all.
- Unlike -OH/=O, a -CHO carbon is always a chain terminus (it has exactly one
  carbon neighbor, being otherwise saturated by =O and one H), so it is
  never a genuine locant choice: whichever end of the principal chain bears
  it is numbered C1 (P-44.4.1.8, suffix locants minimized first, same
  ordering as `_alcohol.py`/`_amine.py`/`_ketone.py`), and that locant is
  always '1' (or '1' and the last position, for a chain with -CHO at both
  ends) - so unlike those modules, this one never actually cites an 'al'
  locant in the name (P-14.3.3's general omission-when-unambiguous applies
  uniformly here, not just to a mononuclear/two-carbon parent): 'propanal',
  not 'propan-1-al'; 'pentanedial', not 'pentane-1,5-dial'. An 'ene'/'yne'
  locant elsewhere on the chain is still cited as usual, e.g. 'hex-4-enal'.
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-aldehyde name reuses the same mechanics
  as `_unsaturated.py`/`_ketone.py` (ending replaces 'ane' entirely,
  multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before
  'yne', euphonic stem 'a' before a multiplied ending); the 'al'/'dial'
  suffix word is appended directly onto the last such segment (eliding a
  trailing vowel the same way 'ene'+'ol' does in `_alcohol.py`), since it
  never carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'al' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that isn't a doubly-bonded, isolated aldehyde carbonyl oxygen
  (ethers, any singly-bonded -OH-shaped oxygen indicating a coexisting
  alcohol, and any oxygen bonded to more than one heavy atom).
- A carbonyl carbon with other than exactly one carbon neighbor: zero (a
  carbon-less carbonyl, e.g. formaldehyde) or two-or-more (a ketone) -
  neither is this module's territory.
- An aromatic carbonyl carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory.
- Any other heteroatom (N, S, ...).
- -CHO on a ring (the 'carbaldehyde' suffix, P-33.3.1.2, a substituent-style
  name rather than this module's parent-chain suffix) - deferred entirely;
  only an acyclic terminal -CHO is supported here.
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
_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _validate_and_collect_aldehydes(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return the set of aldehyde carbonyl-oxygen atom indices."""
    aldehydes = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an aldehyde carbonyl oxygen (P-33.3) "
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
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope; only an isolated aldehyde "
                    "carbonyl is supported (Table 3.3, P-33.3)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "a singly-bonded oxygen indicates a coexisting alcohol "
                    "(-OH) rather than a plain aldehyde, which this module "
                    "does not attempt to disambiguate (Table 3.3's suffix "
                    "seniority order, P-33.3)"
                )
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("an aldehyde carbonyl must be attached to a carbon atom")
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for this "
                    "module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "a carbonyl carbon with other than exactly one carbon "
                    "neighbor (a ketone, or a carbon-less carbonyl such as "
                    "formaldehyde) is not an aldehyde and is out of scope "
                    "for this module (Table 3.3)"
                )
            aldehydes.add(atom.GetIdx())
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
    if not aldehydes:
        raise UnsupportedStructure("no aldehyde (-CHO) group found; this module only handles aldehydes")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return aldehydes


def _multiplied_word(count, base):
    if count == 0:
        return ""
    if count == 1:
        return base
    return numerical_term(count) + base


def _suffix_body(ene_locants, yne_locants, al_count):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'al' ending
    (e.g. '4-enal'). Unlike `_ketone.py`'s '-one', 'al'/'dial' never carries
    its own locant (see module docstring), so it is glued directly onto the
    preceding word instead of getting a hyphenated locant segment of its
    own."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))

    al_word = _multiplied_word(al_count, "al")
    words = [word for _, word in segments] + [al_word]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    locanted_parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words[:-1])
    ]
    body = "-".join(locanted_parts) + words[-1] if locanted_parts else words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, al_count, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, al_count)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, al_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first. The 'al' locant set still drives
    orientation choice even though it is never printed (see module
    docstring)."""
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    al_locant_set = lowest_locant_set(al_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, len(al_locants), ene_locants, yne_locants, grouped)
    return (
        (
            al_locant_set,
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


def _al_locants(position_of, aldehydes, graph):
    locants = []
    for o in aldehydes:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, aldehydes):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in aldehydes]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_aldehyde(mol, aldehydes, bonds):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _al_locants(position_of, aldehydes, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every aldehyde-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more -CHO groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            al_locants = _al_locants(position_of, aldehydes, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, aldehydes)
            key, name = _candidate_key(chain_length, al_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_aldehyde(mol) -> str:
    aldehydes = _validate_and_collect_aldehydes(mol)
    graph = adjacency(mol)
    # Exclude each C=O carbonyl bond itself: `non_single_bonds` reports it as
    # order 2.0 same as a C=C, but it isn't a chain 'ene' bond (one endpoint
    # is the aldehyde oxygen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in aldehydes and b[1] not in aldehydes]
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "an aldehyde on a ring (the 'carbaldehyde' suffix, P-33.3.1.2) "
            "is out of scope for this module; only an acyclic terminal "
            "-CHO is supported"
        )
    return _name_acyclic_aldehyde(mol, aldehydes, bonds)
