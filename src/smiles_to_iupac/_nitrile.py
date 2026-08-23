"""Naming of nitriles (the '-nitrile' suffix, terminal -C#N) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-66.5, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  'nitrile' is the preselected suffix for -C#N, cited in a combined
  chain/suffix name the same way 'al'/'one'/'ol'/'amine' are for their own
  groups (`_aldehyde.py`/`_ketone.py`/`_alcohol.py`/`_amine.py`), e.g.
  'propanenitrile', 'pent-4-enenitrile'. Table 3.3 ranks 'nitrile' junior to
  'amide' and senior to 'al'; this module never attempts suffix-vs-suffix
  seniority competition against a coexisting senior or junior group, since
  only C, N (the nitrile nitrogen itself), and halogen atoms are accepted at
  all (any oxygen routes elsewhere in `core.py`, or is out of scope here).
- Unlike 'al'/'one'/'amine', 'nitrile' begins with a consonant, so the
  euphonic elision `_aldehyde.py` applies before a vowel-initial suffix
  never triggers here: the stem's final 'e' (from the parent alkane/alkene/
  alkyne name) is always kept, e.g. 'ethanenitrile', not 'ethannitrile';
  'pent-4-enenitrile', not 'pent-4-ennitrile'.
- Like -CHO (`_aldehyde.py`), a -C#N carbon is always a chain terminus (it
  has exactly one carbon neighbor, being otherwise saturated by the triple
  bond to N), so it is never a genuine locant choice: whichever end of the
  principal chain bears it is numbered C1 (P-44.4.1.8, suffix locants
  minimized first), and that locant is always '1' and never cited (P-14.3.3):
  'propanenitrile', not 'propane-1-nitrile'.
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-nitrile name reuses the same mechanics as
  `_unsaturated.py`/`_aldehyde.py` (ending replaces 'ane' entirely,
  multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before
  'yne', euphonic stem 'a' before a multiplied ending); the 'nitrile' suffix
  word is appended directly onto the last such segment, since it never
  carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'nitrile' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one nitrile group (a dinitrile, P-66.6.3) - future work.
- Any oxygen - routed to a different module by `core.py`, or out of scope
  entirely if this module is called directly on one.
- A nitrile nitrogen that isn't a plain, isolated -C#N (any degree other
  than 1, or a bond order other than 3.0 to its one carbon neighbor).
- A nitrile carbon with other than exactly one carbon neighbor: zero (bare
  HC#N) or two-or-more (not a valid nitrile shape) - neither is a terminal
  substitutive -C#N.
- An aromatic nitrile carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory.
- -C#N on a ring (the 'carbonitrile' suffix, P-66.5.1.2, a substituent-style
  name rather than this module's parent-chain suffix) - deferred entirely;
  only an acyclic terminal -C#N is supported here.
- Any other heteroatom (O, S, ...), or a nitrile carbon entangled with
  another heteroatom.
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
_NITRILE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def has_nitrile_shape(mol) -> bool:
    """True iff `mol` has at least one nitrogen atom of degree 1, triple-
    bonded to a carbon (a plain -C#N pattern), regardless of whether the
    rest of the molecule is in scope. Used by `core.py` to route ahead of
    the amine dispatch within the no-oxygen nitrogen branch, since a
    nitrile nitrogen would otherwise look like an unhandled shape to
    `name_amine`'s own bond-order check."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        (neighbor,) = atom.GetNeighbors()
        if bond.GetBondTypeAsDouble() == _NITRILE_ORDER and neighbor.GetAtomicNum() == 6:
            return True
    return False


def _validate_and_collect_nitriles(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return the set of nitrile-nitrogen atom indices."""
    nitriles = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a nitrile nitrogen (P-66.5) and "
                "halogen substituents (P-35.2.1) are not supported yet"
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
        elif atomic_num == 7:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a nitrogen bonded to more than one heavy atom is not a "
                    "plain nitrile (-C#N) and is out of scope for this "
                    "module (P-66.5)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a nitrile nitrogen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != _NITRILE_ORDER:
                raise UnsupportedStructure(
                    "a nitrogen not triple-bonded to its one carbon neighbor "
                    "is not a nitrile and is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a nitrile carbon on an aromatic ring is out of scope "
                    "for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "a nitrile carbon with other than exactly one carbon "
                    "neighbor is out of scope for this module (Table 3.3)"
                )
            nitriles.add(atom.GetIdx())
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
    if not nitriles:
        raise UnsupportedStructure("no nitrile (-C#N) group found; this module only handles nitriles")
    if len(nitriles) > 1:
        raise UnsupportedStructure("more than one nitrile group (a dinitrile) is not supported yet")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return nitriles


def _multiplied_word(count, base):
    if count == 0:
        return ""
    if count == 1:
        return base
    return numerical_term(count) + base


def _suffix_body(ene_locants, yne_locants, nitrile_count):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'nitrile'
    ending (e.g. '4-enenitrile'). Like `_aldehyde.py`'s '-al', 'nitrile'
    never carries its own locant (see module docstring), so it is glued
    directly onto the preceding word instead of getting a hyphenated locant
    segment of its own; unlike '-al', it never elides a preceding stem's
    trailing 'e' either, since 'nitrile' starts with a consonant."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))

    nitrile_word = _multiplied_word(nitrile_count, "nitrile")
    words = [word for _, word in segments] + [nitrile_word]

    locanted_parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words[:-1])
    ]
    body = "-".join(locanted_parts) + words[-1] if locanted_parts else words[-1]
    return body


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, nitrile_count, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, nitrile_count)
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, nitrile_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first. The nitrile locant set still drives
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
    nitrile_locant_set = lowest_locant_set(nitrile_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, len(nitrile_locants), ene_locants, yne_locants, grouped)
    return (
        (
            nitrile_locant_set,
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


def _nitrile_locants(position_of, nitriles, graph):
    locants = []
    for n in nitriles:
        (carbon,) = graph[n]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, nitriles):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in nitriles]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_nitrile(mol, nitriles, bonds):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _nitrile_locants(position_of, nitriles, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the nitrile-bearing carbon (and/or multiple bond) does not lie "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing the -C#N group (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            nitrile_locants = _nitrile_locants(position_of, nitriles, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, nitriles)
            key, name = _candidate_key(chain_length, nitrile_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_nitrile(mol) -> str:
    nitriles = _validate_and_collect_nitriles(mol)
    graph = adjacency(mol)
    # Exclude each C#N nitrile bond itself: `non_single_bonds` reports it as
    # order 3.0 same as a C#C, but it isn't a chain 'yne' bond (one endpoint
    # is the nitrile nitrogen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in nitriles and b[1] not in nitriles]
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a nitrile on a ring (the 'carbonitrile' suffix, P-66.5.1.2) is "
            "out of scope for this module; only an acyclic terminal -C#N is "
            "supported"
        )
    return _name_acyclic_nitrile(mol, nitriles, bonds)
