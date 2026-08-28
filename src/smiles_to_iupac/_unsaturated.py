"""Naming of acyclic hydrocarbons containing one or more carbon-carbon double
and/or triple bonds (alkenes, alkynes, and mixed enynes) on the principal
chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-31.1.1.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): the
  presence of a double or triple bond in an otherwise saturated chain is
  denoted by changing the 'ane' ending to 'ene' or 'yne'. Locants as low as
  possible are given to multiple bonds as a set, even though this may at
  times give 'yne' endings lower locants than 'ene' endings; if a choice
  remains, preference for low locants is given to the double bonds. In
  names, 'ene' always precedes 'yne', with elision of the final letter 'e'
  in 'ene' when directly followed by 'yne' (e.g. 'pent-3-en-1-yne (PIN)',
  'pent-1-en-4-yne (PIN)').
- P-31.1.1.2 (Chapter P-3): the multiplying prefixes 'di', 'tri', etc. are
  placed before 'ene'/'yne' to indicate the number of multiple bonds of each
  kind. For euphonic reasons, when 'ene'/'yne' are preceded by a multiplying
  prefix and a locant, the letter 'a' is inserted after the parent stem,
  e.g. 'buta-1,3-diene (PIN)', 'nona-1,3,5,7-tetraene (PIN)'. There is no
  elision of a multiplying prefix's own final 'a' before 'ene'/'yne'
  ('tetraene', not 'tetrene').
- P-31.1.2.2.1 (Chapter P-3): systematic acyclic examples combining the two
  rules above, e.g. 'hexa-2,3-diene (PIN)', 'hexa-1,3-dien-5-yne (PIN)'.
- P-31.1.2.1 (Chapter P-3): 'acetylene' is the retained preferred IUPAC name
  for unsubstituted HC#CH (substitution of any kind is disallowed for this
  retained name); once substituted, chain-extended, or combined with a
  second multiple bond, the systematic '-yne' name is used instead.
- P-14.3.4.2(d) / P-14.3.3 (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf):
  the locant '1' is omitted only for unsubstituted, two- or three-carbon
  alkenes/alkynes carrying a single multiple bond ('ethene (PIN)',
  'acetylene (PIN)', 'propene (PIN)', 'propyne (PIN)'); with a substituent,
  a longer chain, or a second multiple bond, the chain is always long enough
  that more than one position is possible, so the locant(s) are always
  essential and this omission never applies.
- P-44.3.2 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf):
  chain length is chosen first (as in `_acyclic.py`), i.e. the principal
  chain has the greatest number of skeletal atoms.
- P-44.4.1.1 / P-44.4.1.2 (Chapter P-4): among chains tied for length, the
  principal chain has the greater number of multiple bonds, then (if still
  tied) the greater number of double bonds. This module only supports the
  case where some candidate longest chain carries *every* multiple bond in
  the molecule (any multiple bond left off the principal chain would need an
  alkenyl/alkynyl substituent prefix, which is out of scope, see below), so
  in practice this criterion reduces to: only longest chains containing all
  multiple bonds are eligible as the principal chain.
- P-14.4(e) / P-44.4.1.10, P-44.4.1.10.1 (Chapter P-1 / P-4): numbering
  direction is chosen to give the lowest locants to the full set of multiple
  bonds (ene and yne together) ahead of substituent locants; if a choice
  still remains, lower locants go to the double ('ene') bonds specifically
  (same rule restated for the naming/citation format in P-31.1.1.1 above).
- P-45.2 (Chapter P-4) / P-14.4(f,g) (Chapter P-1): any remaining tie is
  broken exactly as in `_acyclic.py`: by substituent count, then lowest
  locant set, then lowest locants in citation order.
- P-29.4 / P-46 (Chapter P-2, P-4): substituents on the chain are named the
  same way as for alkanes/cycloalkanes, via `_substituents.py`, so simple
  *and* branched ("compound") substituents (e.g. isopropyl) are supported
  here too.
- P-35.2.1 (Chapter P-3): halogen substituents (fluoro, chloro, bromo, iodo)
  are never skeletal atoms (P-44.3), so the principal chain is found over
  carbon-carbon connectivity only (`carbon_adjacency`, see `_common.py`)
  while substituent detection still uses the full atom graph.

A multiple bond located in a substituent rather than the principal chain
(i.e. no candidate longest chain carries every multiple bond in the
molecule), and unsaturation in a ring, are out of scope and raise
`UnsupportedStructure`.

- P-93 (Chapter P-9, https://iupac.qmul.ac.uk/BlueBook/P9.html), as of
  `tasks/ez-double-bond-naming.md` (2026-08-25): when the molecule has
  exactly one C=C double bond in total (no triple bond, no second double
  bond) and its geometry is specified in the input (`/`/`\`), a
  "(E)-"/"(Z)-" prefix is added to the whole name, e.g. "(E)-but-2-ene",
  "(Z)-2-chlorobut-2-ene" (both cross-checked against PubChem). No locant
  is included in the prefix (R-7.1.2/P-93: a locanted "(2E)-" form is only
  needed when there's more than one stereogenic double bond to
  distinguish, which never arises in this single-double-bond scope). CIP
  priority computation is delegated entirely to RDKit
  (`_common.single_specified_double_bond_stereo`), mirroring
  `_alcohol.py`'s R/S handling: a non-stereogenic double bond, or one left
  unspecified in the input, is not a new rejection case -- it's named
  exactly as before (no prefix). A specified double bond alongside a
  triple bond, a second double bond, or a tetrahedral stereocenter is out
  of scope and raises `UnsupportedStructure` explicitly.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    path_between,
    single_specified_double_bond_stereo,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_VALID_ORDERS = (_ENE_ORDER, _YNE_ORDER)


def _longest_chains(graph):
    """All maximum-length simple paths in the tree (P-44.3.2: greater number
    of skeletal atoms)."""
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


def _substituents_for_chain(graph, chain, halogens):
    """Return {position (1-based) -> [(name, is_compound), ...]} for a
    candidate chain."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _bond_locant(chain, bond_atoms):
    """1-based locant of the lower-numbered atom of a multiple bond under
    this chain ordering, or None if the bond isn't an edge of this chain."""
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def _bond_locants(chain, bonds):
    """(ene_locants, yne_locants) for every bond under this chain ordering,
    or None if any bond isn't an edge of this chain."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == _ENE_ORDER else yne).append(locant)
    return ene, yne


def _multiplied_word(count, base):
    """P-31.1.1.2: 'ene'/'yne' with a multiplying prefix ('di', 'tri', ...)
    for two or more bonds of the same kind; bare 'ene'/'yne' for exactly
    one; '' for none. P-16.3.3: a multiplying prefix's terminal 'a' is
    elided before a suffix beginning with 'a' or 'o' (see `_common.py`'s
    `multiplied_word` docstring for the confirmed examples this mirrors)."""
    if count == 0:
        return ""
    if count == 1:
        return base
    prefix = numerical_term(count)
    if prefix.endswith("a") and base[:1] in "ao":
        prefix = prefix[:-1]
    return prefix + base


def _unsaturation_suffix(ene_locants, yne_locants):
    """Locant-and-suffix string (e.g. '1,3-dien-5-yne') plus whether the
    parent stem needs its euphonic trailing 'a' (P-31.1.1.2), for a chain's
    full set of multiple bonds. 'ene' is always cited before 'yne'
    (P-31.1.1.1), with its final 'e' elided only when directly followed by
    an unprefixed 'yne' (P-31.1.1.1's own examples; a multiplying-prefixed
    'diyne'/'triyne' begins with a consonant and elides nothing)."""
    ene_locants = sorted(ene_locants)
    yne_locants = sorted(yne_locants)
    ene_count, yne_count = len(ene_locants), len(yne_locants)
    ene_word = _multiplied_word(ene_count, "ene")
    yne_word = _multiplied_word(yne_count, "yne")

    if ene_count and yne_count:
        ene_part = ene_word[:-1] if yne_word[0] in "aeiouy" else ene_word
        ene_loc_str = ",".join(str(loc) for loc in ene_locants)
        yne_loc_str = ",".join(str(loc) for loc in yne_locants)
        body = f"{ene_loc_str}-{ene_part}-{yne_loc_str}-{yne_word}"
    elif ene_count:
        body = f"{','.join(str(loc) for loc in ene_locants)}-{ene_word}"
    else:
        body = f"{','.join(str(loc) for loc in yne_locants)}-{yne_word}"

    needs_stem_a = (ene_count >= 2) if ene_count else (yne_count >= 2)
    return body, needs_stem_a


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-3]
    single_bond = len(ene_locants) + len(yne_locants) == 1
    if single_bond and chain_length <= 3 and not prefix:
        # P-14.3.4.2(d): the locant is omittable only when unsubstituted, a
        # single multiple bond, and the chain is short enough that no other
        # position is possible.
        if chain_length == 2 and yne_locants:
            return "acetylene"
        return stem + ("ene" if ene_locants else "yne")
    if single_bond and chain_length == 2 and len(grouped) == 1:
        name = next(iter(grouped))
        info = grouped[name]
        if len(info["locants"]) == 1 and not info["compound"]:
            # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly
            # one substituent has only one possible structure regardless of
            # numbering direction, so both the multiple bond's and the
            # substituent's locants are omittable, e.g. 'fluoroethyne (PIN)'
            # for fluoroacetylene (P-31.1.2.1).
            return name + stem + ("ene" if ene_locants else "yne")
    body, needs_stem_a = _unsaturation_suffix(ene_locants, yne_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    """Sort key implementing P-14.4(e)/P-44.4.1.10 then P-45.2.1-P-45.2.3,
    most-preferred first."""
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
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
    # The multiple bonds' locants (as a set, then 'ene' locants specifically)
    # outrank every substituent criterion; substituent count and lower
    # locants are preferred, so negate the count to sort every field in
    # ascending "most preferred first" order.
    return (
        (combined_locant_set, ene_locant_set, -total_count, locant_set, citation_locants, name),
        name,
    )


def name_acyclic_unsaturated(mol) -> str:
    validate_atoms_and_bonds(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported by this module (see smiles_to_iupac._cyclic)"
        )

    bonds = non_single_bonds(mol)
    if not bonds or any(order not in _VALID_ORDERS for _, _, order in bonds):
        raise UnsupportedStructure(
            "this module handles only carbon-carbon double and triple bonds; "
            "a bond order other than double or triple is not supported (see "
            "P-31.1.1.1)"
        )

    stereo = single_specified_double_bond_stereo(mol)
    if stereo is not None and len(bonds) != 1:
        raise UnsupportedStructure(
            "a specified double-bond E/Z stereo element combined with any "
            "other multiple bond (a second double bond or a triple bond) "
            "is not supported yet (see P-93)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    # P-44.3.2 / P-44.4.1.1: among the longest chains, only those containing
    # every multiple bond in the molecule can be the principal chain (a bond
    # left off the chain would need an alkenyl/alkynyl substituent prefix,
    # out of scope here).
    chains_with_all_bonds = [c for c in chains if _bond_locants(c, bonds) is not None]
    if not chains_with_all_bonds:
        raise UnsupportedStructure(
            "not every multiple bond lies on a single longest chain; "
            "expressing one in a substituent (an alkenyl/alkynyl prefix) is "
            "not supported yet (see P-29.2, P-32.1)"
        )

    best_key = None
    best_name = None
    for chain in chains_with_all_bonds:
        for candidate in (chain, list(reversed(chain))):
            ene_locants, yne_locants = _bond_locants(candidate, bonds)
            substituents = _substituents_for_chain(graph, candidate, halogens)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    if stereo is not None:
        # P-93/R-7.1.2: no locant in the prefix -- with only one multiple
        # bond in the whole molecule, there's nothing to disambiguate.
        return f"({stereo[1]})-{best_name}"
    return best_name
